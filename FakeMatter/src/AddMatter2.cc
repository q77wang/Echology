#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <cassert>
#include <cmath>

using namespace std;

// Matrix determinant
static CCTK_REAL determinant(const CCTK_REAL (&g)[3][3]) {
  return g[0][0] * (g[1][1] * g[2][2] - g[1][2] * g[2][1]) -
         g[1][0] * (g[0][1] * g[2][2] - g[0][2] * g[2][0]) +
         g[2][0] * (g[0][1] * g[1][2] - g[0][2] * g[1][0]);
}

// Matrix inverse
static void inverse(const CCTK_REAL (&g)[3][3], const CCTK_REAL detg,
                    CCTK_REAL (&gu)[3][3]) {
  gu[0][0] = (g[1][1] * g[2][2] - g[1][2] * g[2][1]) / detg;
  gu[1][0] = (g[0][1] * g[2][2] - g[0][2] * g[2][1]) / detg;
  gu[2][0] = (g[0][1] * g[1][2] - g[0][2] * g[1][1]) / detg;
  gu[1][1] = (g[0][0] * g[2][2] - g[0][2] * g[2][0]) / detg;
  gu[2][1] = (g[0][0] * g[1][2] - g[0][2] * g[1][0]) / detg;
  gu[2][2] = (g[0][0] * g[1][1] - g[0][1] * g[1][0]) / detg;
  gu[0][1] = gu[1][0];
  gu[0][2] = gu[2][0];
  gu[1][2] = gu[2][1];
}

// Length of a vector
static CCTK_REAL length(const CCTK_REAL (&g)[3][3], const CCTK_REAL (&x)[3]) {
  CCTK_REAL len = 0;
  for (int a = 0; a < 3; ++a)
    for (int b = 0; b < 3; ++b)
      len += g[a][b] * x[a] * x[b];
  len = sqrt(len);
  return len;
}

// Normalize a vector (to unit length)
static void normalize(CCTK_REAL (&x)[3], CCTK_REAL len) {
  for (int a = 0; a < 3; ++a)
    x[a] /= len;
}

extern "C" void FakeMatter_AddMatter2(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  // Horizon surface
  const int hi = horizon_index;
  if (!sf_active[hi]) {
    CCTK_WARN(CCTK_WARN_ALERT, "No horizon found at this time");
    return;
  }
  // Origin (centre) of horizon
  const CCTK_REAL hc[3] = {sf_origin_x[hi], sf_origin_y[hi], sf_origin_z[hi]};

  // Constant expansion surface
  const int si = surface_index;
  if (!sf_active[si]) {
    CCTK_WARN(CCTK_WARN_ALERT,
              "No constant expansion surface found at this time");
    return;
  }
  // Origin (centre) of constant expansion surface
  const CCTK_REAL sc[3] = {sf_origin_x[si], sf_origin_y[si], sf_origin_z[si]};

  // Loop over all grid points
  for (int k = 0; k < cctk_lsh[2]; ++k) {
    for (int j = 0; j < cctk_lsh[1]; ++j) {
      for (int i = 0; i < cctk_lsh[0]; ++i) {
        int ijk = CCTK_GFINDEX3D(cctkGH, i, j, k);
        // Are we inside the constant expansion surface?
        if (smask[ijk] < 1) {

          // Current point
          const CCTK_REAL xx[3] = {x[ijk], y[ijk], z[ijk]};

          // Projection from 3d onto surface:
          //   theta in [0; pi]
          //   phi in [0; 2pi]

          //   x = r sin(theta) cos(phi)
          //   y = r sin(theta) sin(phi)
          //   z = r cos(theta)

          //   rho^2 = x^2 + y^2
          //   rho / z = tan(theta)
          //   theta = atan (rho / z) = atan2(z, rho)

          //   x / y = tan(phi)
          //   phi = atan(x / y) = atan2(y, x)

          // Project onto horizon

          // TODO: The constant expansion has a different origin --
          // take this into account
          const CCTK_REAL dx[3] = {xx[0] - hc[0], xx[1] - hc[1], xx[2] - hc[2]};

          const CCTK_REAL radius =
              sqrt(pow(dx[0], 2) + pow(dx[1], 2) + pow(dx[2], 2));
          const CCTK_REAL rho = hypot(dx[0], dx[1]);
          const CCTK_REAL theta = atan2(rho, -dx[2]);
          assert(theta >= 0 && theta <= M_PI);
          const CCTK_REAL phi = fmod(atan2(dx[1], dx[0]) + 2 * M_PI, 2 * M_PI);
          assert(phi >= 0 && phi < 2 * M_PI);

          // Find radius
          // Find nearest grid point on surface
          // TODO: Use (at least) linear interpolation instead
          const int itheta =
              lrint((theta - sf_origin_theta[hi]) / sf_delta_theta[hi]);
          const int jphi = lrint((phi - sf_origin_phi[hi]) / sf_delta_phi[hi]);
          assert(itheta >= 0 && itheta < sf_ntheta[hi]);
          assert(jphi >= 0 && jphi < sf_nphi[hi]);

          // Horizon radius
          const int hij = itheta + maxntheta * (jphi + maxnphi * hi);
          const CCTK_REAL horizon_radius = sf_radius[hij];

          // Constant expansion surface radius
          const int sij = itheta + maxntheta * (jphi + maxnphi * si);
          const CCTK_REAL surface_radius = sf_radius[sij];

          // Approximate expansion at current point
          //   Theta(rH) = 0
          //   Theta(rS) = Theta
          // Linear approximation:
          //   Theta(r) = m*r+b
          //     m*rH+b = 0
          //     m*rS+b = Theta
          //     m*(rS-rH) = Theta
          //     m = Theta / (rS - rH)
          //     b = - m * rH
          const CCTK_REAL m =
              surface_expansion / (surface_radius - horizon_radius);
          const CCTK_REAL b = -m * horizon_radius;

          const CCTK_REAL expansion = m * radius + b;

          // Get metric
          const CCTK_REAL g[3][3] = {{gxx[ijk], gxy[ijk], gxz[ijk]},
                                     {gxy[ijk], gyy[ijk], gyz[ijk]},
                                     {gxz[ijk], gyz[ijk], gzz[ijk]}};
          // Metric determinant
          const CCTK_REAL detg = determinant(g);
          // Inverse metric
          CCTK_REAL gu[3][3];
          inverse(g, detg, gu);

          // Find spacelike vector er^a orthoginal to the surface:
          CCTK_REAL er[3] = {dx[0], dx[1], dx[2]};
          CCTK_REAL len = length(gu, er);
          if (len >= 1.0e-10)
            normalize(er, len);
          else
            for (int d = 0; d < 3; ++d)
              er[d] = 0.0;

          // Calculate surface two-metric q_ab
          //   qu^ab = gu^ab - er^a er^b
          CCTK_REAL qu[3][3];
          for (int a = 0; a < 3; ++a)
            for (int b = 0; b < 3; ++b)
              qu[a][b] = gu[a][b] - er[a] * er[b];
          CCTK_REAL q[3][3];
          for (int a = 0; a < 3; ++a)
            for (int b = 0; b < 3; ++b) {
              q[a][b] = 0;
              for (int c = 0; c < 3; ++c)
                for (int d = 0; d < 3; ++d)
                  q[a][b] += g[a][c] * g[b][d] * qu[c][d];
            }

          // Calculate fake pressure
          const CCTK_REAL press =
              param_a * (1 + 2 * param_a) /
              (64 * M_PI * pow(param_M, 2) *
               pow(param_a + pow(param_M, 2) * pow(expansion, 2), 2));

          // Define fake T_ab
          CCTK_REAL T[3][3];
          for (int a = 0; a < 3; ++a)
            for (int b = 0; b < 3; ++b)
              T[a][b] = 1.0 / 3.0 * press * q[a][b];

          // // Ensure that the constraints hold
          // // McLachlan says:
          // //   rho   ~ T00 - 2 beta[ua] T0[la] + beta[ua] beta[ub] T[la,lb]
          // //   S[la] ~ T0[la] - beta[ub] T[la,lb]
          // // If rho and S_a are zero, then the constraints will hold.
          //
          // const CCTK_REAL beta[3] = {betax[ijk], betay[ijk], betaz[ijk]};
          //
          // CCTK_REAL Tt[3];
          // for (int a = 0; a < 3; ++a) {
          //   Tt[a] = 0;
          //   for (int b = 0; b < 3; ++b)
          //     Tt[a] += beta[b] * T[a][b];
          // }
          // CCTK_REAL Ttt;
          // Ttt = 0;
          // for (int a = 0; a < 3; ++a)
          //   for (int b = 0; b < 3; ++b)
          //     Ttt += 2 * beta[a] * Tt[a] - beta[a] * beta[b] * T[a][b];

          const CCTK_REAL Ttt = 0;
          const CCTK_REAL Tt[3] = {0, 0, 0};

          eTtt[ijk] += Ttt;
          eTtx[ijk] += Tt[0];
          eTty[ijk] += Tt[1];
          eTtz[ijk] += Tt[2];
          eTxx[ijk] += T[0][0];
          eTxy[ijk] += T[0][1];
          eTxz[ijk] += T[0][2];
          eTyy[ijk] += T[1][1];
          eTyz[ijk] += T[1][2];
          eTzz[ijk] += T[2][2];
        }
      }
    }
  }
}
