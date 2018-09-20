#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

void FakeMatter_AddMatter(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  for (int k=0; k<cctk_lsh[2]; ++k) {
    for (int j=0; j<cctk_lsh[1]; ++j) {
      for (int i=0; i<cctk_lsh[0]; ++i) {
        int ijk = CCTK_GFINDEX3D(cctkGH, i, j, k);
        if (smask[ijk] > 0) {
          eTtt[ijk] += T00;
        }
      }
    }
  }
}
