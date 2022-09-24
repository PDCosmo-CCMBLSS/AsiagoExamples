#include "background.h"

#include <cmath>

namespace ACDM
{
//========================================================================================
// Class Background
//========================================================================================
double Background::E2(const double z) const
{
    return Om0*pow(1.+z, 3) + (1.-Om0);
};

} //namespace ACDM
