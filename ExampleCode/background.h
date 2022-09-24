#ifndef ACDM_BACKGROUND
#define ACDM_BACKGROUND

namespace ACDM
{
//========================================================================================
// Class Background
//========================================================================================
class Background
{
  public:
    double H0_in_km_over_s_Mpc;
    double Om0;
    
    double E2(const double z) const;
    
  private:
};

} //namespace ACDM

#endif
