#ifndef GFHypWasaFIELDNEW_H
#define GFHypWasaFIELDNEW_H

#include "AbsBField.h"
#include "TVirtualMagField.h"

namespace genfit {


class GFWasaMap : public AbsBField 
{
 public:
  explicit GFWasaMap(TVirtualMagField* field);
    
  virtual TVector3 get(const TVector3& pos) const;
  virtual void get(const double& posX, const double& posY, const double& posZ, double& Bx, double& By, double& Bz) const ;

private:

  // bool ConstField = false;  
  // double Bz;
  TVirtualMagField* const FieldMap;

};

}
#endif
