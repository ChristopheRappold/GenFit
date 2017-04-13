#include "TMath.h"
#include "GFWasaMap.h"
#include "Riostream.h"

namespace genfit {

// GFWasaMap::GFWasaMap(double F):ConstField(true),Bz(F),FieldMap(nullptr)
// {

// }

GFWasaMap::GFWasaMap(TVirtualMagField* field):FieldMap(field)
{

}

  
TVector3 GFWasaMap::get(const TVector3& Pos) const
{
  TVector3 temp_B;
  Double_t B[3];
  Double_t pos[3];
  
  pos[0]=Pos.X();
  pos[1]=Pos.Y();
  pos[2]=Pos.Z();
  
  FieldMap->Field(pos,B);
  
  temp_B.SetXYZ(B[0],B[1],B[2]);
  return temp_B;
}


void GFWasaMap::get(const double& posX, const double& posY, const double& posZ, double& Bx, double& By, double& Bz) const 
{
  Double_t B_field[3];
  Double_t pos[3]={posX,posY,posZ};
  
  FieldMap->Field(pos,B_field);
  
  Double_t temp_B[3]={B_field[0],B_field[1],B_field[2]};
  
  Bx = temp_B[0];
  By = temp_B[1];
  Bz = temp_B[2];
}

}
