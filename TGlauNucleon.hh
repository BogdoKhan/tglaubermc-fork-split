#pragma once
#include <TVector3.h>
#include <TSystem.h>

class TGlauNucleon : public TObject
{
  protected:
    Double32_t fX;            //Position of nucleon
    Double32_t fY;            //Position of nucleon
    Double32_t fZ;            //Position of nucleon
    Int_t      fType;         //0 = neutron, 1 = proton
    Bool_t     fInNucleusA;   //=1 from nucleus A, =0 from nucleus B
    Int_t      fNColl;        //Number of binary collisions
    Double32_t fEn;           //Energy
  public:
    TGlauNucleon() : fX(0), fY(0), fZ(0), fInNucleusA(0), fNColl(0), fEn(0) {}
    virtual   ~TGlauNucleon() {}
    void       Collide()                                  {++fNColl;}
    Double_t   Get2CWeight(Double_t x) const              {return 2.*(0.5*(1-x)+0.5*x*fNColl);}
    Double_t   GetEnergy()             const              {return fEn;}
    Int_t      GetNColl()              const              {return fNColl;}
    Int_t      GetType()               const              {return fType;}
    Double_t   GetX()                  const              {return fX;}
    Double_t   GetY()                  const              {return fY;}
    Double_t   GetZ()                  const              {return fZ;}
    Bool_t     IsNeutron()             const              {return (fType==0);}
    Bool_t     IsInNucleusA()          const              {return fInNucleusA;}
    Bool_t     IsInNucleusB()          const              {return !fInNucleusA;}
    Bool_t     IsProton()              const              {return (fType==1);}
    Bool_t     IsSpectator()           const              {return !fNColl;}
    Bool_t     IsWounded()             const              {return fNColl>0;}
    void       Reset()                                    {fNColl=0;}
    void       RotateXYZ(Double_t phi, Double_t theta);
    void       RotateXYZ_3D(Double_t psiX, Double_t psiY, Double_t psiZ);
    void       SetEnergy(Double_t en)                     {fEn = en;}
    void       SetInNucleusA()                            {fInNucleusA=1;}
    void       SetInNucleusB()                            {fInNucleusA=0;}
    void       SetNColl(Int_t nc)                         {fNColl = nc;}
    void       SetType(Bool_t b)                          {fType = b;}
    void       SetXYZ(Double_t x, Double_t y, Double_t z) {fX=x; fY=y; fZ=z;}
    ClassDef(TGlauNucleon,4) // TGlauNucleon class
};

