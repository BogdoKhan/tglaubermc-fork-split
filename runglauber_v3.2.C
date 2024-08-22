/*
 $Id: runglauber.C 186 2019-01-13 17:33:43Z loizides $
 -------------------------------------------------------------------------------------
 Latest documentation: https://arxiv.org/abs/1710.07098
 -------------------------------------------------------------------------------------
 To run the code, you need to have the ROOT (http://root.cern.ch/drupal/)
 environment. On the root prompt, then enter
 root [0] gSystem->Load("libMathMore")
 root [1] .L runglauber_X.Y.C+
 (where X.Y denotes the version number).
 If you do not have libMathMore comment out "#define HAVE_MATHMORE" below.
 See the documentation for more information.
 -------------------------------------------------------------------------------------
 v3.2: Incorporates changes from v2.7
 -------------------------------------------------------------------------------------
 v3.1:
  Fixes related to spherical nuclei, as well as consistent set of reweighted profiles 
  for Cu, Au and Xe, see https://arxiv.org/abs/1710.07098v2
 -------------------------------------------------------------------------------------
 v3.0:
  Major update to include separate profile for protons and neutrons, placement of nucleon 
  dof on lattice, as well as reweighted profiles for recentering, 
  see https://arxiv.org/abs/1710.07098v1
 -------------------------------------------------------------------------------------
 v2.7: 
  New macro "runAndOutputLemonTree" for IP-Jazma input (1808.01276), as well as nucleon 
  configurations for He4, C, and O from wavefunction calculations, clarified use of Hulthen
  for deuteron, harmonic oscillator param for O, and new mode to use GlauberGribov also 
  in AA (enable with SetCalcAAGG)
 -------------------------------------------------------------------------------------
 v2.6:
  Includes runAndCalcDens macro, as well as definition for Al, and fixes beta4 for Si2,
  see https://arxiv.org/abs/1408.2549v8
 -------------------------------------------------------------------------------------
 v2.5:
  Include core/corona determination in Npart, and if requested for area from mc and eccentricity,
  as well as various Xe parameterizations including deformation,
  see https://arxiv.org/abs/1408.2549v7
 -------------------------------------------------------------------------------------
 v2.4: 
  Minor update to include Xenon and fix of the TGlauberMC::Draw function, 
  see https://arxiv.org/abs/1408.2549v4
 -------------------------------------------------------------------------------------
 v2.3: 
  Small bugfixes, see https://arxiv.org/abs/1408.2549v3
 -------------------------------------------------------------------------------------
 v2.2:
  Minor update to provide higher harmonic eccentricities up to n=5, and the average
  nucleon--nucleon impact parameter (bNN) in tree output. 
 -------------------------------------------------------------------------------------
 v2.1: 
  Minor update to include more proton pdfs, see https://arxiv.org/abs/1408.2549v2
 -------------------------------------------------------------------------------------
 v2.0: 
  First major update with inclusion of Tritium, Helium-3, and Uranium, as well as the 
  treatment of deformed nuclei and Glauber-Gribov fluctuations of the proton in p+A 
  collisions, see https://arxiv.org/abs/1408.2549v1
 -------------------------------------------------------------------------------------
 v1.1: 
  First public release of the PHOBOS MC Glauber, see https://arxiv.org/abs/0805.4411
 -------------------------------------------------------------------------------------

 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 You should have received a copy of the GNU General Public License
 along with this program. If not, see <http://www.gnu.org/licenses/>
*/

//#define HAVE_MATHMORE

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <Riostream.h>
#include <TBits.h>
#include <TCanvas.h>
#include <TEllipse.h>
#include <TF1.h>
#include <TF2.h>
#include <TFile.h>
#include <TH2.h>
#include <TLine.h>
#include <TMath.h>
#include <TNamed.h>
#include <TNtuple.h>
#include <TObjArray.h>
#include <TRandom.h>
#include <TRotation.h>
#include <TString.h>
#include <TSystem.h>
#include <TVector3.h>
#ifdef HAVE_MATHMORE
 #include <Math/SpecFuncMathMore.h>
#endif
using namespace std;
#endif

#ifndef _runglauber_
#if !defined(__CINT__) || defined(__MAKECINT__)
#define _runglauber_ 3
#endif

#include "TGlauNucleon.hh"
#include "TGlauNucleus.hh"
#include "TGlauberMC.hh"

//---------------------------------------------------------------------------------
TF1 *getNNProf(Double_t snn=67.6, Double_t omega=0.4, Double_t G=1);

//---------------------------------------------------------------------------------
void runAndSaveNtuple(const Int_t n,
                      const char *sysA        = "Pb",
                      const char *sysB        = "Pb",
                      const Double_t signn    = 67.6,
                      const Double_t sigwidth = -1,
                      const Double_t mind     = 0.4,
		      const Double_t omega    = -1,
                      const Double_t noded    = -1,
                      const char *fname       = 0);

//---------------------------------------------------------------------------------
void runAndSaveNucleons(const Int_t n,                    
                        const char *sysA        = "Pb",           
                        const char *sysB        = "Pb",           
                        const Double_t signn    = 67.6,           
                        const Double_t sigwidth = -1,
                        const Double_t mind     = 0.4,
                        const Bool_t verbose    = 0,
			const Double_t bmin     = 0.0,
			const Double_t bmax     = 20.0,
                        const char *fname       = 0);

//---------------------------------------------------------------------------------
void runAndSmearNtuple(const Int_t n,
                       const Double_t sigs  = 0.4,
                       const char *sysA     = "p",
                       const char *sysB     = "Pb",
                       const Double_t signn = 67.6,
                       const Double_t mind  = 0.4,
		       const Double_t bmin  = 0.0,
		       const Double_t bmax  = 20.0,
                       const char *fname    = 0);


//---------------------------------------------------------------------------------
void runAndOutputLemonTree(const Int_t n,
			   const Double_t sigs  = 0.4,
			   const char *sysA     = "p",
			   const char *sysB     = "Pb",
			   const Double_t signn = 67.6,
			   const Double_t mind  = 0.4,
			   const Double_t bmin  = 0.0,
			   const Double_t bmax  = 20.0,
			   const Bool_t   ogrid = 0,
			   const char *fname    = 0);

//---------------------------------------------------------------------------------
void runAndCalcDens(const Int_t n,
		    const Double_t alpha = 0.1,
		    const char *sysA     = "Pb",
		    const char *sysB     = "Pb",
		    const Double_t signn = 67.6,
		    const Double_t mind  = 0.4,
		    const char *fname    = "glau_dens_hists.root");

//---------------------------------------------------------------------------------

//---------------------------------------------------------------------------------
TF1 *getNNProf(Double_t snn, Double_t omega, Double_t G) 
{ // NN collisoin profile from https://arxiv.org/abs/1307.0636
  if ((omega<0) || (omega>1))
    return 0;
  Double_t R2 = snn/10./TMath::Pi();
  TF1 *nnprof = new TF1("nnprofgamma","[2]*(1-TMath::Gamma([0],[1]*x^2))",0,3);
  nnprof->SetParameters(1./omega,G/omega/R2,G);
  return nnprof;
}

//---------------------------------------------------------------------------------
void runAndSaveNtuple(const Int_t n,
                      const char *sysA,
                      const char *sysB,
                      const Double_t signn,
                      const Double_t sigwidth,
                      const Double_t mind,
		      const Double_t omega,
                      const Double_t noded,
                      const char *fname)
{
  TGlauberMC *mcg=new TGlauberMC(sysA,sysB,signn,sigwidth);
  mcg->SetMinDistance(mind);
  mcg->SetNodeDistance(noded);
  mcg->SetCalcLength(0);
  mcg->SetCalcArea(0);
  mcg->SetCalcCore(0);
  mcg->SetDetail(99);
  TString om;
  if ((omega>=0) && (omega<=1)) {
    TF1 *f1 = getNNProf(signn, omega);
    mcg->SetNNProf(f1);
    om=Form("-om%.1f",omega);
  }
  TString name;
  if (fname) 
    name = fname; 
  else {
    TString nd;
    if (noded>0) 
      nd=Form("-nd%.1f",noded);
    name = Form("%s%s%s.root",mcg->Str(),om.Data(),nd.Data());
  }
  mcg->Run(n);
  TFile out(name,"recreate",name,9);
  TNtuple  *nt=mcg->GetNtuple();
  nt->Write();
  out.Close();
}

//---------------------------------------------------------------------------------
void runAndSaveNucleons(const Int_t n,                    
                        const char *sysA,           
                        const char *sysB,           
                        const Double_t signn,
                        const Double_t sigwidth,
                        const Double_t mind,
                        const Bool_t verbose,
			const Double_t bmin,
			const Double_t bmax,
			const char *fname)
{
  TGlauberMC *mcg=new TGlauberMC(sysA,sysB,signn,sigwidth);
  mcg->SetMinDistance(mind);
  mcg->SetBmin(bmin);
  mcg->SetBmax(bmax);
  TFile *out=0;
  if (fname) 
    out=new TFile(fname,"recreate",fname,9);

  for (Int_t ievent=0; ievent<n; ++ievent) {
    //get an event with at least one collision
    mcg->Run(1);
    if (ievent%100==0)
      cout << "\r" << 100.*ievent/n << "% done" << flush;

    //access, save and (if wanted) print out nucleons
    TObjArray* nucleons=mcg->GetNucleons();
    if (!nucleons) 
      continue;
    if (out)
      nucleons->Write(Form("nucleonarray%d",ievent),TObject::kSingleKey);

    if (verbose) {
      cout<<endl<<endl<<"EVENT NO: "<<ievent<<endl;
      cout<<"B = "<<mcg->GetB()<<"  Npart = "<<mcg->GetNpart()<<endl<<endl;
      printf("Nucleus\t X\t Y\t Z\tNcoll\n");
      Int_t nNucls=nucleons->GetEntries();
      for (Int_t iNucl=0; iNucl<nNucls; ++iNucl) {
        TGlauNucleon *nucl=(TGlauNucleon *)nucleons->At(iNucl);
        Char_t nucleus='A';
        if (nucl->IsInNucleusB()) 
	  nucleus='B';
        Double_t x=nucl->GetX();
        Double_t y=nucl->GetY();
        Double_t z=nucl->GetZ();
        Int_t ncoll=nucl->GetNColl();
        printf("   %c\t%2.2f\t%2.2f\t%2.2f\t%3d\n",nucleus,x,y,z,ncoll);
      }
    }
  }
  cout << endl << "Done!" << endl;
  if (out) {
    TNtuple *nt = mcg->GetNtuple();
    nt->Write();
    if (verbose)
      out->ls();
    delete out;
  }
}

//---------------------------------------------------------------------------------
void runAndSmearNtuple(const Int_t n,
                       const Double_t sigs,
                       const char *sysA,
                       const char *sysB,
                       const Double_t signn,
                       const Double_t mind,
		       const Double_t bmin,
		       const Double_t bmax,
                       const char *fname)
{
  // Run Glauber and store ntuple with smeared eccentricities in file.

  TGlauberMC *mcg = new TGlauberMC(sysA,sysB,signn);
  mcg->SetMinDistance(mind);
  mcg->SetBmin(bmin);
  mcg->SetBmax(bmax);
  
  TFile *out = TFile::Open(fname,"recreate",fname,9);
  if (!out)
    return;
  TNtuple *nt = new TNtuple("nt","nt",
      "Npart:Ncoll:B:Psi1P:Ecc1P:Psi2P:Ecc2P:Psi3P:Ecc3P:Psi4P:Ecc4P:Psi5P:Ecc5P:Psi1G:Ecc1G:Psi2G:Ecc2G:Psi3G:Ecc3G:Psi4G:Ecc4G:Psi5G:Ecc5G:Sx2P:Sy2P:Sx2G:Sy2G");
  nt->SetDirectory(out);

  const Int_t NSAMP = 100;
  TF1 *rad = new TF1("rad","x*TMath::Exp(-x*x/(2.*[0]*[0]))",0.0,3*sigs);
  rad->SetParameter(0,sigs);

  for (Int_t ievent=0; ievent<n; ++ievent) {
    while (!mcg->NextEvent()) {}

    const TGlauNucleus *nucA   = mcg->GetNucleusA();
    const TObjArray *nucleonsA = nucA->GetNucleons();
    const Int_t AN             = nucA->GetN();
    const TGlauNucleus *nucB   = mcg-> GetNucleusB();
    const TObjArray *nucleonsB = nucB->GetNucleons();
    const Int_t BN             = nucB->GetN();

    Double_t sinphi[10] = {0};
    Double_t cosphi[10] = {0};
    Double_t rn[10]     = {0};
    Double_t ecc[10]    = {0};
    Double_t psi[10]    = {0};
    Double_t sx2g       = 0;
    Double_t sy2g       = 0;

    for (Int_t s=0; s<NSAMP; ++s) {
      Int_t ni = 0;
      Double_t xvals[1000] = {0};
      Double_t yvals[1000] = {0};
      for (Int_t i = 0; i<AN; ++i) {
        TGlauNucleon *nucleonA=(TGlauNucleon*)(nucleonsA->At(i));
        if (!nucleonA->IsWounded())
          continue;
        Double_t sr = rad->GetRandom();
        Double_t sp = gRandom->Uniform(-TMath::Pi(), +TMath::Pi());
        xvals[ni]   = nucleonA->GetX() + sr*TMath::Cos(sp);
        yvals[ni]   = nucleonA->GetY() + sr*TMath::Sin(sp);
        ++ni;
      }
      for (Int_t i = 0; i<BN; ++i) {
        TGlauNucleon *nucleonB=(TGlauNucleon*)(nucleonsB->At(i));
        if (!nucleonB->IsWounded())
          continue;
        Double_t sr = rad->GetRandom();
        Double_t sp = gRandom->Uniform(-TMath::Pi(), +TMath::Pi());
        xvals[ni]   = nucleonB->GetX() + sr*TMath::Cos(sp);
        yvals[ni]   = nucleonB->GetY() + sr*TMath::Sin(sp);
        ++ni;
      }

      Double_t MeanX  = 0;
      Double_t MeanY  = 0;
      Double_t MeanX2 = 0;
      Double_t MeanY2 = 0;
      for (Int_t i = 0; i<ni; ++i) {
        MeanX  += xvals[i];
        MeanY  += yvals[i];
        MeanX2 += xvals[i]*xvals[i];
        MeanY2 += yvals[i]*yvals[i];
      }
      MeanX  /= ni;
      MeanY  /= ni;
      MeanX2 /= ni;
      MeanY2 /= ni;
      sx2g        += MeanX2-MeanX*MeanX;
      sy2g        += MeanY2-MeanY*MeanY;

      for (Int_t j = 1; j<9; ++j) {
        for (Int_t i = 0; i<ni; ++i) {
          Double_t x   = xvals[i] - MeanX;
          Double_t y   = yvals[i] - MeanY;
          Double_t r   = TMath::Sqrt(x*x+y*y);
          Double_t phi = TMath::ATan2(y,x);
          Double_t w = j;
          if (j==1)
            w = 3; // use r^3 weighting for Ecc1/Psi1
          cosphi[j] += TMath::Power(r,w)*TMath::Cos(j*phi);
          sinphi[j] += TMath::Power(r,w)*TMath::Sin(j*phi);
          rn[j]     += TMath::Power(r,w);
        }
      }
    }
    for (Int_t j = 1; j<9; ++j) {
      psi[j] = (TMath::ATan2(sinphi[j],cosphi[j]) + TMath::Pi())/j;
      ecc[j] = TMath::Sqrt(sinphi[j]*sinphi[j] + cosphi[j]*cosphi[j]) / rn[j];
    }

    Float_t v[27]; Int_t i=0;
    v[i++] = mcg->GetNpart();
    v[i++] = mcg->GetNcoll();
    v[i++] = mcg->GetB();
    v[i++] = mcg->GetPsi(1); // point-like calculation values
    v[i++] = mcg->GetEcc(1);
    v[i++] = mcg->GetPsi(2);
    v[i++] = mcg->GetEcc(2);
    v[i++] = mcg->GetPsi(3);
    v[i++] = mcg->GetEcc(3);
    v[i++] = mcg->GetPsi(4);
    v[i++] = mcg->GetEcc(4);
    v[i++] = mcg->GetPsi(5);
    v[i++] = mcg->GetEcc(5);
    v[i++] = psi[1];         // Gaussian smeared values
    v[i++] = ecc[1];
    v[i++] = psi[2];
    v[i++] = ecc[2];
    v[i++] = psi[3];
    v[i++] = ecc[3];
    v[i++] = psi[4];
    v[i++] = ecc[4];
    v[i++] = psi[5];
    v[i++] = ecc[5];
    v[i++] = mcg->GetSx2();
    v[i++] = mcg->GetSy2();
    v[i++] = sx2g/NSAMP;
    v[i++] = sy2g/NSAMP;
    nt->Fill(v);
  }

  out->Write();
  out->Close();
  delete out;
}

//---------------------------------------------------------------------------------
void runAndOutputLemonTree(const Int_t n,
                       const Double_t sigs,
                       const char *sysA,
                       const char *sysB,
                       const Double_t signn,
                       const Double_t mind,
		       const Double_t bmin,
		       const Double_t bmax,
		       const Bool_t   ogrid,
                       const char *fname)
{
  // Run Glauber and store Lemon TTree in format needed for IP-Jazma input 

  TGlauberMC *mcg = new TGlauberMC(sysA,sysB,signn);
  mcg->SetMinDistance(mind);
  mcg->SetBmin(bmin);
  mcg->SetBmax(bmax);
  
  TFile *out = TFile::Open(fname,"recreate",fname,9);
  if (!out) return;

  // create new TTree with MC Glauber information for input to IP-Jazma
  const Int_t lemonmaxNucleons = 500;
  Int_t       lemonnpart;
  Int_t       lemonncoll;
  Int_t       lemonnparta;
  Int_t       lemonnpartb;  
  Float_t     lemonb;                        // collision impact parameter
  Float_t     lemoneccgaus[10];
  Float_t     lemoneccpoint[10];  
  Int_t       lemonnproj;
  Int_t       lemonntarg;
  Float_t     lemonxproj[lemonmaxNucleons];  // x,y,z coordinates for all nucleons
  Float_t     lemonyproj[lemonmaxNucleons];  // note these must be in the global coordinate frame
  Float_t     lemonzproj[lemonmaxNucleons];
  Float_t     lemonxtarg[lemonmaxNucleons];  // x,y,z coordinates for all nucleons
  Float_t     lemonytarg[lemonmaxNucleons];  // note these must be in the global coordinate frame
  Float_t     lemonztarg[lemonmaxNucleons];
  
  TTree *lemon = new TTree("lemon","lemon");
  lemon->Branch("npart",&lemonnpart,"npart/I");
  lemon->Branch("nparta",&lemonnparta,"nparta/I");
  lemon->Branch("npartb",&lemonnpartb,"npartb/I");  
  lemon->Branch("ncoll",&lemonncoll,"ncoll/I");
  lemon->Branch("b",&lemonb,"b/F");
  lemon->Branch("eccgaus",lemoneccgaus,"eccgaus[10]/F");
  lemon->Branch("eccpoint",lemoneccpoint,"eccpoint[10]/F");  
  lemon->Branch("nproj",&lemonnproj,"nproj/I");
  lemon->Branch("ntarg",&lemonntarg,"ntarg/I");  
  lemon->Branch("xproj",lemonxproj,"xproj[500]/F");
  lemon->Branch("yproj",lemonyproj,"yproj[500]/F");
  lemon->Branch("zproj",lemonyproj,"zproj[500]/F");  
  lemon->Branch("xtarg",lemonxtarg,"xtarg[500]/F");
  lemon->Branch("ytarg",lemonytarg,"ytarg[500]/F");
  lemon->Branch("ztarg",lemonytarg,"ztarg[500]/F");    
  lemon->SetDirectory(out);

  const Int_t NSAMP = 100;
  TF1 *rad = new TF1("rad","x*TMath::Exp(-x*x/(2.*[0]*[0]))",0.0,3*sigs);
  rad->SetParameter(0,sigs);
  TF2* smearing_function = new TF2("smear_tf2", "TMath::Exp(-(x*x+y*y)/(2.*[0]*[0]))/(2*TMath::Pi()*[0]*[0])", 0, 10*sigs, 0, 10*sigs);
  smearing_function->SetParameter(0,sigs);
  
  for (Int_t ievent=0; ievent<n; ++ievent) {
    while (!mcg->NextEvent()) {}

    const TGlauNucleus *nucA   = mcg->GetNucleusA();
    const TObjArray *nucleonsA = nucA->GetNucleons();
    const Int_t AN             = nucA->GetN();
    const TGlauNucleus *nucB   = mcg-> GetNucleusB();
    const TObjArray *nucleonsB = nucB->GetNucleons();
    const Int_t BN             = nucB->GetN();

    Double_t sinphi[10] = {0};
    Double_t cosphi[10] = {0};
    Double_t rn[10]     = {0};
    Double_t ecc[10]    = {0};
    Double_t psi[10]    = {0};
    Double_t sx2g       = 0;
    Double_t sy2g       = 0;

    for (Int_t s=0; s<NSAMP; ++s) {
      Int_t ni = 0;
      Double_t xvals[1000] = {0};
      Double_t yvals[1000] = {0};
      for (Int_t i = 0; i<AN; ++i) {
        TGlauNucleon *nucleonA=(TGlauNucleon*)(nucleonsA->At(i));
        if (!nucleonA->IsWounded())
          continue;
        Double_t sr = rad->GetRandom();
        Double_t sp = gRandom->Uniform(-TMath::Pi(), +TMath::Pi());
        xvals[ni]   = nucleonA->GetX() + sr*TMath::Cos(sp);
        yvals[ni]   = nucleonA->GetY() + sr*TMath::Sin(sp);
        ++ni;
      }
      for (Int_t i = 0; i<BN; ++i) {
        TGlauNucleon *nucleonB=(TGlauNucleon*)(nucleonsB->At(i));
        if (!nucleonB->IsWounded())
          continue;
        Double_t sr = rad->GetRandom();
        Double_t sp = gRandom->Uniform(-TMath::Pi(), +TMath::Pi());
        xvals[ni]   = nucleonB->GetX() + sr*TMath::Cos(sp);
        yvals[ni]   = nucleonB->GetY() + sr*TMath::Sin(sp);
        ++ni;
      }

      Double_t MeanX  = 0;
      Double_t MeanY  = 0;
      Double_t MeanX2 = 0;
      Double_t MeanY2 = 0;
      for (Int_t i = 0; i<ni; ++i) {
        MeanX  += xvals[i];
        MeanY  += yvals[i];
        MeanX2 += xvals[i]*xvals[i];
        MeanY2 += yvals[i]*yvals[i];
      }
      MeanX  /= ni;
      MeanY  /= ni;
      MeanX2 /= ni;
      MeanY2 /= ni;
      sx2g        += MeanX2-MeanX*MeanX;
      sy2g        += MeanY2-MeanY*MeanY;

      for (Int_t j = 1; j<9; ++j) {
        for (Int_t i = 0; i<ni; ++i) {
          Double_t x   = xvals[i] - MeanX;
          Double_t y   = yvals[i] - MeanY;
          Double_t r   = TMath::Sqrt(x*x+y*y);
          Double_t phi = TMath::ATan2(y,x);
          Double_t w = j;
          if (j==1)
            w = 3; // use r^3 weighting for Ecc1/Psi1
          cosphi[j] += TMath::Power(r,w)*TMath::Cos(j*phi);
          sinphi[j] += TMath::Power(r,w)*TMath::Sin(j*phi);
          rn[j]     += TMath::Power(r,w);
        }
      }
    }
    for (Int_t j = 1; j<9; ++j) {
      psi[j] = (TMath::ATan2(sinphi[j],cosphi[j]) + TMath::Pi())/j;
      ecc[j] = TMath::Sqrt(sinphi[j]*sinphi[j] + cosphi[j]*cosphi[j]) / rn[j];
    }

    // fill lemon TTree variables for this event

    lemonnpart   = mcg->GetNpart();
    lemonnparta  = mcg->GetNpartA();
    lemonnpartb  = mcg->GetNpartB();    
    lemonncoll   = mcg->GetNcoll();
    lemonb       = mcg->GetB(); 

    lemoneccpoint[0] = 0.0;
    lemoneccpoint[1] = mcg->GetEcc(1);
    lemoneccpoint[2] = mcg->GetEcc(2);
    lemoneccpoint[3] = mcg->GetEcc(3);
    lemoneccpoint[4] = mcg->GetEcc(4);
    lemoneccpoint[5] = mcg->GetEcc(5);
    lemoneccpoint[6] = mcg->GetEcc(6);    
    lemoneccpoint[7] = mcg->GetEcc(7);    
    lemoneccpoint[8] = mcg->GetEcc(8);    
    lemoneccpoint[9] = mcg->GetEcc(9);    

    lemoneccgaus[0] = 0.0;
    lemoneccgaus[1] = ecc[1];    
    lemoneccgaus[2] = ecc[2];
    lemoneccgaus[3] = ecc[3];
    lemoneccgaus[4] = ecc[4];
    lemoneccgaus[5] = ecc[5];
    lemoneccgaus[6] = ecc[6];
    lemoneccgaus[7] = ecc[7];
    lemoneccgaus[8] = ecc[8];
    lemoneccgaus[9] = ecc[9];

    for (Int_t i = 0; i<AN; ++i) {
      TGlauNucleon *nucleonA=(TGlauNucleon*)(nucleonsA->At(i));
      lemonxproj[i] = nucleonA->GetX();
      lemonyproj[i] = nucleonA->GetY();
      lemonzproj[i] = nucleonA->GetZ();      
    }
    for (Int_t i = 0; i<BN; ++i) {
      TGlauNucleon *nucleonB=(TGlauNucleon*)(nucleonsB->At(i));
      lemonxtarg[i] = nucleonB->GetX();
      lemonytarg[i] = nucleonB->GetY();
      lemonztarg[i] = nucleonB->GetZ();      
    }
    
    lemon->Fill();
    
    //======================================================================
    // also include option to write out nucleon smeared energy density map
    //======================================================================

    if (ogrid) {

      const Int_t nbins = 1000; 
      const Int_t nbinsx = nbins;
      const Int_t nbinsy = nbins;
      const Double_t max_x = 7.5;
      
      // now create an energy density distribution (a.u.)
      TH2D* inited_hist = new TH2D(Form("inited_event%i",ievent), ";x;y;E [a.u.]", nbinsx, -max_x, max_x, nbinsy, -max_x, max_x);
      for (Int_t ybin=1; ybin<=nbinsy; ybin++) {
	for (Int_t xbin=1; xbin<=nbinsx; xbin++) {
	  
	  const Double_t xval = inited_hist->GetXaxis()->GetBinCenter(xbin);
	  const Double_t yval = inited_hist->GetYaxis()->GetBinCenter(ybin);
	  long double content = 0.;  // sum the contributions from all wounded nucleons
	  
	  for (Int_t i = 0; i<AN; ++i) {
	    TGlauNucleon *nucleonA=(TGlauNucleon*)(nucleonsA->At(i));
	    if (!nucleonA->IsWounded()) continue;   // skip non-wounded nucleons
	    content += smearing_function->Eval(nucleonA->GetX() - xval, nucleonA->GetY() - yval);
	  }
	  for (Int_t i = 0; i<BN; ++i) {
	    TGlauNucleon *nucleonB=(TGlauNucleon*)(nucleonsB->At(i));
	    if (!nucleonB->IsWounded()) continue;   // skip non-wounded nucleons
	    content += smearing_function->Eval(nucleonB->GetX() - xval, nucleonB->GetY() - yval);
	  }
	  inited_hist->SetBinContent(xbin, ybin, content);
	  
	}
      }
      inited_hist->Write();
      if (inited_hist) delete inited_hist;
    }
  } // end loop over events

  out->Write();
  out->Close();
  delete out;
}

//---------------------------------------------------------------------------------
void runAndCalcDens(const Int_t n,
		    const Double_t alpha,
		    const char *sysA,
		    const char *sysB,
		    const Double_t signn,
		    const Double_t mind,
		    const char *fname)
{
  // Run Glauber and store per event a density profile in x and y, calculated from participant and binary positions
  // with relative weight given by alpha.

  TGlauberMC *mcg = new TGlauberMC(sysA,sysB,signn);
  mcg->SetMinDistance(mind);

  TFile *out = 0;
  TCanvas *c1 = 0;
  if (fname) {
    out = TFile::Open(fname,"recreate",fname,9);
    if (!out)
      return;
  } else {
    c1 = new TCanvas;
  }

  const Int_t NSAMP = 100;
  const Double_t wp = (1-alpha)/2;
  const Double_t wb = alpha;
  const Double_t sigs = TMath::Sqrt(signn/20/TMath::Pi()); //from arXiv::0711.3724


  TF1 *rad = new TF1("rad","x*TMath::Exp(-x*x/(2.*[0]*[0]))",0.0,3*sigs);
  rad->SetParameter(0,sigs);

  TH2F *h2f = new TH2F("h2f",";x (fm);y (fm)",121,-15.5,15.5,121,-15.5,15.5);
  h2f->SetStats(0);
  
  for (Int_t ievent=0; ievent<n; ++ievent) {
    while (!mcg->NextEvent()) {}
    h2f->Reset();
    h2f->SetName(Form("Event_%d",ievent));
    h2f->SetTitle(Form("Npart=%d, Ncoll=%d",mcg->GetNpart(), mcg->GetNcoll()));

    const TGlauNucleus *nucA   = mcg->GetNucleusA();
    const TObjArray *nucleonsA = nucA->GetNucleons();
    const Int_t AN             = nucA->GetN();
    const TGlauNucleus *nucB   = mcg-> GetNucleusB();
    const TObjArray *nucleonsB = nucB->GetNucleons();
    const Int_t BN             = nucB->GetN();

    for (Int_t i = 0; i<AN; ++i) {
      TGlauNucleon *nucleonA=(TGlauNucleon*)(nucleonsA->At(i));
      if (!nucleonA->IsWounded())
	continue;
      Double_t xA=nucleonA->GetX();
      Double_t yA=nucleonA->GetY();
      for (Int_t s=0; s<NSAMP; ++s) {
	Double_t sr = rad->GetRandom();
	Double_t sp = gRandom->Uniform(-TMath::Pi(), +TMath::Pi());
	h2f->Fill(xA+sr*TMath::Cos(sp),yA+sr*TMath::Sin(sp),wp);
      }
    }

    for (Int_t j = 0; j<BN; ++j) {
      TGlauNucleon *nucleonB=(TGlauNucleon*)(nucleonsB->At(j));
      if (!nucleonB->IsWounded())
	continue;
      Double_t xB=nucleonB->GetX();
      Double_t yB=nucleonB->GetY();
      for (Int_t s=0; s<NSAMP; ++s) {
	Double_t sr = rad->GetRandom();
	Double_t sp = gRandom->Uniform(-TMath::Pi(), +TMath::Pi());
	h2f->Fill(xB+sr*TMath::Cos(sp),yB+sr*TMath::Sin(sp),wp);
      }
    }

    if (alpha>0) {
      for (Int_t i = 0; i<AN; ++i) {
	TGlauNucleon *nucleonA=(TGlauNucleon*)(nucleonsA->At(i));
	if (!nucleonA->IsWounded())
	  continue;
	Double_t xA=nucleonA->GetX();
	Double_t yA=nucleonA->GetY();
	for (Int_t j = 0; j<BN; ++j) {
	  TGlauNucleon *nucleonB=(TGlauNucleon*)(nucleonsB->At(j));
	  if (!mcg->IsBC(i,j))
	    continue;
	  Double_t xB=nucleonB->GetX();
	  Double_t yB=nucleonB->GetY();
	  Double_t dX=(xA+xB)/2;
	  Double_t dY=(yA+yB)/2;
	  for (Int_t s=0; s<NSAMP; ++s) {
	    Double_t sr = rad->GetRandom();
	    Double_t sp = gRandom->Uniform(-TMath::Pi(), +TMath::Pi());
	    h2f->Fill(dX+sr*TMath::Cos(sp),dY+sr*TMath::Sin(sp),wb);
	  }
	}
      }
    }
    h2f->Scale(1./h2f->Integral());
    if (out) {
      h2f->Write();
    } else {
      h2f->Draw("colz");
      c1->Update();
      gSystem->Sleep(1000);
    }
  }
  if (out) {
    out->Write();
    out->Close();
    delete out;
  }
}

#endif  //end of ifndef runglauber