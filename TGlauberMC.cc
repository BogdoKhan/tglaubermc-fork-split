#include "TGlauberMC.hh"

//---------------------------------------------------------------------------------
  ClassImp(TGlauberMC)
  ClassImp(TGlauberMC::Event)
//---------------------------------------------------------------------------------

TGlauberMC::TGlauberMC(const char* NA, const char* NB, Double_t xsect, Double_t xsectsigma) :
  fANucleus(NA),fBNucleus(NB),
  fXSect(xsect),fXSectOmega(0),fXSectLambda(0),fXSectEvent(0),
  fNucleonsA(0),fNucleonsB(0),fNucleons(0),
  fAN(0),fBN(0),fNt(0),
  fEvents(0),fTotalEvents(0),fBmin(0),fBmax(20),fHardFrac(0.65),
  fDetail(99),fCalcArea(0),fCalcLength(0), fDoCore(0), fDoAAGG(1),
  fMaxNpartFound(0),f2Cx(0),fPTot(0),fNNProf(0),
  fEv()
{
  if (xsectsigma>0) {
    fXSectOmega = xsectsigma;
    fXSectLambda = 1;
    fPTot = new TF1("fPTot","((x/[2])/(x/[2]+[0]))*exp(-(((x/[2])/[0]-1 )**2)/([1]*[1]))/[2]",0,300);
    fPTot->SetParameters(fXSect,fXSectOmega,fXSectLambda);
    fPTot->SetNpx(1000);
    fXSectLambda = fXSect/fPTot->GetHistogram()->GetMean();
    cout << "final lambda=" << fXSectLambda << endl;
    fPTot->SetParameters(fXSect,fXSectOmega,fXSectLambda);
    cout << "final <sigma>=" << fPTot->GetHistogram()->GetMean() << endl;
  }

  TString name(Form("Glauber_%s_%s",fANucleus.GetName(),fBNucleus.GetName()));
  TString title(Form("Glauber %s+%s Version",fANucleus.GetName(),fBNucleus.GetName()));
  SetName(name);
  SetTitle(title);
}

Bool_t TGlauberMC::CalcEvent(Double_t bgen)
{
  // calc next event
  if (!fNucleonsA) {
    fNucleonsA = fANucleus.GetNucleons();
    fAN = fANucleus.GetN();
    for (Int_t i = 0; i<fAN; ++i) {
      TGlauNucleon *nucleonA=(TGlauNucleon*)(fNucleonsA->At(i));
      nucleonA->SetInNucleusA();
    }
  }

  if (!fNucleonsB) {
    fNucleonsB = fBNucleus.GetNucleons();
    fBN = fBNucleus.GetN();
    for (Int_t i = 0; i<fBN; ++i) {
      TGlauNucleon *nucleonB=(TGlauNucleon*)(fNucleonsB->At(i));
      nucleonB->SetInNucleusB();
    }
  }

  Double_t xsecA[999] = {0};
  Double_t xsecB[999] = {0};
  if (fPTot) {
    fXSectEvent = fPTot->GetRandom();
    if (fDoAAGG) {
      for (Int_t i = 0; i<fAN; ++i)
	xsecA[i] = fPTot->GetRandom();
      for (Int_t i = 0; i<fBN; ++i)
	xsecB[i] = fPTot->GetRandom();
    }
  } else 
    fXSectEvent = fXSect;

  // "ball" diameter = distance at which two balls interact
  Double_t d2 = (Double_t)fXSectEvent/(TMath::Pi()*10); // in fm^2
  Double_t bh = TMath::Sqrt(d2*fHardFrac);
  if (fNNProf) {
    Double_t xmin=0,xmax=0;
    fNNProf->GetRange(xmin,xmax);
    d2 = xmax*xmax;
  }

  fEv.Reset();
  memset(fBC,0,sizeof(Bool_t)*999*999);
  Int_t nc=0,nh=0;
  for (Int_t i = 0; i<fBN; ++i) {
    TGlauNucleon *nucleonB=(TGlauNucleon*)(fNucleonsB->At(i));
    Bool_t tB=nucleonB->GetType();
    for (Int_t j = 0; j<fAN; ++j) {
      TGlauNucleon *nucleonA=(TGlauNucleon*)(fNucleonsA->At(j));
      Double_t dx = nucleonB->GetX()-nucleonA->GetX();
      Double_t dy = nucleonB->GetY()-nucleonA->GetY();
      Double_t dij = dx*dx+dy*dy;
      if (fDoAAGG && fPTot)
	d2 = 0.5*(xsecA[j]+xsecB[i])/(TMath::Pi()*10);
      if (dij>d2) 
        continue;
      Double_t bij = TMath::Sqrt(dij);
      if (fNNProf) {
        Double_t val = fNNProf->Eval(bij);
        Double_t ran = gRandom->Uniform();
        if (ran>val)
          continue;
      }
      nucleonB->Collide();
      nucleonA->Collide();
      fBC[i][j] = 1;
      fEv.BNN  += bij;
      ++nc;
      if (bij<bh)
        ++nh;
      Bool_t tA=nucleonA->GetType();
      if (tA!=tB)
        ++fEv.Ncollpn;
      else if (tA==1)
        ++fEv.Ncollpp;
      else
        ++fEv.Ncollnn;
      if (nc==1) {
        fEv.X0 = (nucleonA->GetX()+nucleonB->GetX())/2;
        fEv.Y0 = (nucleonA->GetY()+nucleonB->GetY())/2;
      }
    }
  }
  fEv.B = bgen;
  ++fTotalEvents;
  if (nc>0) {
    ++fEvents;
    fEv.Ncoll     = nc;
    fEv.Nhard     = nh;
    fEv.BNN      /= nc;
    return CalcResults(bgen);
  }
  return kFALSE;
}

Bool_t TGlauberMC::CalcResults(Double_t bgen)
{
  // calc results for the given event
  Double_t sumW=0;
  Double_t sumWA=0;
  Double_t sumWB=0;

  Double_t sinphi[10] = {0};
  Double_t cosphi[10] = {0};
  Double_t rn[10]     = {0};

  const Int_t kNc = fDoCore; // used later for core/corona

  for (Int_t i = 0; i<fAN; ++i) {
    TGlauNucleon *nucleonA=(TGlauNucleon*)(fNucleonsA->At(i));
    Double_t xA=nucleonA->GetX();
    Double_t yA=nucleonA->GetY();
    fEv.MeanXSystem  += xA;
    fEv.MeanYSystem  += yA;
    fEv.MeanXA  += xA;
    fEv.MeanYA  += yA;
    if (nucleonA->IsWounded()) {
      Double_t w = nucleonA->Get2CWeight(f2Cx);
      ++fEv.Npart;
      if (nucleonA->GetNColl()==1)
	++fEv.Npart0;
      ++fEv.NpartA;
      sumW   += w;
      sumWA  += w;
      fEv.MeanX  += xA * w;
      fEv.MeanY  += yA * w;
      fEv.MeanX2 += xA * xA * w;
      fEv.MeanY2 += yA * yA * w;
      fEv.MeanXY += xA * yA * w;
    }
  }

  for (Int_t i = 0; i<fBN; ++i) {
    TGlauNucleon *nucleonB=(TGlauNucleon*)(fNucleonsB->At(i));
    Double_t xB=nucleonB->GetX();
    Double_t yB=nucleonB->GetY();
    fEv.MeanXSystem  += xB;
    fEv.MeanYSystem  += yB;
    fEv.MeanXB  += xB;
    fEv.MeanYB  += yB;
    if (nucleonB->IsWounded()) {
      Double_t w = nucleonB->Get2CWeight(f2Cx);
      ++fEv.Npart;
      if (nucleonB->GetNColl()==1)
	++fEv.Npart0;
      ++fEv.NpartB;
      sumW   += w;
      sumWB  += w;
      fEv.MeanX  += xB * w;
      fEv.MeanY  += yB * w;
      fEv.MeanX2 += xB * xB * w;
      fEv.MeanY2 += yB * yB * w;
      fEv.MeanXY += xB * yB * w;
    }
  }
  if (fEv.Npart>0) {
    fEv.MeanX  /= sumW;
    fEv.MeanY  /= sumW;
    fEv.MeanX2 /= sumW;
    fEv.MeanY2 /= sumW;
    fEv.MeanXY /= sumW;
  } else {
    fEv.MeanX = 0;
    fEv.MeanY  = 0;
    fEv.MeanX2 = 0;
    fEv.MeanY2 = 0;
    fEv.MeanXY = 0;
  }

  if (fAN+fBN>0) {
    fEv.MeanXSystem /= (fAN + fBN);
    fEv.MeanYSystem /= (fAN + fBN);
  } else {
    fEv.MeanXSystem = 0;
    fEv.MeanYSystem = 0;
  }
  if (fAN>0) {
    fEv.MeanXA /= fAN;
    fEv.MeanYA /= fAN;
  } else {
    fEv.MeanXA = 0;
    fEv.MeanYA = 0;
  }
  if (fBN>0) {
    fEv.MeanXB /= fBN;
    fEv.MeanYB /= fBN;
  } else {
    fEv.MeanXB = 0;
    fEv.MeanYB = 0;
  }

  fEv.VarX  = fEv.MeanX2-(fEv.MeanX*fEv.MeanX);
  fEv.VarY  = fEv.MeanY2-(fEv.MeanY*fEv.MeanY);
  fEv.VarXY = fEv.MeanXY-fEv.MeanX*fEv.MeanY;
  Double_t tmpa = fEv.VarX*fEv.VarY-fEv.VarXY*fEv.VarXY;
  if (tmpa<0) 
    fEv.AreaW = -1;
  else 
    fEv.AreaW = TMath::Sqrt(tmpa);

  if (fEv.Npart>0) {
    // do full moments relative to meanX and meanY
    for (Int_t n = 1; n<10; ++n) {
      for (Int_t ia = 0; ia<fAN; ++ia) {
        TGlauNucleon *nucleonA=(TGlauNucleon*)(fNucleonsA->At(ia));
	if (nucleonA->GetNColl()<=kNc) 
	  continue;
	Double_t xA=nucleonA->GetX() - fEv.MeanX;
	Double_t yA=nucleonA->GetY() - fEv.MeanY;
	Double_t r = TMath::Sqrt(xA*xA+yA*yA);
	Double_t phi = TMath::ATan2(yA,xA);
	Double_t w = n;
	if (n==1) 
	  w = 3; // use r^3 weighting for Ecc1/Psi1
	Double_t rw = TMath::Power(r,w);
	cosphi[n] += rw*TMath::Cos(n*phi);
	sinphi[n] += rw*TMath::Sin(n*phi);
	rn[n] += rw;
      }
      for (Int_t ib = 0; ib<fBN; ++ib) {
        TGlauNucleon *nucleonB=(TGlauNucleon*)(fNucleonsB->At(ib));
	if (nucleonB->GetNColl()<=kNc) 
	  continue;
	Double_t xB=nucleonB->GetX() - fEv.MeanX;
	Double_t yB=nucleonB->GetY() - fEv.MeanY;
	Double_t r = TMath::Sqrt(xB*xB+yB*yB);
	Double_t phi = TMath::ATan2(yB,xB);
	Double_t w = n;
	if (n==1)
	  w = 3; // use r^3 weighting for Ecc1/Psi1
	Double_t rw = TMath::Power(r,w);
	cosphi[n] += rw*TMath::Cos(n*phi);
	sinphi[n] += rw*TMath::Sin(n*phi);
	rn[n] += rw;
      }
      cosphi[n] /= fEv.Npart;
      sinphi[n] /= fEv.Npart;
      rn[n] /= fEv.Npart;
      if (rn[n]>0) {
	fPsiN[n] = (TMath::ATan2(sinphi[n],cosphi[n]) + TMath::Pi())/n;
	fEccN[n] = TMath::Sqrt(sinphi[n]*sinphi[n]+cosphi[n]*cosphi[n])/rn[n];
      } else {
	fPsiN[n] = -1;
	fEccN[n] = -1;
      }
    }
    if (!kNc) { //silly test but useful to catch errors 
      Double_t t=TMath::Sqrt(TMath::Power(fEv.VarY-fEv.VarX,2)+4.*fEv.VarXY*fEv.VarXY)/(fEv.VarY+fEv.VarX)/fEccN[2];
      if (t<0.99||t>1.01)
        cout << "Error: Expected t=1 but found t=" << t << endl;
    }
  }

  fEv.B      = bgen;
  fEv.PhiA   = fANucleus.GetPhiRot();
  fEv.ThetaA = fANucleus.GetThetaRot();
  fEv.PhiB   = fBNucleus.GetPhiRot();
  fEv.ThetaB = fBNucleus.GetThetaRot();
  fEv.Psi1   = fPsiN[1];
  fEv.Ecc1   = fEccN[1];
  fEv.Psi2   = fPsiN[2];
  fEv.Ecc2   = fEccN[2];
  fEv.Psi3   = fPsiN[3];
  fEv.Ecc3   = fEccN[3];
  fEv.Psi4   = fPsiN[4];
  fEv.Ecc4   = fEccN[4];
  fEv.Psi5   = fPsiN[5];
  fEv.Ecc5   = fEccN[5];

  if (fCalcArea) {
    const Int_t nbins=200;
    const Double_t ell=10;
    const Double_t da=2*ell*2*ell/nbins/nbins;
    const Double_t d2 = (Double_t)fXSectEvent/(TMath::Pi()*10); // in fm^2
    const Double_t r2 = d2/4.;
    const Double_t mx = fEv.MeanX;
    const Double_t my = fEv.MeanY;
    TH2D areaA("hAreaA",";x (fm);y (fm)",nbins,-ell,ell,nbins,-ell,ell);
    TH2D areaB("hAreaB",";x (fm);y (fm)",nbins,-ell,ell,nbins,-ell,ell);
    for (Int_t i = 0; i<fAN; ++i) {
      TGlauNucleon *nucleonA=(TGlauNucleon*)(fNucleonsA->At(i));
      if (!nucleonA->IsWounded())
        continue;
      if (nucleonA->GetNColl()==kNc)
        continue;
      Double_t x = nucleonA->GetX()-mx;
      Double_t y = nucleonA->GetY()-my;
      for (Int_t xi=1; xi<=nbins; ++xi) {
        for (Int_t yi=1; yi<=nbins; ++yi) {
          Int_t bin = areaA.GetBin(xi,yi);
          Double_t val=areaA.GetBinContent(bin);
          if (val>0)
            continue;
          Double_t dx=x-areaA.GetXaxis()->GetBinCenter(xi);
          Double_t dy=y-areaA.GetYaxis()->GetBinCenter(yi);
          if (dx*dx+dy*dy<r2)
            areaA.SetBinContent(bin,1);
        }
      }
    }
    for (Int_t i = 0; i<fBN; ++i) {
      TGlauNucleon *nucleonB=(TGlauNucleon*)(fNucleonsB->At(i));
      if (!nucleonB->IsWounded())
        continue;
      if (nucleonB->GetNColl()==kNc)
        continue;
      Double_t x = nucleonB->GetX()-mx;
      Double_t y = nucleonB->GetY()-my;
      for (Int_t xi=1; xi<=nbins; ++xi) {
        for (Int_t yi=1; yi<=nbins; ++yi) {
          Int_t bin = areaB.GetBin(xi,yi);
          Double_t val=areaB.GetBinContent(bin);
          if (val>0)
            continue;
          Double_t dx=x-areaB.GetXaxis()->GetBinCenter(xi);
          Double_t dy=y-areaB.GetYaxis()->GetBinCenter(yi);
          if (dx*dx+dy*dy<r2)
            areaB.SetBinContent(bin,1);
        }
      }
    }
    Double_t overlap1=0;
    Double_t overlap2=0;
    for (Int_t xi=1; xi<=nbins; ++xi) {
      for (Int_t yi=1; yi<=nbins; ++yi) {
        Int_t bin = areaA.GetBin(xi,yi);
        Double_t vA=areaA.GetBinContent(bin);
        Double_t vB=areaB.GetBinContent(bin);
        if (vA>0&&vB>0)
          ++overlap1;
        if (vA>0||vB>0)
          ++overlap2;
      }
    }
    fEv.AreaO = overlap1*da;
    fEv.AreaA = overlap2*da;
  }

  if (fCalcLength) {
    const Double_t krhs = TMath::Sqrt(fXSectEvent/40./TMath::Pi());
    const Double_t ksg  = krhs/TMath::Sqrt(5);
    const Double_t kDL  = 0.1;
    TF1 rad("rad","2*pi/[0]/[0]*TMath::Exp(-x*x/(2.*[0]*[0]))",0.0,5*ksg); 
    rad.SetParameter(0,ksg);
    const Double_t minval = rad.Eval(5*ksg);
    fEv.Phi0         = gRandom->Uniform(0,TMath::TwoPi());
    Double_t kcphi0  = TMath::Cos(fEv.Phi0);
    Double_t ksphi0  = TMath::Sin(fEv.Phi0);
    Double_t x       = fEv.X0;
    Double_t y       = fEv.Y0;
    Double_t i0a     = 0;
    Double_t i1a     = 0;
    Double_t l       = 0;
    Double_t val     = CalcDens(rad,x,y);
    while (val>minval) {
      x     += kDL * kcphi0;
      y     += kDL * ksphi0;
      i0a   += val;
      i1a   += l*val;
      l+=kDL;
      val    = CalcDens(rad,x,y);
    }
    fEv.Length = 2*i1a/i0a;
  }

  if (fEv.Npart > fMaxNpartFound) 
    fMaxNpartFound = fEv.Npart;

  return kTRUE;
}

Double_t TGlauberMC::CalcDens(TF1 &prof, Double_t xval, Double_t yval) const
{
  Double_t rmin=0,rmax=0;
  prof.GetRange(rmin,rmax);
  Double_t r2max = rmax*rmax;
  Double_t ret = 0;
  for (Int_t i = 0; i<fAN; ++i) {
    TGlauNucleon *nucleonA=(TGlauNucleon*)(fNucleonsA->At(i));
    if (!nucleonA->IsWounded())
      continue;
    Double_t x = nucleonA->GetX();
    Double_t y = nucleonA->GetY();
    Double_t r2=(xval-x)*(xval-x)+(yval-y)*(yval-y);
    if (r2>r2max)
      continue;
    ret += prof.Eval(TMath::Sqrt(r2));
  }
  for (Int_t i = 0; i<fBN; ++i) {
    TGlauNucleon *nucleonB=(TGlauNucleon*)(fNucleonsB->At(i));
    if (!nucleonB->IsWounded())
      continue;
    Double_t x = nucleonB->GetX();
    Double_t y = nucleonB->GetY();
    Double_t r2=(xval-x)*(xval-x)+(yval-y)*(yval-y);
    if (r2>r2max)
      continue;
    ret += prof.Eval(TMath::Sqrt(r2));
  }
  return ret;
}

void TGlauberMC::Draw(Option_t* option)
{
  static TH2F *h2f = new TH2F("hGlauberMC",";x (fm);y(fm)",1,-18,18,1,-12,12);

  h2f->Reset();
  h2f->SetStats(0);
  h2f->Draw();

  TEllipse e;
  e.SetFillColor(0);
  e.SetFillStyle(0);
  e.SetLineColor(1);
  e.SetLineStyle(2);
  e.SetLineWidth(1);
  e.DrawEllipse(GetB()/2,0,fBNucleus.GetR(),fBNucleus.GetR(),0,360,0);
  e.DrawEllipse(-GetB()/2,0,fANucleus.GetR(),fANucleus.GetR(),0,360,0);
  fANucleus.Draw(fXSect, kMagenta, kYellow);
  fBNucleus.Draw(fXSect, kMagenta, kOrange);

  TString opt(option);
  if (opt.IsNull())
    return;

  Double_t sy2 = GetSy2();
  Double_t sx2 = GetSx2();
  Double_t phase = 0;
  if (sy2<sx2) {
    Double_t d = sx2;
    sx2 = sy2;
    sy2 = d;
    phase = TMath::Pi()/2.;
  }
  Double_t x1 = (0.5*(sy2-sx2)+TMath::Sqrt(TMath::Power(sy2-sx2,2.)-4*TMath::Power(GetSxy(),2)));
  Double_t ang = TMath::ATan2(-GetSxy(),x1)+phase;
  TLine l;
  l.SetLineWidth(3);
  l.DrawLine(-10*TMath::Cos(ang),-10*TMath::Sin(ang),10*TMath::Cos(ang),10*TMath::Sin(ang));
}

Double_t TGlauberMC::GetTotXSect() const
{
  return (1.*fEvents/fTotalEvents)*TMath::Pi()*fBmax*fBmax/100;
}

Double_t TGlauberMC::GetTotXSectErr() const
{
  return GetTotXSect()/TMath::Sqrt((Double_t)fEvents) * TMath::Sqrt(Double_t(1.-fEvents/fTotalEvents));
}

TObjArray *TGlauberMC::GetNucleons() 
{
  if (!fNucleonsA || !fNucleonsB) return 0;
  if (fNucleons) return fNucleons;

  fNucleonsA->SetOwner(0);
  fNucleonsB->SetOwner(0);
  TObjArray *allnucleons=new TObjArray(fAN+fBN);
  allnucleons->SetOwner();
  for (Int_t i = 0; i<fAN; ++i) {
    allnucleons->Add(fNucleonsA->At(i));
  }
  for (Int_t i = 0; i<fBN; ++i) {
    allnucleons->Add(fNucleonsB->At(i));
  }
  fNucleons = allnucleons;
  return allnucleons;
}

Bool_t TGlauberMC::NextEvent(Double_t bgen)
{
  if (bgen<0) 
    bgen = TMath::Sqrt((fBmax*fBmax-fBmin*fBmin)*gRandom->Rndm()+fBmin*fBmin);

  fANucleus.ThrowNucleons(-bgen/2.);
  fBNucleus.ThrowNucleons(bgen/2.);
  return CalcEvent(bgen);
}

Bool_t TGlauberMC::ReadNextEvent(Bool_t calc, const char *fname)
{
  static TFile *inf = 0;
  static Int_t iev  = 0;
  if (fname) {
    cout << "ReadNextEvent: Setting up file " << fname << endl;
    delete inf;
    inf = TFile::Open(fname);
    if (!inf) 
      return 0;
    if (!fNucleonsA) {
      fANucleus.ThrowNucleons();
      fNucleonsA = fANucleus.GetNucleons();
      fAN = fANucleus.GetN();
      for (Int_t i = 0; i<fAN; ++i) {
	TGlauNucleon *nucleonA=(TGlauNucleon*)(fNucleonsA->At(i));
	nucleonA->SetInNucleusA();
      }
    }
    if (!fNucleonsB) {
      fBNucleus.ThrowNucleons();
      fNucleonsB = fBNucleus.GetNucleons();
      fBN = fBNucleus.GetN();
      for (Int_t i = 0; i<fBN; ++i) {
	TGlauNucleon *nucleonB=(TGlauNucleon*)(fNucleonsB->At(i));
	nucleonB->SetInNucleusB();
      }
    }
    if (calc)
      return 1;
    fNt = dynamic_cast<TNtuple*>(inf->Get(Form("nt_%s_%s",fANucleus.GetName(),fBNucleus.GetName())));
    if (!fNt) {
      cerr << "ReadNextEvent: Could not find ntuple!" << endl;
      inf->ls();
      return 0;
    }
    fNt->SetBranchAddress("Npart",&fEv.Npart);
    fNt->SetBranchAddress("Ncoll",&fEv.Ncoll);
    fNt->SetBranchAddress("B",&fEv.B);
    fNt->SetBranchAddress("BNN",&fEv.BNN);
    fNt->SetBranchAddress("VarX",&fEv.VarX);
    fNt->SetBranchAddress("VarY",&fEv.VarY);
    fNt->SetBranchAddress("VarXY",&fEv.VarXY);
    fNt->SetBranchAddress("NpartA",&fEv.NpartA);
    fNt->SetBranchAddress("NpartB",&fEv.NpartB);
    fNt->SetBranchAddress("Npart0",&fEv.Npart0);
    fNt->SetBranchAddress("Psi1",&fEv.Psi1);
    fNt->SetBranchAddress("Ecc1",&fEv.Ecc1);
    fNt->SetBranchAddress("Psi2",&fEv.Psi2);
    fNt->SetBranchAddress("Ecc2",&fEv.Ecc2);
    fNt->SetBranchAddress("Psi3",&fEv.Psi3);
    fNt->SetBranchAddress("Ecc3",&fEv.Ecc3);
    fNt->SetBranchAddress("Psi4",&fEv.Psi4);
    fNt->SetBranchAddress("Ecc4",&fEv.Ecc4);
    fNt->SetBranchAddress("Psi5",&fEv.Psi5);
    fNt->SetBranchAddress("Ecc5",&fEv.Ecc5);
    return 1;
  }
  if ((!inf)||(!fNt&&!calc)) {
    cerr << "ReadNextEvent was not initialized" <<endl;
    return 0;
  }
  TObjArray *arr = dynamic_cast<TObjArray*>(inf->Get(Form("nucleonarray%d",iev)));
  if (!arr) {
    if (iev==0) {
      cerr << "ReadNextEvent could not read nucleon array for event " << iev << endl;
      return 0;
    }
    iev = 0;
    cerr << "ReadNextEvent resetting to first event" << endl;
    arr = dynamic_cast<TObjArray*>(inf->Get(Form("nucleonarray%d",iev)));
  }

  Double_t bgenA=0, bgenB=0;
  Int_t inA=0, inB=0;
  const Int_t nNucls = arr->GetEntries();
  for (Int_t iNucl=0; iNucl<nNucls; ++iNucl) {
    TGlauNucleon *nuclinp = static_cast<TGlauNucleon*>(arr->At(iNucl));
    TGlauNucleon *nuclout = 0;
    if (nuclinp->IsInNucleusB()) { 
      nuclout = static_cast<TGlauNucleon*>(fNucleonsB->At(inB));
      bgenB += nuclinp->GetX();
      ++inB;
    } else {
      nuclout = static_cast<TGlauNucleon*>(fNucleonsA->At(inA));
      bgenA += nuclinp->GetX();
      ++inA;
    }
    nuclout->Reset();
    nuclout->SetXYZ(nuclinp->GetX(),nuclinp->GetY(),nuclinp->GetZ());
    nuclout->SetType(nuclinp->GetType());
    nuclout->SetEnergy(nuclinp->GetEnergy());
    if (!calc)
      nuclout->SetNColl(nuclinp->GetNColl());
  }
  delete arr;
  Double_t bgen = bgenB/inB-bgenA/inA;
  if (calc) {
    Bool_t ret = CalcEvent(bgen);
    if (0) 
      cout << iev << ": " << fEv.B << " " << fEv.Npart << " " << fEv.Ncoll << " " << fEv.Npart0 << endl;
    ++iev;
    return ret;
  }
  Int_t ret = fNt->GetEntry(iev);
  if (ret<=0) 
    return 0;
  fEccN[1]=fEv.Ecc1;
  fEccN[2]=fEv.Ecc2;
  fEccN[3]=fEv.Ecc3;
  fEccN[4]=fEv.Ecc4;
  fEccN[5]=fEv.Ecc5;
  if (0) 
    cout << iev << ": " << fEv.B << " " << fEv.Npart << " " << fEv.Ncoll << " " << fEv.Npart0 << endl;
  if (0) { // test ntuple values vs re-calculated values
    Double_t npart = fEv.Npart;
    Double_t ncoll = fEv.Ncoll;
    Double_t ecc2  = fEv.Ecc2;
    CalcEvent(bgen);
    if (npart!=fEv.Npart) 
      cout << iev << " differ in npart " << npart << " " << fEv.Npart << endl;
    if (ncoll!=fEv.Ncoll) 
      cout << iev << " differ in ncoll " << ncoll << " " << fEv.Ncoll << endl;
    if (TMath::Abs(ecc2-fEv.Ecc2)>0.001) 
      cout << iev << " differ in ecc2 " << ecc2 << " " << fEv.Ecc2 << endl;
  }
  ++iev;
  return 1;
}

void TGlauberMC::Run(Int_t nevents, Double_t b)
{
  if (fNt == 0) {
    TString name(Form("nt_%s_%s",fANucleus.GetName(),fBNucleus.GetName()));
    TString title(Form("%s + %s (x-sect = %.1f mb) str %s",fANucleus.GetName(),fBNucleus.GetName(),fXSect,Str()));
    TString vars("Npart:Ncoll:Nhard:B:BNN:Ncollpp:Ncollpn:Ncollnn:VarX:VarY:VarXY:NpartA:NpartB:Npart0:AreaW");
    if (fDetail>1)
      vars+=":Psi1:Ecc1:Psi2:Ecc2:Psi3:Ecc3:Psi4:Ecc4:Psi5:Ecc5";
    if (fDetail>2)
      vars+=":AreaO:AreaA:X0:Y0:Phi0:Length";
    if (fDetail>3)
      vars+=":MeanX:MeanY:MeanX2:MeanY2:MeanXY:MeanXSystem:MeanYSystem:MeanXA:MeanYA:MeanXB:MeanYB";
    if (fDetail>4)
      vars+=":PhiA:ThetaA:PhiB:ThetaB";
    fNt = new TNtuple(name,title,vars);
    fNt->SetDirectory(0);
    TObjArray *l = fNt->GetListOfBranches();
    for (Int_t i=0; i<l->GetEntries(); ++i) {
      TBranch *br = dynamic_cast<TBranch*>(l->At(i));
      if (br)
        br->SetCompressionLevel(9);
    }
  }

  for (Int_t i = 0; i<nevents; ++i) {
    while (!NextEvent(b)) {}
    fNt->Fill((Float_t*)(&fEv.Npart));
    if ((i>0)&&(i%100)==0) 
      cout << "Event # " << i << " x-sect = " << GetTotXSect() << " +- " << GetTotXSectErr() << " b        \r" << flush;
  }
  if (nevents>99)
    cout << endl << "Done!" << endl;
}