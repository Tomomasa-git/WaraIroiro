/* * * * * * * * * * * * * *
 *   DrawColorCont.cc      *
 *    2022. 02. 11 (Sat.)  *
 *    T. FUJIWARA          *
 * * * * * * * * * * * * * */
#include <iostream>
#include <fstream>
#include <math.h>
#include <string>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <dirent.h>
#include <unistd.h>
using namespace std;

#include "TApplication.h"
#include "TApplicationImp.h"
#include "TCanvasImp.h"
#include "TH1.h"
#include "TH2.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TF1.h"
#include "TF1Convolution.h"
#include "TFile.h"
#include "TLeaf.h"
#include "TBranch.h"
#include "TTree.h"
#include "TCut.h"
#include "TCutG.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TFrame.h"
#include "TPad.h"
#include "TVector3.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TGraph2D.h"
#include "TGraph2DErrors.h"
#include "TAxis.h"
#include "TGaxis.h"
#include "TPaletteAxis.h"
#include "TLine.h"
#include "TPolyLine.h"
#include "TPolyLine3D.h"
#include "TCurlyLine.h"
#include "TArrow.h"
#include "TEllipse.h"
#include "TArc.h"
#include "TCurlyArc.h"
#include "TBox.h"
#include "TPave.h"
#include "TPaveStats.h"
#include "TPaveText.h"
#include "TStopwatch.h"
#include "TLatex.h"
#include "TString.h"
#include "TText.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TGraphErrors.h"
#include "TProfile.h"
#include "TSystem.h"
#include "TColor.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "TMath.h"
#include "TMinuit.h"
#include "TSpectrum.h"
#include "TRandom.h"

//#include "../include/Dict.cc"
//#include "../include/Define.hh"

static const int ncnv      = 2;
static const int nrun      = 5;
static const int nparticle = 2;
static const int nItr      = 1E+7;

//enum EColor { kWhite =0,   kBlack =1,   kGray=920,
//              kRed   =632, kGreen =416, kBlue=600, kYellow=400, kMagenta=616, kCyan=432,
//              kOrange=800, kSpring=820, kTeal=840, kAzure =860, kViolet =880, kPink=900 };

static const double PI            = 4.0*atan(1.);
static const double deg_to_rad    = PI / 180.;
static const double rad_to_deg    = 180. / PI;
static const double mrad_to_deg   = 1./1000*180./PI;
static const double sigma_to_fwhm = 2.*sqrt(2.*log(2.));
static const double fwhm_to_sigma = 1./sigma_to_fwhm;
static const double cm_to_barn    = 1e+24;

static const double re = 2.817e-13;           // classical electron radius (cm)
static const double Na = 6.02214129*1e+23;    // Avogadro constant

static const double Me   = 0.510998928; // electron     mass (MeV/c2)
static const double Mmu  = 105.6583715; // muon         mass (MeV/c2)
static const double Mpi  = 139.57018;   // charged pion mass (MeV/c2)
static const double Mpi0 = 134.9766;    // neutral pion mass (MeV/c2)
static const double MK   = 493.677;     // charged Kaon mass (MeV/c2)
static const double Mp   = 938.272046;  // proton       mass (MeV/c2)
static const double Mn   = 939.565379;  // neutron      mass (MeV/c2)
static const double Mu   = 931.494061;  // proton       mass (MeV/c2)
static const double ML   = 1115.683;    // Lambda       mass (MeV/c2)
static const double MS0  = 1192.642;    // Sigma Zero   mass (MeV/c2)
static const double MSm  = 1197.449;    // Sigma Minus  mass (MeV/c2)
static const double MSp  = 1189.37;     // Sigma Plus   mass (MeV/c2)
//++++++++++++++++++++++++//
const double OurGoalofNewToF = 150.;   // Ore-tachi no goal of new NKS2 ToF det. (ps)
//++++++++++++++++++++++++//
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
void   SetTH1     ( TH1 *h1, TString hname="", TString xname="", TString yname="", int LColor=602, int FStyle=3005, int FColor=7 );
void   SetTH2     ( TH2 *h2, TString hname="", TString xname="", TString yname="" );
void   SetGraph   ( TGraph *gr      , TString hname="", TString xname="", TString yname="", int LColor=602, int LStyle=2, int LWidth=2, int MColor=6, int MStyle=29, double Yoffset=0.8, double min=0., double max=100. );
void   SetGrErr   ( TGraphErrors *gr, TString hname="", TString xname="", TString yname="", int LColor=602, int LStyle=2, int LWidth=2, int MColor=6, int MStyle=29, double Yoffset=0.8, double min=0., double max=100. );
void   SetTF1     ( TF1 *f, int LColor=2, int LWidth=2, int LStyle=2, int Npx=10000 );
void   SetLatex   ( TLatex *tlat, int tfont=42, int talign=12, int tcol=603, double tsize=0.05 );
void   SetLine    ( TLine *ln, int lcol=602, int lsty=7, int lwid=3 );
void   SetBox     ( TBox *box, int bcol=602, int bsty=7, int bwid=2 );
void   SetLegend  ( TLegend *leg, int Font=42, int Al=22, int Col=603, double Size=0.065, int TranspFlag=1 );
int    SearchFile ( string path, string file );
void   RedrawFrame( TFrame *frame );
double GetFWHM    ( TH1 *h );
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
void SetTH1(TH1 *h1, TString hname, TString xname, TString yname, int LColor, int FStyle, int FColor){
  h1->SetTitle(hname);
  h1->GetXaxis()->SetTitle(xname);
  h1->GetXaxis()->CenterTitle();
  h1->GetYaxis()->SetTitle(yname);
  h1->GetYaxis()->CenterTitle();
  h1->SetMinimum(0.8);
  h1->SetLineColor(LColor);
  h1->SetFillStyle(FStyle);
  h1->SetFillColor(FColor);
  h1->GetYaxis()->SetTitleOffset(1.2);
  h1->SetLineWidth(1);
  h1->SetTitleSize(0.05,"");
  h1->GetXaxis()->SetTitleFont(62);
  h1->GetXaxis()->SetTitleSize(0.08);
  h1->GetXaxis()->SetTitleOffset(0.90);
  h1->GetXaxis()->SetLabelFont(42);
  h1->GetXaxis()->SetLabelSize(0.050);
  h1->GetYaxis()->SetTitleFont(62);
  h1->GetYaxis()->SetTitleSize(0.08);
  h1->GetYaxis()->SetTitleOffset(0.90);
  h1->GetYaxis()->SetLabelFont(42);
  h1->GetYaxis()->SetLabelSize(0.055);
}

////////////////////////////////////////////////////////////////////////////
void SetTH2(TH2 *h2, TString hname, TString xname, TString yname){
  h2->SetTitle(hname);
  h2->GetXaxis()->SetTitle(xname);
  h2->GetXaxis()->CenterTitle();
  h2->GetYaxis()->SetTitle(yname);
  h2->GetYaxis()->CenterTitle();
  h2->SetMinimum(0.8);
  h2->SetLineWidth(0);
  h2->SetTitleSize(0.05,"");
  h2->GetXaxis()->SetTitleFont(62);
  h2->GetYaxis()->SetTitleFont(62);
  h2->GetXaxis()->SetTitleSize(0.08);
  h2->GetYaxis()->SetTitleSize(0.07);
  h2->GetXaxis()->SetTitleOffset(0.90);
  h2->GetYaxis()->SetTitleOffset(0.80);
  h2->GetXaxis()->SetLabelSize(0.060);
  h2->GetYaxis()->SetLabelSize(0.060);
  h2->GetZaxis()->SetLabelSize(0.060);
}
////////////////////////////////////////////////////////////////////////////
void SetGraph(TGraph *gr, TString hname, TString xname, TString yname, int LColor, int LStyle, int LWidth, int MColor, int MStyle, double Yoffset, double min, double max){
  gr->SetTitle(hname);
  gr->SetName(hname);
  gr->GetXaxis()->SetTitle(xname);
  gr->GetXaxis()->SetTitleFont(62);
  gr->GetXaxis()->SetTitleSize(0.080);
  gr->GetXaxis()->SetTitleOffset(0.90);
  gr->GetXaxis()->SetLabelFont(42);
  gr->GetXaxis()->SetLabelSize(0.060);
  gr->GetXaxis()->CenterTitle();
  gr->GetYaxis()->SetTitle(yname);
  gr->GetYaxis()->SetTitleFont(62);
  gr->GetYaxis()->SetTitleSize(0.080);
  gr->GetYaxis()->SetTitleOffset(0.70);
  gr->GetYaxis()->SetLabelFont(42);
  gr->GetYaxis()->SetLabelSize(0.060);
  gr->GetYaxis()->CenterTitle();

  gr->SetLineColor(LColor);
  gr->SetLineStyle(LStyle);
  gr->SetMarkerStyle(MStyle);
  gr->SetMarkerColor(MColor);
  gr->SetMarkerSize(2.);
//  gr->GetYaxis()->SetTitleOffset(Yoffset);
//  gr->GetYaxis()->SetRangeUser(min,max);
}

////////////////////////////////////////////////////////////////////////////
void SetGrErr(TGraphErrors *gr, TString hname, TString xname, TString yname, int LColor, int LStyle, int LWidth, int MColor, int MStyle, double Yoffset, double min, double max){
  gr->SetTitle(hname);
  gr->SetName(hname);
  gr->GetXaxis()->SetTitle(xname);
  gr->GetXaxis()->CenterTitle();
  gr->GetYaxis()->SetTitle(yname);
  gr->GetYaxis()->CenterTitle();
  gr->SetLineColor(LColor);
  gr->SetMarkerStyle(MStyle);
  gr->SetMarkerColor(MColor);
  gr->SetMarkerSize(0.5);
  gr->GetYaxis()->SetTitleOffset(Yoffset);

  gr->SetLineColor(LColor);
  gr->SetLineStyle(LStyle);
  gr->SetMarkerStyle(MStyle);
  gr->SetMarkerColor(MColor);

  gr->SetMarkerSize(0.5);
//  gr->GetYaxis()->SetRangeUser(min,max);
}
////////////////////////////////////////////////////////////////////////////

void SetTF1(TF1 *f, int LColor, int LWidth, int LStyle, int Npx){
  f->SetLineColor(LColor);
  f->SetLineWidth(LWidth);
  f->SetLineStyle(LStyle);
  f->SetNpx(Npx);
}

////////////////////////////////////////////////////////////////////////////
void SetLatex( TLatex *tlat, int tfont, int talign, int tcol, double tsize ){
  tlat->SetTextFont(tfont);
  tlat->SetTextAlign(talign);
  tlat->SetTextColor(tcol);
  tlat->SetTextSize(tsize);
}

////////////////////////////////////////////////////////////////////////////
void SetLine( TLine *ln, int lcol, int lsty, int lwid ){
  ln->SetLineColor(lcol);
  ln->SetLineWidth(lwid);
  ln->SetLineStyle(lsty);
}

////////////////////////////////////////////////////////////////////////////
void SetBox( TBox *box, int bcol, int bsty, int bwid ){
  box->SetLineColor(bcol);
  box->SetLineWidth(bwid);
  box->SetLineStyle(bsty);
  
  box->SetFillColor(0);
  box->SetFillStyle(0);
}

////////////////////////////////////////////////////////////////////////////
void SetLegend( TLegend *leg, int Font, int Al, int Col, double Size, int TranspFlag ){
  leg->SetTextSize(Size);
  leg->SetTextFont(Font);
  leg->SetTextAlign(Al);
  leg->SetTextColor(Col);
  leg->SetBorderSize(4);
  
  
  if(TranspFlag==1){
  	leg->SetLineColor(0);
  	leg->SetLineWidth(0);
  	leg->SetFillColor(0);
  	leg->SetFillStyle(0);
  	leg->SetBorderSize(0);
  }else;
}

////////////////////////////////////////////////////////////////////////////
int SearchFile(string path, string file){
  const char *SearchPath = path.c_str();
  string SearchTarget=file;
  DIR *DirPtr;
  dirent *Entry;
  
  int id=0;
  
  DirPtr = opendir(SearchPath);
  if(DirPtr==NULL)exit(1);
  do{
  	Entry = readdir(DirPtr);
  	if(Entry!=NULL){
  		if( (Entry->d_name)==SearchTarget )id++;
  	}
  		
  }while(Entry!=NULL);
  
  return id;
}		

////////////////////////////////////////////////////////////////////////////
void RedrawFrame(TFrame *frame){
  gPad->Update();
  gPad->Modified(1);
  
  frame=gPad->GetFrame();
  frame->SetFillStyle(0);
  frame->Draw();
  
  gPad->Update();
  gPad->Modified(1);
  
  return;
}//RedrawFrame

///////////////////////
double GetFWHM(TH1 *h){
  double fwhm=0.;
  double MaximumThisHisto=0.;
  double HalfMaximum=0.;
  double Half[2];
  
  // The case that "h" does not exist.
  if(!h){
  	cout<<"There is no object: "<<h->GetName()<<endl;
  	fwhm=-9999.;
  	return fwhm;
  }else;
  
  MaximumThisHisto = h->GetMaximum();	
  HalfMaximum = MaximumThisHisto/2.;
  Half[0] = h->FindFirstBinAbove(HalfMaximum);
  Half[1] = h->FindLastBinAbove(HalfMaximum);
  
  // Always positive value
  fwhm = fabs(Half[1]-Half[0]);	
  
  return fwhm;
}// GetFWHM
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
class MyDraw{
public:
   MyDraw(int ccontnnumber=0);
  ~MyDraw();
  void SetRoot();
  void SetIO();
  void DefineCanv();
  void DefineObj();
  void Fill();
  void DrawHist();
  void Export();

private:
  TCanvas *ca[ncnv];
  TH2D *h_cont[2];
  double gaMean;
  double gaSigma;
  double bwMean;
  double bwGamma;

  TFile *ofp;

  string pdfpath="./pdf/";
  string pdfname;    
  string pdfname_[ncnv];    
  string pdf_op;    
  string pdf_cl;    
  string rootpath="./output/";
  string rootname;    
  string txtpath="./txt/";
  string tfname;

};
//_____________________________________________
MyDraw::MyDraw(int ccontnnumber){

  //////////////////////////////
  // General setting
  gErrorIgnoreLevel = kError;
  gROOT -> Reset();
  
  gROOT->SetBatch(kTRUE);
  
  gStyle->SetOptDate(0);
  gStyle->SetHistFillStyle(3002);
  gStyle->SetHistFillColor(0);
  // Grid
  gStyle->SetGridColor(kGray);
  gStyle->SetGridWidth(1);
  gStyle->SetPadGridX(1);
  gStyle->SetPadGridY(1);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  // Margin
  gStyle->SetPadRightMargin(.050);
  gStyle->SetPadLeftMargin(.160);
  gStyle->SetPadTopMargin(.05);
  gStyle->SetPadBottomMargin(.18);
  gStyle->SetFrameLineWidth(1);
  gStyle->SetLineWidth(1);
  gStyle->SetOptDate(0);
  gStyle->SetOptStat("ei");
  //gStyle->SetStatW(0.15);
  gStyle->SetStatFontSize(0.03);
  gStyle->SetStatTextColor(1);
  gStyle->SetTitleX(0.30);
  gStyle->SetTitleFontSize(0.05);
  gStyle->SetTitleTextColor(1);
  gStyle->SetNdivisions(510); // tertiary*10000 + secondary*100 + first
  
  const Int_t NRGBs = 5;
  const Int_t NCont = 99;
  Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
  Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
  Double_t green[NRGBs] = { 0.00, 0.80, 1.00, 0.20, 0.00 };
  Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
  TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
  gStyle->SetNumberContours(NCont);
  //////////////////////////////

  runnum[0] = runnumber;
  runnum[1] = runnumber+1;
  runnum[2] = runnumber+2;
  runnum[3] = runnumber+3;
  runnum[4] = runnumber+4;
}

//_____________________________________________
MyDraw::~MyDraw(){
}

//_____________________________________________
