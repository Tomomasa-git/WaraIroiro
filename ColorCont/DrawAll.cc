/* * * * * * * * * * * * * *
 *    DrawAll.cc           *
 *    2022. 02. 23 (Wed.)  *
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

static const int nECP = 64;
static const int ncnv = 128;
static const int nItr = 1E+7;

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
//  h2->SetMinimum(0.8);
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
  h2->GetZaxis()->SetLabelSize(0.050);
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
     int ItNum=0;
     bool fIsDebugFlag=kFALSE;
     int fDebugModeLevel=0;
   
     TCanvas *ca[ncnv];
     TH2D *h_cont[nECP][2];
     TFrame *fr_cont[nECP][2];
     double gaMean=0.;
     double gaSigma=2.0;
     double bwMean=0.;
     double bwGamma=3.0;
   
     TFile *ifp;
     TFile *ofp;
   
     string pdfpath="./pdf/Separated/";
     string pdfname;    
     string pdfname_[ncnv];    
     string pdf_op;    
     string pdf_cl;    
     string inputpath="./output/";
     string inputname;    
     string rootpath="./output/Separated/";
     string rootname;    
     string txtpath="./txt/";
     string tfname;
   
     int pltcol[nECP] = {
                          0,                         kDeepSea,                  kGreyScale,            kDarkBodyRadiator,  kBlueYellow,     kRainBow, 
                          kInvertedDarkBodyRadiator, kBird,                 kCubehelix,         kGreenRedViolet, kBlueRedYellow,
                          kOcean,                    kColorPrintableOnGrey, kAlpine,            kAquamarine,     kArmy,
                          kAtlantic,                 kAurora,               kAvocado,           kBeach,          kBlackBody,
                          kBlueGreenYellow,          kBrownCyan,            kCMYK,              kCandy,          kCherry,
                          kCoffee,                   kDarkRainBow,          kDarkTerrain,       kFall,           kFruitPunch,
                          kFuchsia,                  kGreyYellow,           kGreenBrownTerrain, kGreenPink,      kIsland,
                          kLake,                     kLightTemperature,     kLightTerrain,      kMint,           kNeon,
                          kPastel,                   kPearl,                kPigeon,            kPlum,           kRedBlue,
                          kRose,                     kRust,                 kSandyTerrain,      kSienna,         kSolar,
                          kSouthWest,                kStarryNight,          kSunset,            kTemperatureMap, kThermometer,
                          kValentine,                kVisibleSpectrum,      kWaterMelon,        kCool,           kCopper,
                          kGistEarth,                kViridis,              kCividis
                        };
   
   
     string pltcolname[nECP] = {
                                 "kNagaoPalette",             "kDeepSea",              "kGreyScale",         "kDarkBodyRadiator",  "kBlueYellow",     "kRainBow", 
                                 "kInvertedDarkBodyRadiator", "kBird",                 "kCubehelix",         "kGreenRedViolet",    "kBlueRedYellow",
                                 "kOcean",                    "kColorPrintableOnGrey", "kAlpine",            "kAquamarine",        "kArmy",
                                 "kAtlantic",                 "kAurora",               "kAvocado",           "kBeach",             "kBlackBody",
                                 "kBlueGreenYellow",          "kBrownCyan",            "kCMYK",              "kCandy",             "kCherry",
                                 "kCoffee",                   "kDarkRainBow",          "kDarkTerrain",       "kFall",              "kFruitPunch",
                                 "kFuchsia",                  "kGreyYellow",           "kGreenBrownTerrain", "kGreenPink",         "kIsland",
                                 "kLake",                     "kLightTemperature",     "kLightTerrain",      "kMint",              "kNeon",
                                 "kPastel",                   "kPearl",                "kPigeon",            "kPlum",              "kRedBlue",
                                 "kRose",                     "kRust",                 "kSandyTerrain",      "kSienna",            "kSolar",
                                 "kSouthWest",                "kStarryNight",          "kSunset",            "kTemperatureMap",    "kThermometer",
                                 "kValentine",                "kVisibleSpectrum",      "kWaterMelon",        "kCool",              "kCopper",
                                 "kGistEarth",                "kViridis",              "kCividis"
                               };
   
     //enum EColorPalette{
     //                    kDeepSea=51,          kGreyScale=52,    kDarkBodyRadiator=53,
     //                    kBlueYellow= 54,      kRainBow=55,      kInvertedDarkBodyRadiator=56,
     //                    kBird=57,             kCubehelix=58,    kGreenRedViolet=59,
     //                    kBlueRedYellow=60,    kOcean=61,        kColorPrintableOnGrey=62,
     //                    kAlpine=63,           kAquamarine=64,   kArmy=65,
     //                    kAtlantic=66,         kAurora=67,       kAvocado=68,         
     //                    kBeach=69,            kBlackBody=70,    kBlueGreenYellow=71,
     //                    kBrownCyan=72,        kCMYK=73,         kCandy=74,
     //                    kCherry=75,           kCoffee=76,       kDarkRainBow=77,
     //                    kDarkTerrain=78,      kFall=79,         kFruitPunch=80,
     //                    kFuchsia=81,          kGreyYellow=82,   kGreenBrownTerrain=83,
     //                    kGreenPink=84,        kIsland=85,       kLake=86,
     //                    kLightTemperature=87, kLightTerrain=88, kMint=89,
     //                    kNeon=90,             kPastel=91,       kPearl=92,
     //                    kPigeon=93,           kPlum=94,         kRedBlue=95,
     //                    kRose=96,             kRust=97,         kSandyTerrain=98,
     //                    kSienna=99,           kSolar=100,       kSouthWest=101,
     //                    kStarryNight=102,     kSunset=103,      kTemperatureMap=104,
     //                    kThermometer=105,     kValentine=106,   kVisibleSpectrum=107,
     //                    kWaterMelon=108,      kCool=109,        kCopper=110,
     //                    kGistEarth=111,       kViridis=112,     kCividis=113
     //};       
   
   
   
};
//_____________________________________________
//_____________________________________________
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
  gStyle->SetPadRightMargin(.150);
  gStyle->SetPadLeftMargin(.160);
  gStyle->SetPadTopMargin(.050);
  gStyle->SetPadBottomMargin(.16);
  gStyle->SetFrameLineWidth(1);
  gStyle->SetLineWidth(1);
  gStyle->SetOptDate(0);
  gStyle->SetOptStat("ei");
  //gStyle->SetStatW(0.15);
  gStyle->SetStatFontSize(0.03);
  gStyle->SetStatTextColor(1);
  gStyle->SetStatX(0.50);
  gStyle->SetStatY(0.90);
  gStyle->SetTitleX(0.325);
  gStyle->SetTitleFontSize(0.04);
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
}

//_____________________________________________
MyDraw::~MyDraw(){
}

//_____________________________________________
void MyDraw::SetIO(){
  pdfname=Form("DrawColorCont_Separated_%03d.pdf",ItNum);
  while( SearchFile(rootpath,rootname)!=0 ){
    ItNum++;
    rootname = Form("DrawColorCont_Separated_%03d.root",ItNum);
    pdfname  = Form("DrawColorCont_Separated_%03d.pdf" ,ItNum);
  }
  rootname =  Form("DrawColorCont_Separated_%03d.root",ItNum);
  pdfname  =  Form("DrawColorCont_Separated_%03d.pdf" ,ItNum);
  rootname = rootpath+rootname;
  pdfname  = pdfpath+pdfname;
  pdf_op = pdfname + "[";
  pdf_cl = pdfname + "]";

  for(int i=0;i<nECP;i++){
    pdfname_[i]      = pdfpath+Form("DrawColorCont_Separated_gaus_%03d_%03d_%s.pdf",ItNum,i,pltcolname[i].c_str() );
    pdfname_[i+nECP] = pdfpath+Form("DrawColorCont_Separated_bw_%03d_%03d_%s.pdf"  ,ItNum,i,pltcolname[i].c_str() );
  }

  inputname = "DrawColorCont_Single_000.root";
  inputname = inputpath + inputname;

  SetRoot();
}

//_____________________________________________
void MyDraw::SetRoot(){
  ifp = new TFile( inputname.c_str(), "READONLY" );
}

//_____________________________________________

void MyDraw::DefineCanv(){
  for(int i=0; i<ncnv; i++){
    ca[i] = new TCanvas( Form("ca[%d]",i), Form("ca[%d]",i), 1202,1024);
  }// for i
}

//_____________________________________________
void MyDraw::DefineObj(){
  TString htitle="";

  ifp->cd();
  for(int i=0;i<nECP;i++){
    for(int j=0;j<2;j++){
      h_cont[i][j] = (TH2D*)gROOT->FindObject(Form("h_cont_%02d",j));
      if(!h_cont[i][j])cout<<"\033[1;31m"<<"bbbbbbb"<<"\033[m"<<endl;
      h_cont[i][j]->SetStats(kFALSE);
    }// for j

    htitle = pltcolname[i]+" Gaus dist.";  
    SetTH2(h_cont[i][0], htitle, "x", "y");

    htitle = pltcolname[i]+" B.W. dist";   
    SetTH2(h_cont[i][1], htitle, "x", "y");
  }

}

//_____________________________________________
void MyDraw::Fill(){
}

//_____________________________________________
void MyDraw::DrawHist(){
  cout<<endl;
  cout<<"Draw start."<<endl;
  for(int i=0; i<nECP; i++){
    cout<<setw(4)<<i<<"/"<<nECP<<"\r"<<flush;

    ca[i]->cd();
    gPad->SetXstat(0.50);
    gPad->SetYstat(0.90);
    gPad->SetLogz(kTRUE);
    if(i>0)gStyle->SetPalette(pltcol[i]);
    h_cont[i][0]->Draw("colz");
    RedrawFrame(fr_cont[i][0]);
    ca[i]->Print(pdfname_[i].c_str());
    
    ca[i+nECP]->cd();
    gPad->SetLogz(kTRUE);
    if(i>0)gStyle->SetPalette(pltcol[i]);
    h_cont[i][1]->Draw("colz");
    RedrawFrame(fr_cont[i][1]);
    ca[i+nECP]->Print(pdfname_[i+nECP].c_str());
  }//for i

  cout<<endl;
  cout<<"Draw end."<<endl;
}

//_____________________________________________
void MyDraw::Export(){

  cout<<"\033[1;33m"<<endl;
  cout<<"==============="<<endl;
  cout<<"  Info  of MyDraw::Export"<<endl;
  cout<<"    This run data was successfully written in "<<pdfname<<"  "<<endl;
  cout<<"==============="<<endl;
  cout<<"\033[m"<<endl;
  cout<<endl;

  ofp = new TFile( rootname.c_str(), "recreate" );
  ofp->cd();
  for(int i=0; i<ncnv; i++){ca[i]->Write(Form("ca_%02d",i));}
  for(int i=0; i<nECP; i++){
   h_cont[i][0]->Write(Form("h_cont_%02d_00",i));
   h_cont[i][1]->Write(Form("h_cont_%02d_01",i));
  }
  ofp->Close();

  cout<<"\033[1;33m"<<endl;
  cout<<"==============="<<endl;
  cout<<"  Info  of MyDraw::Export"                                <<endl;
  cout<<"    This run data was successfully written in "<<rootname<<"  "<<endl;
  cout<<"==============="<<endl;
  cout<<"\033[m"<<endl;
  cout<<endl;
}

////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
//////////////////////////////////  main  //////////////////////////////////
////////////////////////////////////////////////////////////////////////////
int main(int argc, char** argv){
  int ch;
  string varg;
  string aarg;
  //int anamode=0;
  int vflag=0;
  extern char *optarg;

  MyDraw *ana;
  TStopwatch *Stw;

  string file="tmp.root";
  //int runnum=6;
  //bool debug=kFALSE;
	
  cout<<"Processing..."<<endl;
  while( (ch=getopt(argc,argv,"hdr:f:"))!=-1 ){
    cout<<"---"<<endl;
    switch(ch){

      case 'h':
      cout<<endl;
      cout<<"\033[1;36m============ HELP ============\033[m"<<endl;
      cout<<endl;
      cout<<"\033[1;36m==============================\033[m"<<endl;
      cout<<endl;
      return 0;
      break;
      
      case '?':
      cout<<"unknown option...."<<endl;
      return 0;
      break;
      
      default:
      cout<<"\033[1;31mtype -h to see help!!\033[m"<<endl;
      return 0;
    }
    cout<<"---"<<endl;
  }
  cout<<"bbbbbbbbb"<<endl;

  Stw = new TStopwatch();
  cout<<"Start analysis..."<<endl;
  Stw->Start();

  TApplication theApp("App", &argc, argv);

  ana = new MyDraw();
  ana->SetIO();
  ana->DefineCanv();
  ana->DefineObj();
  ana->Fill();
  ana->DrawHist();
  ana->Export();

  Stw->Stop();
  cout<<"---"<<endl;
  cout<<"End analysis..."<<endl;
  cout<<"---"<<endl;
  cout<<endl;
  cout<<"=============================="<<endl;
  cout<<"\033[1;m  Total Time: "<<Stw->RealTime()<<" sec \033[m"<<endl;
  cout<<"=============================="<<endl;
  cout<<endl;

  gSystem->Exit(-1);
  theApp.Run();

  return 0;

} 
////////////////////////////////////////////////////////////////////////////
