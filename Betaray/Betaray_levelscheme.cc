/* * * * * * * * * * * * * * 
 *  Betaray_levelscheme.cc *
 *    2022. 03. 13 (Sun.)  *
 *    T. FUJIWARA          *
 * * * * * * * * * * * * * */
#include <iostream>
#include <fstream>
#include <math.h>
#include <string>
 //
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
/*
*/

static const int ncnv = 99;

//enum EColor { kWhite =0,   kBlack =1,   kGray=920,
//              kRed   =632, kGreen =416, kBlue=600, kYellow=400, kMagenta=616, kCyan=432,
//              kOrange=800, kSpring=820, kTeal=840, kAzure =860, kViolet =880, kPink=900 };

static const double PI = 4.0*atan(1.);
static const double deg_to_rad  = PI / 180.;
static const double rad_to_deg  = 180. / PI;
static const double mrad_to_deg = 1./1000*180./PI;
static const double sigma_to_fwhm = 2.*sqrt(2.*log(2.));
static const double fwhm_to_sigma = 1./sigma_to_fwhm;
static const double cm_to_barn = 1e+24;

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

const Int_t NRGBs = 5;
const Int_t NCont = 99;

//++++++++++++++++++++++++//
const double OurGoalofNewToF = 150.;   // Ore-tachi no goal of new NKS2 ToF det. (ps)
//++++++++++++++++++++++++//

////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
void SetTH1( TH1 *h1, TString hname, TString xname, TString yname, int LColor=4, int FStyle=3005, int FColor=7 );
void SetTH2( TH2 *h2, TString hname, TString xname, TString yname );
void SetGrErr( TGraphErrors *gr, TString hname, TString xname, TString yname, int tFont, int lFont, int lColor, int lStyle, int lWidth, int mColor, int mStyle, double mSize );
void SetTF1( TF1 *f, int LColor=2, int LWidth=2, int LStyle=2, int Npx=10000);
void SetLatex( TLatex *tlat, int tfont=42, int talign=12, int tcol=603, double tsize=0.05 );
void SetLine( TLine *ln, int lcol=602, int lsty=7, int lwid=3 );
void SetArrow( TArrow *ln, int lcol=602, int lsty=7, int lwid=3 );
void SetBox( TBox *box, int bcol=602, int bsty=7, int bwid=2 );
void SetLegend( TLegend *leg, int Font=42, int Al=22, int Col=603, double Size=0.065, int TranspFlag=0 );

int SearchFile(string path, string file);
void RedrawFrame(TFrame *frame);
double GetFWHM(TH1 *h);
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
  h1->GetXaxis()->SetTitleSize(0.08);
  h1->GetYaxis()->SetTitleSize(0.08);
  h1->GetXaxis()->SetTitleOffset(0.90);
  h1->GetYaxis()->SetTitleOffset(0.90);
  h1->GetXaxis()->SetLabelSize(0.0);
  h1->GetYaxis()->SetLabelSize(0.0);
//  h1->GetXaxis()->SetLabelSize(0.050);
//  h1->GetYaxis()->SetLabelSize(0.055);
//
  h1->SetStats(kFALSE);
  h1->GetXaxis()->SetTickSize(0.);
  h1->GetYaxis()->SetTickSize(0.);
  h1->GetXaxis()->SetTickLength(0.);
  h1->GetYaxis()->SetTickLength(0.);
}

////////////////////////////////////////////////////////////////////////////
void SetTH2(TH2 *h2, TString hname, TString xname, TString yname){
  h2->SetTitle(hname);
  h2->GetXaxis()->SetTitle(xname);
  h2->GetXaxis()->CenterTitle();
  h2->GetYaxis()->SetTitle(yname);
  h2->GetYaxis()->CenterTitle();
  h2->SetMinimum(0.8);
  h2->GetYaxis()->SetTitleOffset(1.2);
  h2->SetLineWidth(0);
  h2->SetTitleSize(0.05,"");
  h2->GetXaxis()->SetTitleSize(0.08);
  h2->GetYaxis()->SetTitleSize(0.08);
  h2->GetXaxis()->SetTitleOffset(0.90);
  h2->GetYaxis()->SetTitleOffset(0.90);
  h2->GetXaxis()->SetLabelSize(0.045);
  h2->GetYaxis()->SetLabelSize(0.045);
  h2->GetZaxis()->SetLabelSize(0.040);
}

////////////////////////////////////////////////////////////////////////////
void SetGrErr(TGraphErrors *gr, TString hname, TString xname, TString yname, int tFont, int lFont, int lColor, int lStyle, int lWidth, int mColor, int mStyle, double mSize  ){
  gr->SetTitle(hname);
  gr->SetName(hname);
  gr->GetXaxis()->SetTitleSize(0.08);
  gr->GetXaxis()->SetTitleOffset(0.90);
  gr->GetXaxis()->SetLabelSize(0.050);
  gr->GetXaxis()->SetTitle(xname);
  gr->GetXaxis()->CenterTitle();

  gr->GetYaxis()->SetTitleSize(0.08);
  gr->GetYaxis()->SetTitleOffset(0.90);
  gr->GetYaxis()->SetLabelSize(0.050);
  gr->GetYaxis()->SetTitle(yname);
  gr->GetYaxis()->CenterTitle();

  gr->SetLineColor(lColor);
  gr->SetLineStyle(lStyle);
  gr->SetLineWidth(lWidth);
  gr->SetMarkerColor(mColor);
  gr->SetMarkerStyle(mStyle);
  gr->SetMarkerSize(mSize);
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
void SetArrow( TArrow *ln, int lcol, int lsty, int lwid ){
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
}// SetBox

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
}// SetLegend

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
}// SearchFile		

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
class MyDraw{
  public:
     MyDraw();
    ~MyDraw();
    void SetRoot();
    void SetIO();
    void DefineCanv();
    void DefineObj();
    void DrawObj();
    void Export();
  
  private:
    int ItNum=0;

    TChain *tree;
    TCanvas *ca; 
    
    TFile *ifp;
    TFile *ifp_;
    TH1D *h_fr;
    TFrame *fr;

    TLine *ln_EneLevel[3];
    TLine *ln_base;
    TBox *b_base;
    TArrow *ArTrans[2];

    TLatex *Lat_nucl;
    TLatex *Lat;
    TLegend *Leg;
    
    TFile *ofp;
    
    ofstream ofs;
    
    // Path and filename
    string inputname;
    string inputname_;
    string pdfpath="./pdf/";
    string pdfname;    
    string pdfname_[ncnv];    
    string pdf_op;    
    string pdf_cl;    
  
    // tree branch
    int Ent;
  
    double HalfLife[2] = {28.79, 64.10};
    string HalfLife_str[2];
    
    double Qval[2] = {0.546,2.28};
    double Elevel[3];
    double Elevel_[3];
    double Ratio[2] = {100.,99.9885};

double LnRSize[3] = {0.50, 0.75, 0.75};
double LnLSize[3] = {0.75, 0.75, 0.75};

double ArStart[2] = {0.30, 0.50};
double ArEnd[2]   = {0.30, 0.50};
    double ArDelta = 0.03;
    int ArCol[2] = {603, 420};
    int ArSty = 1;
    int ArWid = 3;
    double ArAng = 40.;
    double ArSize = 0.040;
    
    int lc[3]={602};
    int fc[3]={7};
    int fs[3]={3445};
    
    bool fIsLogY = kFALSE;
    bool fIsLogZ = kFALSE;
    
    bool fIsDebugMode = kTRUE;
};
//////////////////////////////////////////////////
MyDraw::MyDraw(){
  //////////////////////////////
  // General setting
  gErrorIgnoreLevel = kError;
  gROOT->Reset();
  
  gROOT->SetBatch(kTRUE);
  
  gStyle->SetOptDate(0);
  gStyle->SetHistFillStyle(3002);
  gStyle->SetHistFillColor(0);

  // Grid
  gStyle->SetGridWidth(1);
  gStyle->SetPadGridX(0);
  gStyle->SetPadGridY(0);
  gStyle->SetPadTickX(0);
  gStyle->SetPadTickY(0);

  // Margin
  gStyle->SetPadRightMargin(.0020);
  gStyle->SetPadLeftMargin(.0010);
  gStyle->SetPadTopMargin(.070);
  gStyle->SetPadBottomMargin(.020);
  gStyle->SetFrameLineWidth(0);
  gStyle->SetFrameLineColor(0);
  gStyle->SetLineWidth(1);
  gStyle->SetOptDate(0);
  gStyle->SetOptStat("ei");
  gStyle->SetStatFontSize(0.03);
  gStyle->SetStatTextColor(1);
  gStyle->SetTitleX(0.15);
  gStyle->SetTitleFontSize(0.05);
  gStyle->SetTitleTextColor(1);
  gStyle->SetNdivisions(510); // tertiary*10000 + secondary*100 + first
  
  Double_t cStops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
  Double_t cRed[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.40 };
  Double_t cGreen[NRGBs] = { 0.00, 0.80, 1.00, 0.20, 0.00 };
  Double_t cBlue[NRGBs]  = { 0.40, 1.00, 0.12, 0.00, 0.00 };
  float mycol[4][3] = { 
                       { .800, .900, 1.   },
  					   { 1.  , .800, .975 },
  					   { .800, 1.  , .925 }, 
  					   { 1.  , 1.  , .950 } 
                      };
  TColor::CreateGradientColorTable(NRGBs, cStops, cRed, cGreen, cBlue, NCont);
  gStyle->SetNumberContours(NCont);

  for( int n=0;n<4;n++ ){gROOT->GetColor(11+n)->SetRGB( mycol[n][0], mycol[n][1], mycol[n][2] );}
  //////////////////////////////
      
}// MyDraw::MyDraw

//_______________________________________________
MyDraw::~MyDraw(){
}// ~MyDraw::MyDraw

//_______________________________________________
void MyDraw::SetRoot(){
}// MyDraw::SetRoot

//_______________________________________________
void MyDraw::SetIO(){
  pdfname=Form("Betaray_levelscheme_%03d.pdf",ItNum);
  pdfname=pdfpath+pdfname;
  pdf_op=pdfname+"[";
  pdf_cl=pdfname+"]";
}// MyDraw::SetIO

//_______________________________________________
void MyDraw::DefineCanv(){
ca = new TCanvas("ca","ca", 1202, 924);
}// MyDraw::DefineCanv

//_______________________________________________
void MyDraw::DefineObj(){
  Elevel[0] = 0.;
  Elevel[1] = Elevel[0]+Qval[1];
  Elevel[2] = Elevel[1]+Qval[0];

  HalfLife_str[0] = Form("T_{1/2} = %.2lf yr",HalfLife[0]);  
  HalfLife_str[1] = Form("T_{1/2} = %.2lf h",HalfLife[1]);

  // Frame
  h_fr = new TH1D("h_fr","h_fr",1,0,3.8);
  SetTH1(h_fr,"","","");
  h_fr->GetYaxis()->SetRangeUser( 0.,1.10*(Qval[0]+Qval[1]) );
  h_fr->GetXaxis()->SetAxisColor(0);
  h_fr->GetYaxis()->SetAxisColor(0);
  
  // Energy Level
  b_base = new TBox(0,-1.,4.,1.);
  SetBox(b_base, 0, 1, 1);
  b_base->SetFillStyle(1001);
  b_base->SetFillColor(0);

  ln_base = new TLine(0.,0.,4.,0.);
  SetLine(ln_base,0,1,3);

  for(Int_t i_ln=0;i_ln<3;i_ln++){
    //ln_EneLevel[i_ln] = new TLine( ( (Double_t)i_ln+1. )-0.75, Elevel[2-i_ln], ( (Double_t)i_ln+1. )+0.6, Elevel[2-i_ln]);
    ln_EneLevel[i_ln] = new TLine( ( (Double_t)i_ln+1. )-LnLSize[i_ln], Elevel[2-i_ln], ( (Double_t)i_ln+1. )+LnRSize[i_ln], Elevel[2-i_ln]);
    SetLine(ln_EneLevel[i_ln],1,1,3);
  }// for i_ln

  // Arrow
  for(Int_t i_ar=0;i_ar<2;i_ar++){
    //ArTrans[i_ar] = new TArrow( ((Double_t)i_ar+1.)+0.30, Elevel[2-i_ar]-ArDelta,
    //                            ((Double_t)i_ar+1.)+0.30, Elevel[1-i_ar]+2.*ArDelta);
    ArTrans[i_ar] = new TArrow( ((Double_t)i_ar+1.)+ArStart[i_ar], Elevel[2-i_ar]-ArDelta,
                                ((Double_t)i_ar+1.)+ArEnd[i_ar], Elevel[1-i_ar]+2.*ArDelta);
    ArTrans[i_ar]->SetLineWidth(ArWid);
    ArTrans[i_ar]->SetLineStyle(ArSty);
    ArTrans[i_ar]->SetLineColor(ArCol[i_ar]);
    ArTrans[i_ar]->SetFillColor(ArCol[i_ar]);
    ArTrans[i_ar]->SetAngle(ArAng);
    ArTrans[i_ar]->SetArrowSize(ArSize);
  }
  
  // Text
  Lat      = new TLatex();
  Lat_nucl = new TLatex();
  SetLatex(Lat,     42,32,600,0.065);
  SetLatex(Lat_nucl,62,12,602,0.120);

}// MyDraw::DefineObj

//_______________________________________________
void MyDraw::DrawObj(){
  ca->cd();
  h_fr->Draw();
  RedrawFrame(fr);
  
  //ln_base->Draw();
  b_base->Draw();
  for(Int_t i_ln=0;i_ln<3;i_ln++){ln_EneLevel[i_ln]->Draw();}

  // Latex
                         Lat->DrawLatex(1.20,Elevel[2]-0.150, HalfLife_str[0].c_str() );
  Lat->SetTextFont(42);  Lat->DrawLatex(1.20,Elevel[2]-0.390, Form("Q = %.3lf MeV",Qval[0]));
  Lat->SetTextFont(42);  Lat->DrawLatex(1.20,Elevel[2]-0.630, Form("%.0lf %%",Ratio[0]));
  
  Lat->SetTextColor(417);
                         Lat->DrawLatex(2.45,Elevel[1]-0.150, HalfLife_str[1].c_str() );
  Lat->SetTextFont(42);  Lat->DrawLatex(2.45,Elevel[1]-0.390, Form("Q = %.2lf MeV",Qval[1]));
  Lat->SetTextFont(42);  Lat->DrawLatex(2.45,Elevel[1]-0.630, Form("%.4lf %%",Ratio[1]));

  // Nuclei
                                Lat_nucl->DrawLatex(0.20, Elevel[2]+0.26, "{}^{90}_{38}Sr");
  Lat_nucl->SetTextColor(419);  Lat_nucl->DrawLatex(1.55, Elevel[1]+0.26, "{}^{90}_{39}Y");
  Lat_nucl->SetTextColor(1);    Lat_nucl->DrawLatex(2.50, Elevel[0]+0.26, "{}^{90}_{40}Zr");

  // spin-parity
  SetLatex(Lat_nucl,42,12,1,0.06);
  Lat_nucl->DrawLatex(0.90, Elevel[2]+0.15, "J^{#pi}_{G.S.} = 0^{+}");
  Lat_nucl->DrawLatex(2.15, Elevel[1]+0.15, "J^{#pi}_{G.S.} = 2^{-}");
  Lat_nucl->DrawLatex(3.20, Elevel[0]+0.15, "J^{#pi}_{G.S.} = 0^{+}");

  // Arrow
  for(Int_t i_ar=0;i_ar<2;i_ar++){
    ArTrans[i_ar]->Draw("|>");
  }

}// MyDraw::DrawObj

//_______________________________________________
void MyDraw::Export(){
  ca->cd();
  ca->Print( pdf_op.c_str(), "pdf" );
  ca->Print( pdfname.c_str(), "pdf" );
  ca->Print( pdf_cl.c_str(), "pdf" );
  if(fIsDebugMode){
    cout<<"\033[1;32m"<<endl;
    cout<<"==============="<<endl;
    cout<<"  Info  of MyDraw::Export"<<endl;
    cout<<"    This run data was successfully written in "<<pdfname<<"  "<<endl;
    cout<<"    MyDraw::Export was end."<<endl;
    cout<<"==============="<<endl;
    cout<<"\033[m"<<endl;
    cout<<endl;
  }else{
  }
}// MyDraw::Export

//_______________________________________________
/////////////////////////////////////////////////
///////////////////// main //////////////////////
/////////////////////////////////////////////////
int main(int argc, char** argv){
  int ch;
  string varg;
  string aarg;
  //int anamode=0;
  //int vflag=0;
  //extern char *optarg;

  MyDraw *ana;

  cout<<"Processing..."<<endl;
  while( (ch=getopt(argc,argv,"h"))!=-1 ){
    cout<<"---"<<endl;
    switch(ch){
    //  case 't':
    //  break;

      case 'h':
        cout<<endl;
        cout<<"\033[1;33m+---- HELP -----+\033[m"<<endl;
        cout<<"\033[1;33m+---------------+\033[m"<<endl;
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

  TApplication theApp("App", &argc, argv);

  ana = new MyDraw();
  ana->SetIO();
  ana->DefineCanv();
  ana->DefineObj();
  ana->DrawObj();
  ana->Export();

  gSystem->Exit(-1);
  theApp.Run();

  return 0;

} 
