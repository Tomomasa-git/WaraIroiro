/* * * * * * * * * * * * 
 *  BaryonOctet.cc     *
 *                     *
 *  T. Fujiwara        *
 *  2020. 12. 25 (Fri) *
 * * * * * * * * * * * */
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cmath>
using namespace std;

#include "TApplication.h"
#include "TApplicationImp.h"
#include "TCanvasImp.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TRint.h"
#include "TRandom.h"
#include "TStyle.h"
#include "TObject.h"
#include "TList.h"
#include "TFile.h"
#include "TTree.h"
#include "TColor.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TFrame.h"
#include "TCut.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TF1.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TAxis.h"
#include "TGaxis.h"
#include "TMath.h"
#include "TSpectrum.h"
#include "TString.h"
#include "TMathText.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TArrow.h"
#include "TArc.h"
#include "TBox.h"
#include "TPaveStats.h"
#include "TPaveText.h"

#include "./Setting.h"

int main( int argc, char** argv ){
	int ItNum=2;

	const int NofBaryon = 8;
	double Xpos[NofBaryon] = { -0.5, 0.5, -1.0, -0.100,  0.100, 1.0, -0.5,  0.5 };
	double Ypos[NofBaryon] = {  1.0, 1.0,  0.0,  0.200, -0.200, 0.0, -1.0, -1.0 };
	double Xpos_Correct[NofBaryon] = { 0., 0.,  0.   , 0.   , 0., 0.010,  0.025, 0.030 };
	double Ypos_Correct[NofBaryon] = { 0., 0., -0.020, 0.015, 0., 0.   , -0.020, 0.020 };
	double LimX[2] = {-2.0, 2.0};
	double LimY[2] = {-2.3, 1.7};
	double ArrLim[2][4];
	double ArcRad = 0.100*sqrt(5.);
	double LnX1[10] = { -0.50, -0.50, -1.60, -0.75,  0.25, -1.00, -0.50,  0.50, -0.50,  0.50 };
	double LnX2[10] = {  0.50,  0.50,  0.00,  1.00,  1.75, -0.50,  0.50,  1.00, -0.50,  0.50 };
	double LnY1[10] = {  1.00, -1.00,  1.20,  1.50,  1.50,  0.00, -1.00, -1.00,  1.00,  1.00 };
	double LnY2[10] = {  1.00, -1.00, -2.00, -2.00, -1.50,  1.00,  1.00,  0.00, -1.00, -1.00 };
//	TString Name[NofBaryon] = { "n", "p", "#Sigma^{-}", "#Lambda", "#Sigma^{0}", "#Sigma^{+}", "#Xi^{-}", "#Xi^{0}" };
	TString Name[NofBaryon] = { "n", "p", "\\Sigma^{-}", "\\Lambda", "\\Sigma^{0}", "\\Sigma^{+}", "\\Xi^{-}", "\\Xi^{0}" };
	
	double TickLength = 0.060;


	TApplication *theApp;

	Setting* Set;
	TCanvas* Ca;
	TH1D* h_frame;
	TLatex* Lat;
	TLatex* AxText;
	TMathText* Mat;
	TArc* Arc[NofBaryon];
	TArrow* Arr[2];
	TLine* Line[10];
	TLine* Ticks[4];
	TPaveText* Pt;

	theApp = new TApplication( "App", &argc, argv );
	Set = new Setting();
	Set -> Setting_Gene();
	gStyle->SetFrameLineColor(0);
	gStyle->SetPadGridX(0);
	gStyle->SetPadGridY(0);
	gStyle->SetPadTickX(0);
	gStyle->SetPadTickY(1);

	Ca = new TCanvas( "Ca", "Ca", 802, 824 );
	Ca->SetMargin( .025, .025, .010, .040 );
	
//	Lat    = new TLatex();	Set->Setting_Latex( Lat   , 42, 22 , 1, .0725 );
	Lat    = new TLatex();	Set->Setting_Latex( Lat   , 42, 22 , 1, .0725 );
	AxText = new TLatex();	Set->Setting_Latex( AxText, 42, 22 , 1, .060  );
	Mat    = new TMathText();
	Mat->SetTextSize(0.0725);
	Mat->SetTextAlign(22);
	Mat->SetTextColor(1);
	h_frame = new TH1D( "h_frame", "h_frame", 1, LimX[0], LimX[1] );
	Set->Setting_Hist1D( h_frame, "", "", "", 0, 0, 42, 0, 0 );
	h_frame->SetStats(kFALSE);
	h_frame->GetXaxis()->SetNdivisions(0);
	h_frame->GetYaxis()->SetNdivisions(0);
	h_frame->GetXaxis()->SetAxisColor(0);
	h_frame->GetYaxis()->SetAxisColor(0);
	h_frame->GetYaxis()->SetRangeUser(LimY[0], LimY[1]);

	Arr[0] = new TArrow( LimX[0]+0.25, 0.          , LimX[1]-0.25, 0.     , 0.040, "|>");
	Arr[1] = new TArrow( 0.          , LimY[0]+0.25, 0.          , LimY[1], 0.040, "|>");
	for(int i=0; i<2; i++){
		Arr[i]->SetLineWidth(3);
		Arr[i]->SetAngle(40.);
		ArrLim[i][0] = Arr[i]->GetX1();
		ArrLim[i][1] = Arr[i]->GetX2();
		ArrLim[i][2] = Arr[i]->GetY1();
		ArrLim[i][3] = Arr[i]->GetY2();
	}

	for(int i=0; i<NofBaryon; i++){
		Arc[i] = new TArc( Xpos[i], Ypos[i], ArcRad );
		Arc[i]->SetFillColor(11);
		Arc[i]->SetLineColor(4);
		Arc[i]->SetLineWidth(3);
	}

	for(int i=0; i<10; i++){
		Line[i] = new TLine( LnX1[i],LnY1[i], LnX2[i], LnY2[i] );
		Set->Setting_Line( Line[i], 603, 1, 2 );
	//	Set->Setting_Line( Line[i], 602, 1, 2 );
	}
//	Line[2]->SetLineColor(1);
//	Line[2]->SetLineStyle(2);
//	Line[3]->SetLineColor(1);
//	Line[3]->SetLineStyle(2);
//	Line[4]->SetLineColor(1);
//	Line[4]->SetLineStyle(2);

	Ticks[0] = new TLine( -0.5          ,  0.-TickLength, -0.5          ,  0.+TickLength );
	Ticks[1] = new TLine(  0.5          ,  0.-TickLength,  0.5          ,  0.+TickLength );
	Ticks[2] = new TLine(  0.-TickLength,  1.           ,  0.+TickLength,  1.            );
	Ticks[3] = new TLine(  0.-TickLength, -1.           ,  0.+TickLength, -1.            );
	for(int i=0; i<4; i++){Set->Setting_Line( Ticks[i], 1, 3, 1);}


	Pt = new TPaveText( -2.00, LimY[1]-0.40, -1.10, LimY[1]+0.10);
	Pt->SetLineColor(4);
	Pt->SetLineWidth(2);
	Pt->SetFillColor(11);
	Pt->SetTextFont(52);
	Pt->SetTextSize(0.0600);
	Pt->AddText("J^{#it{#pi}}= #frac{1}{2}^{+}");

	h_frame->Draw("");
	
	for(int i=0; i<2; i++){
		Arr[i]->Draw();
	}
	
	for(int i=0; i<10; i++){Line[i]->Draw("");}

	for(int i=0; i<4; i++){Ticks[i]->Draw("");}
	
	for(int i=0; i<NofBaryon; i++){
		Arc[i]->Draw();
		Lat->DrawLatex( Xpos[i]+Xpos_Correct[i], Ypos[i]+Ypos_Correct[i], Name[i] );
	//	Mat->DrawMathText( Xpos[i]+Xpos_Correct[i], Ypos[i]+Ypos_Correct[i], Name[i] );
	}
	
	AxText->SetTextFont(52);
	AxText->DrawLatex(  0.20        , LimY[1]-0.05, "Y"    );
	AxText->DrawLatex(  LimX[1]-0.35, 0.20   , "I#lower[0.25]{#kern[-0.25]{#scale[0.75]{z}}}");
	AxText->SetTextColor(602);
	AxText->SetTextSize(.050);
	AxText->DrawLatex(  LnX2[2]-0.35, LnY2[2]+0.10, "Q=-1" );
	AxText->DrawLatex(  LnX2[3]-0.35, LnY2[3]+0.10, "Q=0" );
	AxText->DrawLatex(  LnX2[4]-0.35, LnY2[4]+0.10, "Q=#lower[-0.18]{+}1" );
	
	AxText->SetTextFont(42);
	AxText->SetTextColor(1);
	AxText->SetTextSize(.050);
	AxText->DrawLatex(  0.22,  0.125       , "O"            );
	AxText->DrawLatex( -1.00,  ArcRad+0.125, "-1"          );
	AxText->DrawLatex(  1.00,  ArcRad+0.125, "#lower[-0.18]{+}1"          );
	AxText->DrawLatex(  0.12,  1.10        , "#lower[-0.18]{+}1"           );
	AxText->DrawLatex(  0.12, -0.90        , "-1"           );
	AxText->SetTextSize(.040);
	AxText->DrawLatex( -0.50,  ArcRad+0.125, "- #frac{#kern[-0.35]{1}}{2}");
	AxText->DrawLatex(  0.50,  ArcRad+0.125, "+ #frac{#kern[-0.35]{1}}{2}");
	
	Pt->Draw();

	Ca->Print(Form("../fig/BaryonOctet_%02d.pdf", ItNum));
//	Ca->Print(Form("../fig/BaryonOctet_%02d.eps", ItNum));
	delete Set;

	gSystem->Exit(-1);

	theApp->Run();

	return 0;
}
