/* * * * * * * * * * * * 
 *  BaryonDecaplet.cc  *
 *                     *
 *  T. Fujiwara        *
 *  2020. 12. 30 (Wed) *
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
	int ItNum=1;

	const int NofBaryon = 10;
	double Xpos[NofBaryon] = { -1.5, -0.5, 0.5, 1.5, -1.0, 0., 1.0, -0.5,  0.5,  0.0 };
	double Ypos[NofBaryon] = {  1.0,  1.0, 1.0, 1.0,  0.0, 0., 0.0, -1.0, -1.0, -2.0 };
	double Xpos_Correct[NofBaryon] = { 0., 0.010, 0.010, 0.040,  0.    , 0.   , 0.   , 0.020,  0.025, 0. };
	double Ypos_Correct[NofBaryon] = { 0., 0.010, 0.010, 0.   ,  0.010 , 0.010, 0.010, 0.010,  0.010, 0. };
	double LimX[2] = {-2.0, 2.0};
	double LimY[2] = {-2.3, 1.7};
	double ArrLim[2][4];
	double ArcRad = 0.100*sqrt(5.);
	double LnX1[12] = { -1.50, -0.50, -1.50, -0.50,  0.50, -1.00, -0.50,  0.00, -0.50,  0.50, -1.50, 1.50 };
	double LnX2[12] = {  1.50,  0.50,  0.00,  0.50,  1.00, -0.50,  0.50,  1.50, -0.50,  0.50, -1.50, 1.50 };
	double LnY1[12] = {  1.00, -1.00,  1.00,  1.00,  1.00,  0.00, -1.00, -2.00,  1.00,  1.00,  1.00, 1.00 };
	double LnY2[12] = {  1.00, -1.00, -2.00, -1.00,  0.00,  1.00,  1.00,  1.00, -1.00, -1.00,  0.00, 0.00 };
	TString Name[NofBaryon] = { 
	                           "#Delta^{-}", "#Delta^{0}", "#Delta^{+}", "#Delta^{++}", 
	                           "#Sigma^{*-}", "#Sigma^{*0}", "#Sigma^{*+}",
	                           "#Xi^{*-}", "#Xi^{*0}",
	                           "#Omega" 
	                          };
	
	double TickLength = 0.060;


	TApplication *theApp;

	Setting* Set;
	TCanvas* Ca;
	TH1D* h_frame;
	TLatex* Lat;
	TLatex* AxText;
	TArc* Arc[NofBaryon];
	TArrow* Arr[2];
	TLine* Line[12];
	TLine* Ticks[6];
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
	
	Lat    = new TLatex();	Set->Setting_Latex( Lat   , 42, 22 , 1, .0725 );
	AxText = new TLatex();	Set->Setting_Latex( AxText, 42, 22, 1, .0600 );

	h_frame = new TH1D( "h_frame", "h_frame", 1, LimX[0], LimX[1] );
	Set->Setting_Hist1D( h_frame, "", "", "", 0, 0, 42, 0, 0 );
	h_frame->SetStats(kFALSE);
	h_frame->GetXaxis()->SetNdivisions(0);
	h_frame->GetYaxis()->SetNdivisions(0);
	h_frame->GetXaxis()->SetAxisColor(0);
	h_frame->GetYaxis()->SetAxisColor(0);
	h_frame->GetYaxis()->SetRangeUser(LimY[0], LimY[1]);

	Arr[0] = new TArrow( LimX[0], 0.     , LimX[1], 0.     , 0.040, "|>");
	Arr[1] = new TArrow( 0.     , LimY[0], 0.     , LimY[1], 0.040, "|>");
	for(int i=0; i<2; i++){
		Arr[i]->SetLineWidth(2);
		Arr[i]->SetAngle(40.);
		ArrLim[i][0] = Arr[i]->GetX1();
		ArrLim[i][1] = Arr[i]->GetX2();
		ArrLim[i][2] = Arr[i]->GetY1();
		ArrLim[i][3] = Arr[i]->GetY2();
	}

	for(int i=0; i<NofBaryon; i++){
		Arc[i] = new TArc( Xpos[i], Ypos[i], ArcRad );
		Arc[i]->SetFillColor(12);
		Arc[i]->SetLineColor(6);
		Arc[i]->SetLineWidth(3);
	}

	for(int i=0; i<12; i++){
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
	Ticks[4] = new TLine( -1.5          ,  0.-TickLength, -1.5          ,  0.+TickLength );
	Ticks[5] = new TLine(  1.5          ,  0.-TickLength,  1.5          ,  0.+TickLength );
	for(int i=0; i<6; i++){Set->Setting_Line( Ticks[i], 1, 3, 1);}


	Pt = new TPaveText( -2.00, LimY[1]-0.40, -1.10, LimY[1]+0.10);
	Pt->SetLineColor(6);
	Pt->SetLineWidth(2);
	Pt->SetFillColor(12);
	Pt->SetTextFont(52);
	Pt->SetTextSize(0.0600);
	Pt->AddText("J^{#it{#pi}}= #frac{3}{2}^{+}");

	h_frame->Draw("");
	
	for(int i=0; i<2; i++){
		Arr[i]->Draw();
	}
	
	for(int i=0; i<12; i++){Line[i]->Draw("");}

	for(int i=0; i<6; i++){Ticks[i]->Draw("");}
	
	for(int i=0; i<NofBaryon; i++){
		Arc[i]->Draw();
		Lat->DrawLatex( Xpos[i]+Xpos_Correct[i], Ypos[i]+Ypos_Correct[i], Name[i] );
	}
	
	AxText->SetTextFont(52);
	AxText->DrawLatex(  0.20        , LimY[1]-0.05, "Y"    );
	AxText->DrawLatex(  LimX[1]-0.15, 0.20   , "I#lower[0.25]{#kern[-0.25]{#scale[0.75]{z}}}");
//	AxText->DrawLatex(  LnX2[2]-0.50, LnY2[2]+0.20, "Q=-1" );
//	AxText->DrawLatex(  LnX2[3]     , LnY2[3]-0.10, "Q= 0" );
//	AxText->DrawLatex(  LnX2[4]     , LnY2[4]-0.10, "Q=+1" );
	
	AxText->SetTextFont(42);
	AxText->SetTextSize(.050);
	AxText->DrawLatex(  0.30,  0.125       , "O"                 );
	AxText->DrawLatex( -1.00,  ArcRad+0.125, "-1"                );
	AxText->DrawLatex(  1.00,  ArcRad+0.125, "#lower[-0.18]{+}1" );
	AxText->DrawLatex(  0.12,  1.10        , "#lower[-0.18]{+}1" );
	AxText->DrawLatex(  0.12, -0.90        , "-1"                );
	AxText->DrawLatex(  0.40, -2.00        , "-2"                );
	AxText->SetTextSize(.040);
	AxText->DrawLatex( -0.50,  ArcRad+0.125, "- #frac{#kern[-0.35]{1}}{2}");
	AxText->DrawLatex(  0.50,  ArcRad+0.125, "+ #frac{#kern[-0.35]{1}}{2}");
	AxText->DrawLatex( -1.50,  ArcRad+0.125, "- #frac{3}{2}");
	AxText->DrawLatex(  1.50,  ArcRad+0.125, "+ #frac{3}{2}");
	
	Pt->Draw();

	Ca->Print(Form("../fig/BaryonDecaplet_%02d.pdf", ItNum));
	delete Set;

	gSystem->Exit(-1);

	theApp->Run();

	return 0;
}
