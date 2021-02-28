/* * * * * * * * * * * * 
 *  main.cc            *
 *                     *
 *  T. Fujiwara        *
 *  2020. 11. 27 (Fri) *
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
#include "TBox.h"
#include "TPaveStats.h"
#include "TPaveText.h"

#include "./HyperTritonLifetime.hh"

#include "./Setting.h"


int main( int argc, char** argv){
	const int NofData  = 11;
	const int NofDummy = 5;
	int ItNum = 3;

	TApplication *theApp;

	TCanvas* Ca;
	TH1D* h_frame;
	TFrame* Fr;
	HyperTritonLifetime* HtlLT[NofData];
	Setting* Set;
	TLegend* Leg;
	TLine*  Ln[3];
		//0: Free Lambda
		//1: Kamada et. al., Phys. Rev. C57 (1998) 1595
		//2: Gal, Garcilazo, Phys. Lett. B791 (2019) 48-53
	TGraphErrors* gr_dummy[NofDummy];

	string figure;
	string figure_op;
	string figure_cl;

	double HMin = 0.;
	double HMax = 500.;
	double WMin = -0.25;
	double WMax = 13.5;

	double LT_Free=263.2;
	double LT_Val[3];
	int LCol[3] = { 602, 6, 870 };
	int LSty[3] = { 7  , 5, 5   };

	int DummyMCol[NofDummy] = { 4 , 3 , 2 , 2 , 2  };
	int DummyMSty[NofDummy] = { 23, 20, 29, 33, 34 };

	LT_Val[0] = LT_Free;
	LT_Val[1] = 256.;
	LT_Val[2] = 213.;	//

	figure = Form("../fig/HTL_Lifetime_%02d.pdf", ItNum);
	figure_op = figure+"[";
	figure_cl = figure+"]";


	theApp = new TApplication( "App", &argc, argv );

	Set = new Setting();
	Set -> Setting_Gene();
	gStyle->SetPadGridX(0);
	gStyle->SetPadGridY(0);
	gStyle->SetPadTickX(0);
	gStyle->SetPadTickY(1);

//	Ca = new TCanvas( "Ca", "Ca", 1200, 680 );
	Ca = new TCanvas( "Ca", "Ca", 1202, 744 );
	Ca -> SetMargin( .100, .020, .040, .020 );
	//Left, Right, Bottom, Top

	Leg = new TLegend( .665, .690, .990, .990 );
	Set -> Setting_Legend( Leg, 42, 12, 1, .0225 );
//	Leg -> SetFillColor(kYellow-10);
	Leg -> SetFillColor(14);

	h_frame = new TH1D( "h_frame", "h_frame", 1, WMin, WMax );
	h_frame -> GetXaxis()->SetNdivisions(1);
	h_frame -> GetYaxis()->SetNdivisions(515);
	h_frame -> GetYaxis()->SetTitle("Lifetime of ^{3}_{#Lambda}H (ps)");
	h_frame -> GetYaxis()->SetTitleSize(.055);
	h_frame -> GetYaxis()->SetTitleOffset(.75);
	h_frame -> GetYaxis()->SetTitleFont(42);
	h_frame -> GetYaxis()->SetRangeUser(HMin, HMax);
	h_frame -> GetXaxis()->SetLabelSize(0.);
	h_frame -> GetYaxis()->SetLabelFont(42);
	h_frame -> SetTitle("");
	h_frame -> SetStats(kFALSE);

	h_frame->Draw("");

	for(int i=0; i<3; i++){
		Ln[i] = new TLine( WMin, LT_Val[i], WMax, LT_Val[i]);
		Set -> Setting_Line( Ln[i], LCol[i], 2, LSty[i] );
		Ln[i] -> Draw();
	}
	Leg->AddEntry( Ln[0], "Free #Lambda (PDG Val.): 263.2 #pm 2.0 ps"      , "l" );
	Leg->AddEntry( Ln[1], "H. Kamada #it{et. al.}, PRC #bf{57} (1998) 1595", "l" );
	Leg->AddEntry( Ln[2], "A. Gal, H. Garcilazo, PLB #bf{791} (2019) 48-53", "l" );

	for(int i=0; i<NofData; i++){
		HtlLT[i] = new HyperTritonLifetime(i);
		HtlLT[i] -> SetDT();
		HtlLT[i] -> SetDataPointFromDT();
		HtlLT[i] -> LTDraw();
	}

	for(int i=0; i<NofDummy; i++){
		gr_dummy[i] = new TGraphErrors();
		gr_dummy[i] -> SetLineColor(DummyMCol[i]);
		gr_dummy[i] -> SetMarkerStyle(DummyMSty[i]);
		gr_dummy[i] -> SetMarkerColor(DummyMCol[i]);
		gr_dummy[i] -> SetMarkerSize(2.00);
	}
	gr_dummy[2] -> SetMarkerSize(2.50);
	gr_dummy[3] -> SetMarkerSize(2.30);
	gr_dummy[4] -> SetMarkerSize(2.30);
	Leg -> AddEntry( gr_dummy[0], "Emulsion"                   , "pl" );
	Leg -> AddEntry( gr_dummy[1], "Bubble chamber"             , "pl" );
	Leg -> AddEntry( gr_dummy[2], "STAR" , "pl" );
	Leg -> AddEntry( gr_dummy[3], "HypHI", "pl" );
	Leg -> AddEntry( gr_dummy[4], "ALICE", "pl" );

	Fr = gPad->GetFrame();
	Fr -> SetFillStyle(0);
	Fr -> SetLineWidth(1);
	Fr -> Draw();
	
	Leg -> Draw();

	Ca->Print( figure_op.c_str(), "pdf");
	Ca->Print( figure.c_str() );
	Ca->Print( figure_cl.c_str(), "pdf");
	
	delete Set;
	for(int i=0; i<NofData; i++)delete  HtlLT[i];

	gSystem->Exit(-1);

	theApp->Run();

	return 0;
}
