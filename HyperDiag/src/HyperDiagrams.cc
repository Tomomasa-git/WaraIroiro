/* * * * * * * * * * * * 
 *  HyperDiagrams.cc   *
 *                     *
 *  T. Fujiwara        *
 *  2021. 01. 09 (Sat) *
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
#include "TPolyLine.h"
#include "TCurlyLine.h"
#include "TArrow.h"
#include "TEllipse.h"
#include "TArc.h"
#include "TBox.h"
#include "TPave.h"
#include "TPaveStats.h"
#include "TPaveText.h"

#include "./Setting.h"
#include "./MaruiChohokei.h"

const string figurepath="../fig/";

void Kpi_diag();
void piK_diag();
void eeK_diag();
void gamK_diag();


int main( int argc, char** argv ){
	
	TApplication *theApp;

	theApp = new TApplication( "App", &argc, argv );

	gROOT->SetBatch(kTRUE);

	Kpi_diag();
	piK_diag();
	eeK_diag();

	gSystem->Exit(-1);
	theApp->Run();

	return 0;
}



void Kpi_diag(){
	int ItNum=0;
	const int NofParticles=4;

	double QuarkPos[NofParticles][2] = {
	                          {0.240, 0.250},	//Neutron
	                          {0.240, 0.700},	//K-
	                          {0.760, 0.250},	//Lambda
	                          {0.760, 0.700}	//pi-
	                        };
	int ParType[NofParticles] = {3, 2, 3, 2};
	int Pid[NofParticles]     = {1, 3, 2, 4};

	Setting* MySetting;
	TCanvas *Ca;
	TLatex* Lat;
	TPaveText* Pt;
	MaruiChohokei* MyBox[NofParticles];
	TLine* Ln[5];

	TString figure=Form( "Hyp_Kpidiag_%03d.pdf", ItNum );
	figure = figurepath+figure;

	MySetting = new Setting();
	MySetting->Setting_Gene();

	Ca = new TCanvas("Ca","Ca",1202,924);
	Ca -> cd();

	Pt = new TPaveText( .040, .870, .490, .990 );
	MySetting->Setting_Pave( Pt, 52, 1, 22, 0.072, 1, 1, 1, 13, 1001 );
	Pt->AddText("n(K^{-}, #it{#pi}^{-})#Lambda reaction");

	Lat = new TLatex();
	MySetting->Setting_Latex(Lat,12,22,1, 0.200);

	for(int i=0; i<NofParticles; i++){MyBox[i] = new MaruiChohokei( QuarkPos[i][0], QuarkPos[i][1], ParType[i], Pid[i] );}
	
	Ln[0] = new TLine( MyBox[1]->GetQuarkPosX(1)+0.040, MyBox[1]->GetQuarkPosY(1), 
	                   MyBox[2]->GetQuarkPosX(0)-0.040, MyBox[2]->GetQuarkPosY(0) 
	                 );
	Ln[1] = new TLine( MyBox[0]->GetQuarkPosX(0)+0.040, MyBox[0]->GetQuarkPosY(0), 
	                   MyBox[3]->GetQuarkPosX(1)-0.040, MyBox[3]->GetQuarkPosY(1) 
	                 );
	Ln[2] = new TLine( MyBox[0]->GetQuarkPosX(1)+0.040, MyBox[0]->GetQuarkPosY(1), 
	                   MyBox[2]->GetQuarkPosX(1)-0.040, MyBox[2]->GetQuarkPosY(1) 
	                 );
	Ln[3] = new TLine( MyBox[0]->GetQuarkPosX(2)+0.040, MyBox[0]->GetQuarkPosY(2), 
	                   MyBox[2]->GetQuarkPosX(2)-0.040, MyBox[2]->GetQuarkPosY(2) 
	                 );
	Ln[4] = new TLine( MyBox[1]->GetQuarkPosX(0)+0.040, MyBox[1]->GetQuarkPosY(0), 
	                   MyBox[3]->GetQuarkPosX(0)-0.040, MyBox[3]->GetQuarkPosY(0) 
	                 );

	MySetting->Setting_Line( Ln[0], 1, 3, 1 );
	MySetting->Setting_Line( Ln[1], 1, 3, 1 );
	MySetting->Setting_Line( Ln[2], 1, 1, 1 );
	MySetting->Setting_Line( Ln[3], 1, 1, 1 );
	MySetting->Setting_Line( Ln[4], 1, 1, 1 );

	Ca -> cd();
	Pt->Draw();	

	for(int i=0; i<NofParticles; i++){
		MyBox[i]->Draw();
		MyBox[i]->DrawQuarks();
	}

	for(int i=0; i<5; i++){Ln[i]->Draw();}

	Lat->SetTextColor(MyBox[0]->GetParCol());	Lat->DrawLatex( QuarkPos[0][0]-0.160, QuarkPos[0][1], MyBox[0]->GetParName(999) );
	Lat->SetTextColor(MyBox[1]->GetParCol());	Lat->DrawLatex( QuarkPos[1][0]-0.160, QuarkPos[1][1], MyBox[1]->GetParName(999) );
	Lat->SetTextColor(MyBox[2]->GetParCol());	Lat->DrawLatex( QuarkPos[2][0]+0.160, QuarkPos[2][1], MyBox[2]->GetParName(999) );
	Lat->SetTextColor(MyBox[3]->GetParCol());	Lat->DrawLatex( QuarkPos[3][0]+0.160, QuarkPos[3][1], MyBox[3]->GetParName(999) );

	gPad->Update();
	Ca->Print(figure, "pdf");
	delete MySetting;
	for(int i=0; i<NofParticles; i++){delete MyBox[i];}
	Ca->Destructor();

	return;
}

void piK_diag(){
	int ItNum=0;
	const int NofParticles=4;

	double QuarkPos[NofParticles][2] = {
	                          {0.240, 0.250},	//Neutron
	                          {0.240, 0.700},	//K-
	                          {0.760, 0.250},	//Lambda
	                          {0.760, 0.700}	//pi-
	                        };
	int ParType[NofParticles] = {3, 2, 3, 2};
	int Pid[NofParticles]     = {1, 5, 2, 6};

	double ArcCent[2][2];
	double ArcRad;
	double ArcPhi[2]={90. , 270.};
	double ArcTheta[2]={180. , 0.};
	double CulPos[2][2];
	double CulAmp=0.0200;
	double CulWL=0.035;

	Setting* MySetting;
	TCanvas *Ca;
	TLatex* Lat;
	TPaveText* Pt;
	MaruiChohokei* MyBox[NofParticles];
	TLine* Ln[7];
	TArc* Arc[2];
	TCurlyLine* Cul;

	MySetting = new Setting();
	MySetting->Setting_Gene();

	TString figure=Form( "Hyp_piKdiag_%03d.pdf", ItNum );
	figure = figurepath+figure;

	Ca = new TCanvas("Ca","Ca",1202,924);
	Ca -> cd();

	Pt = new TPaveText( .040, .870, .490, .990 );
	MySetting->Setting_Pave( Pt, 52, 1, 22, 0.072, 1, 1, 1, 12, 1001 );
	Pt->AddText("n(#it{#pi^{+}}, K^{+})#Lambda reaction");

	Lat = new TLatex();
	MySetting->Setting_Latex(Lat,12,22,1, 0.200);

	for(int i=0; i<NofParticles; i++){MyBox[i] = new MaruiChohokei( QuarkPos[i][0], QuarkPos[i][1], ParType[i], Pid[i] );}
	
	Ln[0] = new TLine( 
	                   MyBox[1]->GetQuarkPosX(1)+0.040, MyBox[1]->GetQuarkPosY(1), 
	                   MyBox[1]->GetQuarkPosX(0)+0.080, MyBox[1]->GetQuarkPosY(1) 
	                 );
	Ln[1] = new TLine( 
	                   MyBox[0]->GetQuarkPosX(0)+0.040, MyBox[0]->GetQuarkPosY(0), 
	                   MyBox[0]->GetQuarkPosX(1)+0.080, MyBox[0]->GetQuarkPosY(0) 
	                 );
	Ln[2] = new TLine( 
	                   MyBox[3]->GetQuarkPosX(1)-0.040, MyBox[3]->GetQuarkPosY(1), 
	                   MyBox[3]->GetQuarkPosX(1)-0.080, MyBox[3]->GetQuarkPosY(1) 
	                 );
	Ln[3] = new TLine( 
	                   MyBox[2]->GetQuarkPosX(0)-0.040, MyBox[2]->GetQuarkPosY(0), 
	                   MyBox[2]->GetQuarkPosX(0)-0.080, MyBox[2]->GetQuarkPosY(0) 
	                 );
	Ln[4] = new TLine( 
	                   MyBox[0]->GetQuarkPosX(1)+0.040, MyBox[0]->GetQuarkPosY(1), 
	                   MyBox[2]->GetQuarkPosX(1)-0.040, MyBox[2]->GetQuarkPosY(1) 
	                 );
	Ln[5] = new TLine( 
	                   MyBox[0]->GetQuarkPosX(2)+0.040, MyBox[0]->GetQuarkPosY(2), 
	                   MyBox[2]->GetQuarkPosX(2)-0.040, MyBox[2]->GetQuarkPosY(2) 
	                 );
	Ln[6] = new TLine( 
	                   MyBox[1]->GetQuarkPosX(0)+0.040, MyBox[1]->GetQuarkPosY(0), 
	                   MyBox[3]->GetQuarkPosX(0)-0.040, MyBox[3]->GetQuarkPosY(0) 
	                 );

	for(int i=0; i<7; i++){
		int Wid;
		if(i<4){
			Wid=3;
		}else{
			Wid=1;
		}
		MySetting->Setting_Line( Ln[i], 1, Wid, 1 );
	}

	for(int i=0; i<2; i++){
		ArcCent[i][0]=Ln[2*i]->GetX2();
		ArcCent[i][1]=0.5*( (Ln[2*i]->GetY1())+(Ln[2*i+1]->GetY1()) );
	}
	ArcRad=0.5*( (Ln[0]->GetY1())-(Ln[1]->GetY1()) );

	for(int i=0; i<2; i++){
		Arc[i] = new TArc( ArcCent[i][0], ArcCent[i][1], ArcRad, ArcPhi[0], ArcPhi[1] );
		MySetting->Setting_Arc(Arc[i], 1, 3, 1);
		Arc[i]->SetNoEdges(kTRUE);
		Arc[i]->SetTheta(ArcTheta[i]);
		Arc[i]->SetR2(Arc[i]->GetR1());
		Arc[i]->SetR1(0.80*Arc[i]->GetR1());
		CulPos[i][0] = (Arc[i]->GetX1())+pow(-1.,(double)i)*(Arc[i]->GetR1());
		CulPos[i][1] =  Arc[i]->GetY1();
	}

	Cul = new TCurlyLine( CulPos[0][0], CulPos[0][1], CulPos[1][0], CulPos[1][1] );
	MySetting->Setting_Line( Cul, 1, 3, 1, CulAmp, CulWL );
//	Cul->SetAmplitude(CulAmp);
//	Cul->SetWaveLength(CulWL);
//	cout<<"**********"<<endl;
//	cout<<Cul->GetAmplitude()<<endl;
//	cout<<Cul->GetWaveLength()<<endl;
//	cout<<"**********"<<endl;
	
	Ca -> cd();
	Pt->Draw();	

	for(int i=0; i<NofParticles; i++){
		MyBox[i]->Draw();
		MyBox[i]->DrawQuarks();
	}

	for(int i=0; i<7; i++){Ln[i]->Draw();}
	for(int i=0; i<2; i++){Arc[i]->Draw();}
//	Cul->Draw();

	Lat->SetTextColor(MyBox[0]->GetParCol());	Lat->DrawLatex( QuarkPos[0][0]-0.160, QuarkPos[0][1], MyBox[0]->GetParName(999) );
	Lat->SetTextColor(MyBox[1]->GetParCol());	Lat->DrawLatex( QuarkPos[1][0]-0.160, QuarkPos[1][1], MyBox[1]->GetParName(999) );
	Lat->SetTextColor(MyBox[2]->GetParCol());	Lat->DrawLatex( QuarkPos[2][0]+0.160, QuarkPos[2][1], MyBox[2]->GetParName(999) );
	Lat->SetTextColor(MyBox[3]->GetParCol());	Lat->DrawLatex( QuarkPos[3][0]+0.160, QuarkPos[3][1], MyBox[3]->GetParName(999) );

	gPad->Update();
	Ca->Print(figure, "pdf");

	delete MySetting;
	for(int i=0; i<NofParticles; i++){delete MyBox[i];}
	Ca->Destructor();

	return;
}

void eeK_diag(){
	int ItNum=1;
	const int NofParticles=3;

	double QuarkPos[NofParticles][2] = {
	                          {0.240, 0.250},	//Neutron
	                          {0.760, 0.250},	//Lambda
	                          {0.760, 0.700}	//pi-
	                        };
	int ParType[NofParticles] = {3, 3, 2};
	int Pid[NofParticles]     = {0, 2, 6};

	double TSize;

	double ElPos[2][3];
	double GamPos[2];
	double ArcCent[2];
	double ArcRad;
	double ArcPhi[2]={90. , 270.};
	double ArcTheta=0.;
	double CulPos[2][2];
	double CulAmp=0.0200;
	double CulWL=0.025;

	Setting* MySetting;
	TCanvas* Ca;
	TH1D* h_Draw;
//	TPad* Pad_Draw;
	TLatex* Lat;
	TPaveText* Pt;
	MaruiChohokei* MyBox[NofParticles];
	TLine* Ln[5];
	TPolyLine* Pl;
	TArc* Arc;
	TCurlyLine* Cul;

	MySetting = new Setting();
	MySetting->Setting_Gene();

	TString figure=Form( "Hyp_eeKdiag_%03d.pdf", ItNum );
	figure = figurepath+figure;

	Ca        = new TCanvas( "Ca",       "Ca"       , 1202, 1014 );
	Ca->SetMargin(0. ,0., 0., 0.);
	Ca->SetFrameLineColor(0);
	Ca->SetFrameLineWidth(0);
	Ca->SetFrameFillColor(0);
	Ca->SetFrameFillStyle(0);
	h_Draw = new TH1D( "h_Draw", "h_Draw", 1, 0., 1. );
	MySetting->Setting_Hist1D( h_Draw, "", "", "", 0, 0, 42, 0, 0 );
	h_Draw->SetStats(kFALSE);
	h_Draw->GetXaxis()->SetNdivisions(0);
	h_Draw->GetYaxis()->SetNdivisions(0);
	h_Draw->GetXaxis()->SetAxisColor(0);
	h_Draw->GetYaxis()->SetAxisColor(0);
//	h_Draw->GetYaxis()->SetRangeUser(0., 4./3.);
	h_Draw->GetYaxis()->SetRangeUser(0., 1.1 );
	h_Draw->Draw();

//	Pad_Draw  = new TPad( "Pad_Draw","Pad_Draw", 0., 0., 1., 0.75  );
//	Pad_Draw -> SetFillColor(0);
//	Pad_Draw -> SetFillStyle(0);
//	Pad_Draw -> cd();

	Pt = new TPaveText( .040, 0.971, .510, 1.080 );
	MySetting->Setting_Pave( Pt, 52, 1, 22, 0.0654, 1, 1, 1, 11, 1001 );
	Pt->AddText("p(e, e'K^{+})#Lambda reaction");

	Lat = new TLatex();
	MySetting->Setting_Latex(Lat,12,22,1, 0.165);

	for(int i=0; i<NofParticles; i++){
		MyBox[i] = new MaruiChohokei( QuarkPos[i][0], QuarkPos[i][1], ParType[i], Pid[i] );
		TSize = MyBox[i]->GetTextSize();
		MyBox[i]->SetTextSize( (10./11.)*TSize);
	}

	ElPos[0][0]	= MyBox[0]->GetQuarkPosX(0)-0.120;
	ElPos[0][1] = 0.350;
	ElPos[0][2] = MyBox[2]->GetQuarkPosX(0)+0.075;
	ElPos[1][0]	= MyBox[2]->GetQuarkPosY(0)+0.100;
	ElPos[1][1] = ElPos[1][0];
	ElPos[1][2] = MyBox[2]->GetQuarkPosY(0)+0.200;

	Pl = new TPolyLine(3, ElPos[0], ElPos[1]);
	MySetting->Setting_Line( Pl, 1, 1, 1 );
	
	Ln[0] = new TLine( 
	                   MyBox[0]->GetQuarkPosX(0)+0.040, MyBox[0]->GetQuarkPosY(0), 
	                   MyBox[2]->GetQuarkPosX(0)-0.040, MyBox[2]->GetQuarkPosY(0) 
	                 );
	Ln[1] = new TLine( 
	                   MyBox[2]->GetQuarkPosX(0)-0.040, MyBox[2]->GetQuarkPosY(1), 
	                   MyBox[2]->GetQuarkPosX(1)-0.080, MyBox[2]->GetQuarkPosY(1) 
	                 );
	Ln[2] = new TLine( 
	                   MyBox[1]->GetQuarkPosX(0)-0.040, MyBox[1]->GetQuarkPosY(0), 
	                   MyBox[1]->GetQuarkPosX(0)-0.080, MyBox[1]->GetQuarkPosY(0) 
	                 );
	Ln[3] = new TLine( 
	                   MyBox[0]->GetQuarkPosX(1)+0.040, MyBox[0]->GetQuarkPosY(1), 
	                   MyBox[1]->GetQuarkPosX(1)-0.040, MyBox[1]->GetQuarkPosY(1) 
	                 );
	Ln[4] = new TLine( 
	                   MyBox[0]->GetQuarkPosX(2)+0.040, MyBox[0]->GetQuarkPosY(2), 
	                   MyBox[1]->GetQuarkPosX(2)-0.040, MyBox[1]->GetQuarkPosY(2) 
	                 );

	for(int i=0; i<5; i++){
		int Wid;
		if(i<3){
			Wid=3;
		}else{
			Wid=1;
		}
		MySetting->Setting_Line( Ln[i], 1, Wid, 1 );
	}
	
	ArcCent[0]=Ln[2]->GetX2();
	ArcCent[1]=0.5*( (Ln[1]->GetY1())+(Ln[2]->GetY1()) );
	
	ArcRad=0.5*( (Ln[1]->GetY1())-(Ln[2]->GetY1()) );
	if(ArcRad<0.){
		ArcRad = abs(ArcRad);
	}else;

	Arc = new TArc( ArcCent[0], ArcCent[1], ArcRad, ArcPhi[0], ArcPhi[1] );
	MySetting->Setting_Arc(Arc, 1, 3, 1);
	Arc->SetNoEdges(kTRUE);
	Arc->SetTheta(ArcTheta);
	Arc->SetR2(Arc->GetR1());
	Arc->SetR1(0.80*Arc->GetR1());
	CulPos[1][0] = (Arc->GetX1())-(Arc->GetR1());
	CulPos[1][1] =  Arc->GetY1();
	CulPos[0][0] = ElPos[0][1];
	CulPos[0][1] = ElPos[1][1];

	//Gamma*
	GamPos[0] = CulPos[0][0]-0.025;//0.5*(CulPos[0][0]+CulPos[1][0]);
	GamPos[1] = 0.5*(CulPos[0][1]+CulPos[1][1]);

	Cul = new TCurlyLine( CulPos[0][0], CulPos[0][1], CulPos[1][0], CulPos[1][1] );
	Cul -> SetWavy();
	MySetting->Setting_Line( Cul, 602, 3, 1, CulAmp, CulWL );

	for(int i=0; i<NofParticles; i++){
		MyBox[i]->Draw();
		MyBox[i]->DrawQuarks();
	}
	
	Lat->SetTextColor(MyBox[0]->GetParCol());	Lat->DrawLatex( QuarkPos[0][0]-0.160, QuarkPos[0][1]   , MyBox[0]->GetParName(999) );
	Lat->SetTextColor(MyBox[1]->GetParCol());	Lat->DrawLatex( QuarkPos[1][0]+0.160, QuarkPos[1][1]   , MyBox[1]->GetParName(999) );
	Lat->SetTextColor(MyBox[2]->GetParCol());	Lat->DrawLatex( QuarkPos[2][0]+0.160, QuarkPos[2][1]   , MyBox[2]->GetParName(999) );
	Lat->SetTextColor(870);                    	Lat->DrawLatex( QuarkPos[0][0]-0.160, ElPos[1][0]      , "e"                       );
	Lat->SetTextColor(870);                    	Lat->DrawLatex( QuarkPos[1][0]+0.160, ElPos[1][2]+0.025, "e'"                      );
	Lat->SetTextSize(0.135);
	Lat->SetTextColor(602);                    	Lat->DrawLatex( GamPos[0]           , GamPos[1]        , "#it{#gamma}*"         );
	
	for(int i=0; i<5; i++){Ln[i]->Draw();}
	Arc->Draw();
	Cul->Draw();

	gPad->Update();

	Ca->cd();
//	Pad_Draw->Draw();
	Pl->Draw();
	
	Pt->Draw();

	Ca->Print(figure, "pdf");
//	cout<<"**********"<<endl;
//	cout<<Ca->GetFrameLineColor()<<endl;
//	cout<<"**********"<<endl;

	delete MySetting;
	for(int i=0; i<NofParticles; i++){delete MyBox[i];}
	Ca->Destructor();

	return;
}
