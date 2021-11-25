/* * * * * * * * * * * *
 * TransDecayExp.cc    *
 *	2021. 06. 05 (Sat) *
 *	T. FUJIWARA        *
 * * * * * * * * * * * */
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <iomanip>
#include <cmath>
#include <vector>
using namespace std;

#include "TApplication.h"
#include "TApplicationImp.h"
#include "TCanvasImp.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TRint.h"
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
#include "TGraph2D.h"
#include "TGraph2DErrors.h"
#include "TAxis.h"
#include "TGaxis.h"
#include "TRandom.h"
#include "TMath.h"
#include "TSpectrum.h"
#include "TString.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TPolyLine.h"
#include "TBox.h"
#include "TPave.h"
#include "TPaveStats.h"
#include "TPaveText.h"
#include "TStopwatch.h"

#include "./Setting.hh"

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

double TransExp(double *x, double *par);
double DecayExp(double *x, double *par);
double TransDecayExp(double *x, double *par);

double TransExp(double x, double *par);
double DecayExp(double x, double *par);
double TransDecayExp(double x, double *par);

int SetHistGrLineCol(int i=-999);
void RedrawFrame(TFrame *frame);
TString DrawOption(int Num=0, int Type=0);

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

int main( int argc, char** argv ){
	int ItNum=3;
	int const NofCanv = 3;
	int const NofGr   = 6;
	int NofGr_Draw[NofCanv] = {1, NofGr, NofGr};

	TApplication *theApp;
	Setting *MySetting;
	TCanvas *Ca[NofCanv];
	TF1 *f[3][NofGr];
	TGraph *g[3][NofGr];
	TLegend *Leg[3];
	TFrame *fr[NofCanv];
	TFile *ifp;
	TStopwatch *Stw;
	
//	double par00=1.;
	double par00[NofGr]={4.00,2.60,1.75,1.00,1.00,1.00};
	double par01=-0.10;
	double par02=0.;
	double par03=0.;
	double par10=1.;
	double par11[NofGr]={-0.10, -0.20, -0.50, -1.00, -2.00, -5.00};
	double par12=0.;
	double par13=0.;
	double OptimumPos[NofGr];
	
	double param[NofGr][8];
	
	TString LegEntry;
	TString LegEntry_;

	TString fig;
	TString root;
	
	int lcol;
	int lsty[NofGr]={1,2,5,8,6,3};
	int Npx=9E+4;

	double x1[NofCanv]={.64, .64, .64};
	double x2[NofCanv]={.94, .94, .94};
	double y1[NofCanv]={.20, .70, .20};
	double y2[NofCanv]={.50, .90, .50};

	fig  = Form("./fig/TransDecayExp_%03d.pdf" ,ItNum);
	root = Form("./fig/TransDecayExp_%03d.root",ItNum);

	theApp = new TApplication("App", &argc, argv );

	MySetting = new Setting();
	MySetting -> Setting_Gene(1);
	gStyle->SetStatW(0.25);
	gStyle->SetStatBorderSize(3);
	gStyle->SetPadBottomMargin(.130);

	Stw = new TStopwatch();
	Stw->Start();

	for(int i=0; i<NofGr; i++){
		par00[i] = 1./( (par11[i]/(par01+par11[i]))*pow( (par01/(par11[i]+par01)),par01/par11[i]) );

		param[i][0]=par00[i];
		param[i][1]=par01;
		param[i][2]=par02;
		param[i][3]=par03;
		param[i][4]=par10;
		param[i][5]=par11[i];
		param[i][6]=par12;
		param[i][7]=par13;

		OptimumPos[i] = (-1.*param[i][1])*log( param[i][1]/(param[i][1]+param[i][5]) );
		cout<<OptimumPos[i]<<" "<<par00[i]<<endl;
	}

	for(int i=0; i<NofCanv; i++)Ca[i] = new TCanvas(Form("Ca[%d]",i),Form("Ca[%d]",i),1282,744);

	for(int i=0; i<NofCanv; i++){
		Leg[i] = new TLegend(.64, .94, .20, .50);
		MySetting->Setting_Legend(Leg[i],62,12,602,0.035);
		Leg[i]->SetLineWidth(0);
		Leg[i]->SetFillStyle(0);
		Leg[i]->SetBorderSize(0);
		Leg[i]->SetX1(x1[i]);
		Leg[i]->SetX2(x2[i]);
		Leg[i]->SetY1NDC(y1[i]);
		Leg[i]->SetY2NDC(y2[i]);
	}

	cout<<NofCanv<<endl;
	cout<<NofGr  <<endl;

	for(int i=0; i<NofGr; i++){
		lcol=SetHistGrLineCol(i);

		f[0][i] = new TF1(Form("f[0][%d]",i),TransExp,-2.,10.,4);
		MySetting->Setting_Func(f[0][i],lcol,lsty[i]);
		f[0][i]->SetNpx(Npx);
		f[0][i]->SetParameters(par00[i],par01,par02,par03);

		g[0][i] = new TGraph(f[0][i]);
		MySetting->Setting_Graph(g[0][i],"","x","f(x)",lcol,lsty[i]);
		g[0][i]->GetXaxis()->SetNdivisions(9,2,5);
		g[0][i]->GetXaxis()->SetLimits(-1.,8.);
		g[0][i]->GetYaxis()->SetRangeUser(0.,1.1);
		g[0][i]->SetLineWidth(2);

		LegEntry = Form("#color[%d]{#it{p_{0}}=%.3lf, #it{p_{1}}=%.3lf}", lcol,par00[i],par01);
		if(i==0)Leg[0]->AddEntry(g[0][i], LegEntry, "l");

		f[1][i] = new TF1(Form("f[1][%d]",i),DecayExp,-2.,10.,4);
		MySetting->Setting_Func(f[1][i],lcol,lsty[i]);
		f[1][i]->SetNpx(Npx);
		f[1][i]->SetParameters(par10,par11[i],par12,par13);

		g[1][i] = new TGraph(f[1][i]);
		MySetting->Setting_Graph(g[1][i],"","x","f(x)",lcol,lsty[i]);
		g[1][i]->GetXaxis()->SetNdivisions(9,2,5);
		g[1][i]->GetXaxis()->SetLimits(-1.,8.);
		g[1][i]->GetYaxis()->SetRangeUser(0.,1.1);
		g[1][i]->SetLineWidth(2);

		LegEntry = Form("#color[%d]{#it{p_{1}}=%.3lf}", lcol,par11[i]);
		Leg[1]->AddEntry(g[1][i], LegEntry, "l");

		f[2][i] = new TF1(Form("f[2][%d]",i),TransDecayExp,-2.,10.,8);
		MySetting->Setting_Func(f[2][i],lcol,lsty[i]);
		f[2][i]->SetNpx(Npx);
		f[2][i]->SetParameters(param[i]);

		g[2][i] = new TGraph(f[2][i]);
		MySetting->Setting_Graph(g[2][i],"","x","f(x)",lcol,lsty[i]);
		g[2][i]->GetXaxis()->SetNdivisions(9,2,5);
		g[2][i]->GetXaxis()->SetLimits(-1.,8.);
		g[2][i]->GetYaxis()->SetRangeUser(0.,1.1);
		g[2][i]->SetLineWidth(2);

		LegEntry = Form("#color[%d]{#it{p_{1}}=%.3lf, #it{p_{5}}=%.3lf}", lcol,par01, par11[i]);
		Leg[2]->AddEntry(g[1][i], LegEntry, "l");
	}

	for(int i=0; i<NofCanv; i++){
		Ca[i]->cd();
		for(int j=0; j<NofGr_Draw[i]; j++)g[i][j]->Draw(DrawOption(j));
		RedrawFrame(fr[i]);
		Leg[i]->Draw();
	}

	ifp = new TFile(root, "recreate");
	ifp->cd();
	for(int i=0; i<NofCanv; i++)Ca[i]   -> Write( Form("Ca_%02d",i) );
	for(int i=0; i<3; i++){ for(int j=0; j<NofGr; j++)f[i][j] -> Write( Form("f_%d%d"  ,i,j) ); }
	for(int i=0; i<3; i++){ for(int j=0; j<NofGr; j++)g[i][j] -> Write( Form("g_%d%d"  ,i,j) ); }
	ifp->Close();

	Ca[0]->Print(fig+"[","pdf");
	for(int i=0; i<NofCanv; i++)Ca[i]->Print(fig);
	Ca[NofCanv-1]->Print(fig+"]","pdf");

	Stw->Stop();
	cout<<"Toral time: "<<Stw->RealTime()<<" sec"<<endl;

	delete MySetting;
	gSystem->Exit(-1);

	theApp->Run();

	return 0;
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

double TransExp(double *x, double *par){
	double f;
	if(x[0]>par[2])f=par[0]*( 1.-exp( (x[0]-par[2])/par[1] ))+par[3];
	else f=par[3];

	return f;
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

double DecayExp(double *x, double *par){
	double f;
	if(x[0]>par[2])f=par[0]*exp( (x[0]-par[2])/par[1] )+par[3];
	else f=par[3];

	return f;
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

double TransExp(double x, double *par){
	double f;
	if(x>par[2])f=par[0]*( 1.-exp( (x-par[2])/par[1] ))+par[3];
	else f=par[3];

	return f;
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

double DecayExp(double x, double *par){
	double f;
	if(x>par[2])f=par[0]*exp( (x-par[2])/par[1] )+par[3];
	else f=par[3];

	return f;
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

double TransDecayExp(double *x, double *par){
	double f;
	f = TransExp(x[0],par)*DecayExp(x[0],&par[4]);
	return f;
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

int SetHistGrLineCol(int i){
	int LCol=1;

	if(i<0)LCol=602;
	else if( (i<=9)&&(i>=0) ){
		switch(i){
			case 0:
			LCol=602;
			break;

			case 1:
			LCol=2;
			break;

			case 2:
			LCol=3;
			break;

			case 3:
			LCol=4;
			break;

			case 4:
			LCol=6;
			break;

			case 5:
			LCol=870;
			break;
			
			case 6:
			LCol=800;
			break;

			case 7:
			LCol=8;
			break;
		
			case 8:
			LCol=9;
			break;
			
			case 9:
			LCol=880;
			break;
		
			default:
			LCol=i;
			break;
		}
	}else{
		int i_=i-9;
		switch(i_){
			case 0:
			LCol=602;
			break;

			case 1:
			LCol=2;
			break;

			case 2:
			LCol=3;
			break;

			case 3:
			LCol=4;
			break;

			case 4:
			LCol=6;
			break;

			case 5:
			LCol=870;
			break;
			
			case 6:
			LCol=800;
			break;

			case 7:
			LCol=8;
			break;
		
			case 8:
			LCol=9;
			break;
		
			case 9:
			LCol=880;
			break;
		
			default:
			LCol=i_;
			break;
		}
	}

	return LCol;
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void RedrawFrame(TFrame *frame){
	gPad->Update();
	gPad->Modified();

	frame=gPad->GetFrame();
	frame->SetFillStyle(0);
	frame->Draw();

	gPad->Update();
	gPad->Modified();
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

TString DrawOption(int Num, int Type){
	//Type:
	//0 >> Only Line
	//1 >> Only Marker
	//2 >> Line and Marker
	int DrawOrder=0;
	TString Main="L";
	TString opt;

	if(Num<0)DrawOrder=0;
	else     DrawOrder=Num;

	switch(Type){
		case 0:
		Main="L";
		break;

		case 1:
		Main="P";
		break;

		case 2:
		Main="PL";
		break;

		default:
		Main="L";
		break;
	}

	if(DrawOrder==0)opt="A"+Main;
	else            opt="same"+Main;
	
	return opt;	
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
