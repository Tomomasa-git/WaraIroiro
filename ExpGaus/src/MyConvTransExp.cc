/* * * * * * * * * * * *
 * MyConvTransExp.cc   *
 *	2021. 06. 06 (Sun) *
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
double TransExpGaus(double *x, double *par);
double TransDecayExpGaus(double *x, double *par);

double TransExp(double x, double *par);
double DecayExp(double x, double *par);
double TransDecayExp(double x, double *par);

int SetHistGrLineCol(int i=-999);
int SetHistGrFillCol(int i=-999);
void RedrawFrame(TFrame *frame);
TString DrawOption(int Num=0, int Type=0);

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

int main( int argc, char** argv ){
	int ItNum=2;
	int const NofCanv = 5;
	int const NofHist = 5;
	int const MaxItr=1E+7;

	TApplication *theApp;
	Setting *MySetting;
	TCanvas *Ca[NofCanv];
	TH1D *hRaw[NofHist], *hConv[NofHist];
//	TF1 *f[NofHist];
	TFrame *fr[NofHist][2];
	TLegend *Leg[NofHist][2];
	TFile *ifp;
	TStopwatch *Stw;

	//Exp
//	double par0= 1.;
	double par0= 1.;
	double par1=-1.;
	double par2= 0.;
	double par3= 0.;
	double param[4]={
	                 par0,
	                 par1,
	                 par2,
	                 par3
	                };

	//Gaus
	double sigma[NofHist]={0.1, 0.25, 0.50, 1.00, 2.00};
	
	double x;
	double y;
	double Border;
	double delta;

	//Hist Setting
	double HistoRange[2]={-2., 20.};
	int HistoNbin=440;
	double HistoDiv;
	double HistoLCol;
	double HistoFCol;
	double HistoFSty=1001;
	float HistoFTransp=0.50;

	TString fig;
	TString root;
	
	theApp = new TApplication("App", &argc, argv );
	MySetting = new Setting();
	MySetting -> Setting_Gene(1);
	Stw = new TStopwatch();

	gRandom->SetSeed(0);

	fig=Form("./fig/MyConvTransExp_%03d.pdf",ItNum);
	root=Form("./fig/MyConvTransExp_%03d.root",ItNum);

	//////////////////
	//    Define    //
	//////////////////
	//Canvavs
	for(int i=0; i<NofCanv; i++){
		Ca[i] = new TCanvas( Form("Ca[%d]",i), Form("Ca[%d]",i), 1282,1282 );
		Ca[i]->Divide(1,2);
	}

	//Histo
	HistoDiv = (HistoRange[1]-HistoRange[0])/(double)HistoNbin;
	for(int i=0; i<NofHist; i++){
		HistoLCol=SetHistGrLineCol(i);
		HistoFCol=SetHistGrFillCol(i);

		for(int j=0; j<2; j++){
			Leg[i][j] = new TLegend(.20,.80,.40,.90);
			MySetting->Setting_Legend(Leg[i][j],62,12,HistoLCol,.045);
			Leg[i][j]->SetLineWidth(0);
			Leg[i][j]->SetFillStyle(0);
			Leg[i][j]->SetBorderSize(0);
		}

		hRaw[i]  = new TH1D(Form("hRaw[%d]" ,i) ,Form("hRaw[%d]" ,i) , HistoNbin,-2.,20.);
		hConv[i] = new TH1D(Form("hConv[%d]",i) ,Form("hConv[%d]",i) , HistoNbin,-2.,20.);

		MySetting->Setting_Hist1D( hRaw[i],"","x",Form("Count/bin #scale[0.60]{(bin=%.3lf)}",HistoDiv),HistoLCol,1,62,HistoFCol,HistoFSty );
		hRaw[i]->SetFillColorAlpha(HistoFCol,HistoFTransp);
		hRaw[i]->GetXaxis()->SetLimits(-2.,10.);
		hRaw[i]->GetXaxis()->SetNdivisions(12,2,5);
		hRaw[i]->GetYaxis()->SetNoExponent();
		Leg[i][0]->AddEntry(hRaw[i],Form("#tau=%.3lf", par1),"fl");

		MySetting->Setting_Hist1D( hConv[i],"","x",Form("Count/bin #scale[0.60]{(bin=%.3lf)}",HistoDiv),HistoLCol,1,62,HistoFCol,HistoFSty );
		hConv[i]->SetFillColorAlpha(HistoFCol,HistoFTransp);
		hConv[i]->GetXaxis()->SetLimits(-2.,10.);
		hConv[i]->GetXaxis()->SetNdivisions(12,2,5);
		hConv[i]->GetYaxis()->SetNoExponent();
		Leg[i][1]->AddEntry(hConv[i],Form("#tau=%.3lf, #sigma=%.3lf", par1,sigma[i]),"fl");
	}


	//////////////////
	//	Fill value	//
	//////////////////
	Stw->Start();
	for(int i=0; i<MaxItr; i++){
		x = gRandom->Uniform(HistoRange[0], HistoRange[1]);
		y = gRandom->Uniform(0.,1.);
		Border=TransExp(x,param);
		
		if(y<=Border){
			for(int j=0; j<NofHist; j++){
				hRaw[j]->Fill(x);
			}
			for(int j=0; j<NofHist; j++){
				gRandom->SetSeed(0);
				delta=gRandom->Gaus(0.,sigma[j]);
				hConv[j]->Fill(x+delta);
				if( (i%100000)==0 )cout<<i<<" "<<j<<endl;
				delta=0.;
			}
		}else;
		
	}
	Stw->Stop();
	cout<<"=========="<<endl;
	cout<<"Total Time: "<<Stw->RealTime()<<" sec"<<endl;
	cout<<"=========="<<endl;


	//////////////////
	//	Draw Histo	//
	//////////////////
	for(int i=0; i<NofCanv; i++){
		Ca[i]->cd(1);
		hRaw[i]->Draw();
		RedrawFrame(fr[i][0]);
		Leg[i][0]->Draw();

		Ca[i]->cd(2);
		hConv[i]->Draw();
		RedrawFrame(fr[i][1]);
		Leg[i][1]->Draw();
	}
	
	////////////
	//  Save  //
	////////////
	//As PDF
	Ca[0]->Print(fig+"[", "pdf");
	for(int i=0; i<NofCanv; i++)Ca[i]->Print(fig);
	Ca[NofCanv-1]->Print(fig+"]", "pdf");

	//As Root file
	ifp = new TFile(root,"recreate");
	ifp->cd();
	for(int i=0; i<NofCanv; i++)Ca[i]->Write( Form("Ca_%02d",i) );
	for(int i=0; i<NofHist; i++)hRaw[i]  -> Write( Form("hRaw_%02d" ,i) );
	for(int i=0; i<NofHist; i++)hConv[i] -> Write( Form("hConv_%02d",i) );
	ifp->Close();

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

	if(i<0)LCol=870;
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

int SetHistGrFillCol(int i){
	int FCol=1;

	if(i<0)FCol=602;
	else if( (i<=9)&&(i>=0) ){
		switch(i){
			case 0:
			FCol=870;
			break;

			case 1:
			FCol=6;
			break;

			case 2:
			FCol=3;
			break;

			case 3:
			FCol=7;
			break;

			case 4:
			FCol=6;
			break;

			case 5:
			FCol=7;
			break;
			
			case 6:
			FCol=800;
			break;

			case 7:
			FCol=840;
			break;
		
			case 8:
			FCol=856;
			break;
			
			case 9:
			FCol=871;
			break;
		
			default:
			FCol=i;
			break;
		}
	}else{
		int i_=i-9;
		switch(i_){
			case 0:
			FCol=870;
			break;

			case 1:
			FCol=6;
			break;

			case 2:
			FCol=3;
			break;

			case 3:
			FCol=7;
			break;

			case 4:
			FCol=6;
			break;

			case 5:
			FCol=7;
			break;
			
			case 6:
			FCol=800;
			break;

			case 7:
			FCol=840;
			break;
		
			case 8:
			FCol=856;
			break;
			
			case 9:
			FCol=871;
			break;
		
			default:
			FCol=i_;
			break;
		}
	}

	return FCol;
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

