/* * * * * * * * * * * * 
 *  MaruiChohokei.cc   *
 *                     *
 *  T. Fujiwara        *
 *  2021. 01. 10 (Sun) *
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
#include "TArc.h"
#include "TBox.h"
#include "TPave.h"
#include "TPaveStats.h"
#include "TPaveText.h"

#include "./Setting.hh"
#include "./MaruiChohokei.h"

MaruiChohokei::MaruiChohokei( double gx, double gy, int N, int pid ){
	MySetting = new Setting();
	MySetting -> Setting_Gene();

	gPos[0] = gx;
	gPos[1] = gy;

	ParticleID=pid;
	NofQuarks=N;
	
	switch(N){
		case 2:
			Vsize=0.280;
			break;

		case 3:
			Vsize=0.420;
			break;

		default:
			NofQuarks=3;
			Vsize=0.30;
			break;
	}

	Rad   = 0.075;
//	Rad_V = 0.060;
//	for(int i=0; i<4; i++){
//		AngCenter[i][0] = gx+pow(-1.,(double)i+1.)*0.5*Hsize;
//	}
	AngCenter[0][0] = gx-0.5*Hsize;
	AngCenter[1][0] = gx+0.5*Hsize;
	AngCenter[2][0] = gx-0.5*Hsize;
	AngCenter[3][0] = gx+0.5*Hsize;
	AngCenter[0][1] = gy-0.5*Vsize+Rad;
	AngCenter[1][1] = gy-0.5*Vsize+Rad;
	AngCenter[2][1] = gy+0.5*Vsize-Rad;
	AngCenter[3][1] = gy+0.5*Vsize-Rad;

	Theta[0] = 180.;
	Theta[1] = 270.;
	Theta[2] =  90.;
	Theta[3] =   0.;

	//Line positions
	LnPos[0][0][0] = AngCenter[0][0];
	LnPos[0][0][1] = AngCenter[0][1]-Rad;
	LnPos[0][1][0] = AngCenter[1][0];
	LnPos[0][1][1] = AngCenter[1][1]-Rad;

	LnPos[1][0][0] = AngCenter[1][0]+Rad;
	LnPos[1][0][1] = AngCenter[1][1];
	LnPos[1][1][0] = AngCenter[3][0]+Rad;
	LnPos[1][1][1] = AngCenter[3][1];

	LnPos[2][0][0] = AngCenter[2][0];
	LnPos[2][0][1] = AngCenter[2][1]+Rad;
	LnPos[2][1][0] = AngCenter[3][0];
	LnPos[2][1][1] = AngCenter[3][1]+Rad;

	LnPos[3][0][0] = AngCenter[0][0]-Rad;
	LnPos[3][0][1] = AngCenter[0][1];
	LnPos[3][1][0] = AngCenter[2][0]-Rad;
	LnPos[3][1][1] = AngCenter[2][1];

	VerticalRange = LnPos[2][0][1]-LnPos[0][0][1];

	for(int i=0; i<4; i++){
		Ell[i] = new TEllipse( AngCenter[i][0], AngCenter[i][1], Rad, Rad, 0., 90., Theta[i] );
		Ln[i] = new TLine( LnPos[i][0][0], LnPos[i][0][1], LnPos[i][1][0], LnPos[i][1][1] );
	}

	Lat = new TLatex();
	MySetting->Setting_Latex(Lat, 132, 22, 1, 0.120 );

	SetQuarks();
	SetAttribute( LColor, LStyle, LWidth );
}

MaruiChohokei::~MaruiChohokei(){
	for(int i=0; i<NofAngle; i++){
		delete Ln[i];
		delete Ell[i];
	}
	delete Lat;
	delete MySetting;
}


void MaruiChohokei::SetAttribute( int LCol, int LSty, int LWid ){
	for(int i=0; i<NofAngle; i++){
		MySetting->Setting_Line( Ln[i], LCol, LWid, LSty );

		Ell[i]->SetLineColor(LCol);
		Ell[i]->SetLineStyle(LSty);
		Ell[i]->SetLineWidth(LWid);
		Ell[i]->SetNoEdges(kTRUE);
	}
}


void MaruiChohokei::Draw(){
	for(int i=0; i<NofAngle; i++){
		Ln[i] -> Draw();
		Ell[i]-> Draw();
	}
}

void MaruiChohokei::DrawQuarks(){
	
	for(int i=0; i<NofQuarks; i++){
		Lat->DrawLatex(CharPos[i][0], CharPos[i][1], Quark[i]);
	}
	
}

void MaruiChohokei::SetQuarks(){
	switch(ParticleID){
		//Proton
		case 0:
		Quark[0] = "#it{u}";
		Quark[1] = "#it{u}";
		Quark[2] = "#it{d}";
		ParName  = "p";
		LColor = 2;
		break;

		//Neutron
		case 1:
		Quark[0] = "#it{d}";
		Quark[1] = "#it{u}";
		Quark[2] = "#it{d}";
		ParName  = "n";
		LColor = 4;
		break;

		//Lambda
		case 2:
		Quark[0] = "#it{s}";
		Quark[1] = "#it{u}";
		Quark[2] = "#it{d}";
		ParName  = "#Lambda";
		LColor = 6;
		break;

		//K-
		case 3:
		Quark[0] = "#bar{#it{u}}";
		Quark[1] = "#it{s}";
		ParName  = "K^{-}";
	//	LColor = 800;
		LColor = 843;
		break;
		
		//pi-
		case 4:
		Quark[0] = "#bar{#it{u}}";
		Quark[1] = "#it{d}";
		ParName  = "#it{#pi}^{-}";
		LColor = 841;
		break;

		//pi+
		case 5:
		Quark[0] = "#it{u}";
		Quark[1] = "#bar{#it{d}}";
		ParName  = "#it{#pi}^{+}";
		LColor = 880;
		break;

		//K+
		case 6:
		Quark[0] = "#it{u}";
		Quark[1] = "#bar{#it{s}}";
		ParName  = "K^{+}";
		LColor = 3;
		break;

		//Default
		default:
		Quark[0] = "#it{s}";
		Quark[1] = "#it{u}";
		Quark[2] = "#it{d}";
		ParName  = "#Lambda";
		LColor = 6;
		break;
	}
	
	switch(NofQuarks){
		case 2:
		for(int i=0; i<NofQuarks; i++){
			CharPos[i][0] = gPos[0];
			CharPos[i][1] = gPos[1]+pow(-1., (double)i)*0.25*VerticalRange;
		}
		break;

		case 3:
		for(int i=0; i<NofQuarks; i++){
			CharPos[i][0] = gPos[0];
			CharPos[i][1] = gPos[1]-(1./3.)*(double)i*VerticalRange+(1./3.)*VerticalRange;
		}
		break;

		default:
		for(int i=0; i<NofQuarks; i++){
			CharPos[i][0] = gPos[0];
			CharPos[i][1] = gPos[1]-(1./3.)*(double)i*VerticalRange+(1./3.)*VerticalRange;
		}
		break;
	}
}

double MaruiChohokei::GetQuarkPosX(int i){
	double ret;
	if(i>NofQuarks){
		cout<<"*************************"<<endl;
		cout<<"!!! Warning MaruiChohokei::GetQuarkPosX !!!"<<endl;
		cout<<i<<" is invalid number that is larger than Num. of quark."<<endl;
		cout<<"*************************"<<endl;
		ret=0.;
		return ret;
	}else;

	return CharPos[i][0];
}

double MaruiChohokei::GetQuarkPosY(int i){
	double ret;
	if(i>NofQuarks){
		cout<<"*************************"<<endl;
		cout<<"!!! Warning MaruiChohokei::GetQuarkPosY !!!"<<endl;
		cout<<i<<" is invalid number that is larger than Num. of quark."<<endl;
		cout<<"*************************"<<endl;
		ret=0.;
		return ret;
	}else;

	return CharPos[i][1];
}


TString MaruiChohokei::GetParName(int i){
	TString ret;
	if(i==999){
		ret=ParName;
	}else{
		switch(i){
			case 0:
			ret="p";
			break;

			case 1:
			ret="n";
			break;

			case 2:
			ret="#Lambda";
			break;

			case 3:
			ret="K^{-}";
			break;

			case 4:
			ret="#it{#pi}^{-}";
			break;

			case 5:
			ret="#pi^{+}";
			break;

			case 6:
			ret="K^{+}";
			break;

			default:
			ret="#Lambda";
			break;
		}
	}

	return ret;
}


int MaruiChohokei::GetParCol(int i){
	int ret;
	if(i==999){
		ret=LColor;
	}else{
		switch(i){
			case 0:
			ret=2;
			break;

			case 1:
			ret=4;
			break;

			case 2:
			ret=6;
			break;

			case 3:
			ret=800;
			break;

			case 4:
			ret=841;
			break;

			case 5:
			ret=880;
			break;

			case 6:
			ret=3;
			break;

			default:
			ret=6;
			break;
		}
	}

	return ret;
}


void MaruiChohokei::SetTextSize(double TSize){Lat->SetTextSize(TSize);}
double MaruiChohokei::GetTextSize(){
	double ret;
	ret = Lat->GetTextSize();
	return ret;
}
