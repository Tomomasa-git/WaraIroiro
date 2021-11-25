/* * * * * * * * * * * * * * * 
 *  HyperTritonLifitime.cc   *
 *                           *
 *  T. Fujiwara              *
 *  2020. 11. 27 (Fri)       *
 * * * * * * * * * * * * * * */
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

HyperTritonLifetime::HyperTritonLifetime(int ID){
	Set = new Setting();

	gr_lt = new TGraphAsymmErrors();
	SystErr = new TBox();
	Lat = new TLatex();
	Pt = new TPaveText( .0, .0, 1., 1. );
	Pt -> SetFillColor(0);
	Pt -> SetFillStyle(0);
	Pt -> SetLineColor(0);
	Pt -> SetLineStyle(0);
	Pt -> SetLineWidth(0);
	Pt -> SetTextColor(1);
	Pt -> SetTextFont(42);
//	Pt -> SetTextSize(.0325);
	Pt -> SetTextSize(.0400);
	Pt -> SetBorderSize(0);

	DataID = ID;
}

HyperTritonLifetime::~HyperTritonLifetime(){
	delete Set;
}


void HyperTritonLifetime::LTDraw(){
	if(DataType==2)gr_lt->SetMarkerSize(2.5);
	if(
		   DataID==6
		|| DataID==9
	){
		gr_lt->SetMarkerSize(3.50);
	}else if( DataID==7 ){
		gr_lt->SetMarkerSize(3.00);
	}else;
	gr_lt->Draw("samePL");
	if(DataType==2){
		SystErr->Draw();
	}

	RefPos[0] = Xpos;
	if(DataType==2){	//Heavy ion Collision
		if( (DataID%2) ){
			RefPos[1] = Lifetime_Val-(Lifetime_StatErr[0]+50.);
		}else{
			RefPos[1] = Lifetime_Val+(Lifetime_StatErr[1]+20.);
		}
	}else if( (DataID%2) ){
		RefPos[1] = Lifetime_Val-(Lifetime_StatErr[0]+25.);
	}else{
		RefPos[1] = Lifetime_Val+(Lifetime_StatErr[1]+25.);
	}

	if( DataID==1 )RefPos[1] = Lifetime_Val-(Lifetime_StatErr[0]+25.);
	if( DataID==4 )RefPos[1] = Lifetime_Val-(Lifetime_StatErr[0]+25.);
	if( DataID==5 )RefPos[1] = Lifetime_Val+(Lifetime_StatErr[1]+25.);
	if( DataID==6 )RefPos[1] = Lifetime_Val+(Lifetime_StatErr[1]+30.);
	if( DataID==8 )RefPos[1] = Lifetime_Val+(Lifetime_StatErr[1]+60.);
	if( DataID==10)RefPos[1] = Lifetime_Val+(Lifetime_StatErr[1]+30.);
	if( DataID==11)RefPos[1] = Lifetime_Val-(Lifetime_StatErr[0]+60.);

	Pt -> SetX1(RefPos[0]-0.4);
	Pt -> SetX2(RefPos[0]+0.4);
	Pt -> SetY1(RefPos[1]-30.);
	Pt -> SetY2(RefPos[1]+30.);
	if(DataType==2){
		Pt -> SetY2(RefPos[1]+40.);
//		Pt -> SetTextSize(.02);
	}else;
	Pt -> Draw();
//	Lat->DrawLatex( RefPos[0], RefPos[1], RefText.c_str() );
}


void HyperTritonLifetime::SetDT(){
	int val00    = 0;
	int val01    = 0; 
	double val02 = 0.;
	double val03 = 0.;
	double val04 = 0.;
	double val05 = 0.;
	double val06 = 0.;
	string val07 = "";
	string val08 = "";
	string val09 = "";
	int val10 = 0;
	int val11 = 0;
	int Count=0;

	datamname = "../dat/htl_datatable.dat";
	ifs.open( datamname.c_str(), ios_base::in );
	while(
		ifs>>val00>>val01>>val02>>val03>>val04>>val05>>val06>>val07>>val08>>val09>>val10>>val11
	){
		LT_ID[Count]         = val00;
		LT_DataType[Count]   = val01;
		LT_Val[Count]        = val02;
		LT_StatErr[Count][0] = val03;
		LT_StatErr[Count][1] = val04;
		LT_SystErr[Count][0] = val05;
		LT_SystErr[Count][1] = val06;
		LT_Ref[Count][0]     = val07;
		LT_Ref[Count][1]     = val08;
		LT_Ref[Count][2]     = val09;
		LT_Ref[Count][3]     = Form("(%d)", val10);
		if( Count==9 ){
			LT_Ref[Count][4] = Form("0%d", val11);
		}else{
			LT_Ref[Count][4] = Form("%d", val11);
		}

		Count ++;
	}
}
/*
#ID Type Lifetime StatErr.Low StatErr.Upp SystErr.Low SystErr.Upp Ref1, Ref2, Ref3, Ref4
Type:
	0>>Emulsion
	1>>Bubble Chamber
	2>>Heavy Ion Collision

*/

int HyperTritonLifetime::SetMCol(int Type){
	int ret;

	switch(Type){
		//Emulsion
		case 0:
		ret = 4;
		break;

		//Bubble chamber
		case 1:
		ret = 3;
		break;

		//Heavy ion collision
		case 2:
		ret = 2;
		break;

		default:
		ret = 4;
		break;
	}

	return ret;
}

int HyperTritonLifetime::SetMSty(int ID, int Type){
	int ret;

	if(Type==0){		//Emulsion
		ret=23;
	}else if(Type==1){	//Bubble chamber
		ret=20;
	}else{
		if( (ID==6)||(ID==9) ){			//STAR
			ret=29;
		}else if(ID==7){				//HypHI
			ret=33;
		}else if( (ID==8)||(ID==10)||(ID==11) ){	//ALICE
			ret=34;
		}else ret=20;
	}

	return ret;
}


void HyperTritonLifetime::SetDataPointFromDT(){
	DataType = LT_DataType[DataID];

	RefText="";
//	for(int i=0; i<5; i++)RefText=RefText+LT_Ref[DataID][i];
//	RefText = "#splitline{#splitline{"+LT_Ref[DataID][0]+"}{"+LT_Ref[DataID][1]+" #bf{"+LT_Ref[DataID][2]+"}}}{"+LT_Ref[DataID][3]+" "+LT_Ref[DataID][4]+"}";
	if( (DataType==0)||(DataType==1) ){
		RefPaveText[0] = LT_Ref[DataID][0];									cout<<RefPaveText[0]<<endl;
		RefPaveText[1] = LT_Ref[DataID][1]+" "+LT_Ref[DataID][2];			cout<<RefPaveText[1]<<endl;
		RefPaveText[2] = "";
	}else if(DataType==2){
		RefPaveText[0] = LT_Ref[DataID][0];									cout<<RefPaveText[0]<<endl;
		RefPaveText[1] = LT_Ref[DataID][1]+" #bf{"+LT_Ref[DataID][2]+"}";	cout<<RefPaveText[1]<<endl;
		RefPaveText[2] = LT_Ref[DataID][3]+" "+LT_Ref[DataID][4];			cout<<RefPaveText[2]<<endl;
	}else{
		RefPaveText[0] = LT_Ref[DataID][0];									cout<<RefPaveText[0]<<endl;
		RefPaveText[1] = LT_Ref[DataID][1]+" "+LT_Ref[DataID][2];			cout<<RefPaveText[1]<<endl;
	}
//	cout<<RefText<<endl;
	if( (DataType==0)||(DataType==1) ){
		Pt -> AddText(RefPaveText[0].c_str());
		Pt -> AddText(RefPaveText[1].c_str());
	}else{
		Pt -> AddText(RefPaveText[0].c_str());
		Pt -> AddText(RefPaveText[1].c_str());
		Pt -> AddText(RefPaveText[2].c_str());
	}

	Lifetime_Val = LT_Val[DataID];
	for(int i=0; i<2; i++){
		Lifetime_StatErr[i]=LT_StatErr[DataID][i];
		Lifetime_SystErr[i]=LT_SystErr[DataID][i];
	}

	if(DataType==2){
	//	Xpos = (double)DataID+2.;
		Xpos = (double)DataID+2. + 0.30*( (double)DataID-5.);
	}else{
		Xpos = (double)DataID+1.;
	}

	for(int i=0; i<2; i++){BoxHeight[i] = Lifetime_Val+pow(-1.,i+1)*Lifetime_SystErr[i];}
	BoxEdgePos[0] = Xpos-BoxWidth;
	BoxEdgePos[1] = Xpos+BoxWidth;
	BoxEdgePos[2] = BoxHeight[0];
	BoxEdgePos[3] = BoxHeight[1];

	gr_lt -> SetPoint( 0, Xpos, Lifetime_Val );
	gr_lt -> SetPointError( 0, 0., 0., Lifetime_StatErr[0], Lifetime_StatErr[1] );
	Col=SetMCol(DataType);
	Sty=SetMSty(DataID, DataType);
	Set -> Setting_GError( gr_lt, "", "", "", Col, 1, 42, Col, Sty, 2.50 );
	gr_lt -> SetLineWidth(1);

	SystErr -> SetX1(BoxEdgePos[0]);
	SystErr -> SetX2(BoxEdgePos[1]);
	SystErr -> SetY1(BoxEdgePos[2]);
	SystErr -> SetY2(BoxEdgePos[3]);
	Set -> Setting_Box(SystErr, Col);
	
//	Set -> Setting_Latex( Lat, 42, 22, 1, .025 );
}

