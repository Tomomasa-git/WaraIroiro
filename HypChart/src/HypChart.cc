/* * * * * * * * * * * * 
 * HypChart.cc         *
 *                     *
 *	T. Fujiwara        *
 *	2020. 10. 13 (Tue) *
 * * * * * * * * * * * */
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
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
#include "TPolyLine.h"
#include "TBox.h"
#include "TArc.h"
#include "TArrow.h"
#include "TPave.h"
#include "TPaveStats.h"
#include "TPaveText.h"

#include "./Setting.h"

#define ExamplePos_X 10.0
#define ExamplePos_Y 3.

typedef struct{
	int A;	//Mass number
	int Z;	//Atomic number (Num. of proton)
	int N;	//Neutron number
	string Name;
	TPave *Box;
	TPolyLine *Ang[3];
	TArc *Arc[2];
}Hyp;

const int NofReactionType = 3;

class LamHyp{
	public:
		 LamHyp(
		         string Name = "H",									//
		         int Mass    = 3,									//
		         int NP      = 1,									//
		         int NN      = 1,									//
		         int Kpi=1, int piK=1, int eeK=1, int Gamlay=0,
		         int Measure=1,
		         int Emulsion=1
		       );
		~LamHyp();
		void SetHypANum(int i);
		void SetHypZNum(int i);
		void SetHypNNum(int i);
		void SetHypName(string str){ HypName=str; }
		void SetNameTable();
		void GetNuclNameFromNT();
		void HypDraw();
		
	private:
		TPave *HypBox;
		TLatex *HypSymbol;
		TPolyLine *HypReaction[NofReactionType];
		TPolyLine *EmulsionMearsure;
		TGraph *GamLaySpec[2];
		Setting *Set;

		int A;	//Num. of Baryons
		int Z;	//Num. of Protons
		int N;	//Num. of Neutrons
		string HypName;

		int ReactFlag[NofReactionType] = {1, 1, 1};
		//0: Kpi
		//1: piK
		//2: ee'K

		int KpiFlag      = 1;
		int piKFlag      = 1;
		int eeKFlag      = 1;
		int GamlayFlag   = 1;
		int MeasureFlag  = 1;
		int EmulsionFlag = 1;
	
		double BoxGPos[2] = {1., 1.};
		double BoxEdgePos[4] = { 0.5, 0.5, 1.5, 1.5 };
		//                        |    |    |    |
		//                        |    |    |    Y2
		//                        |    |    X2
		//                        |    Y1
		//                        X1
		int BoxLCol   = 1;                        
		int BoxLSty   = 1;
		int BoxLWid   = 1;
		int BoxFCol   = 0;

		int AngCol[4] = {3, 2, 4, 18};	//Kpi, piK, eeK
		double GamlayPosX[2];
		double GamlayPosY[2];

//		double AngEdgeX[3][4];
//		double AngEdgeY[3][4];
		
//		int GamLayCol = 609;
		int GamLayCol = 602;
		double GamLaySize = 1.50;

		ifstream fin_NameTable;
		string NameTable;
		string NuclName[118];
};

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////

LamHyp::LamHyp(string Name, int Mass, int NP, int NN, int Kpi, int piK, int eeK, int Gamlay, int Measure, int Emulsion){
	Set = new Setting();

	HypName = Name;
	A = Mass;
	Z = NP;
	N = NN;

	if( A==777 ){
		HypName = "{}^{A}_{#Lambda}Z";
	}else if( (Z==0)&&(N==0) ){
		HypName = HypName;
	}else{
		HypName = Form("{}^{%d}_{#Lambda}",A)+HypName;
	}

	KpiFlag      = Kpi;
	piKFlag      = piK;
	eeKFlag      = eeK;
	GamlayFlag   = Gamlay;
	MeasureFlag  = Measure;
	EmulsionFlag = Emulsion;

	ReactFlag[0] = KpiFlag;
	ReactFlag[1] = piKFlag;
	ReactFlag[2] = eeKFlag;

	switch(MeasureFlag){
		case 1:
		BoxLSty = 1;
		BoxLWid = 1;
		break;

		case 0:
		BoxLSty = 2;
		BoxLWid = 1;
		break;

		case 2:
		BoxLCol = 800;
		BoxLSty = 3;
		BoxLWid = 1;
		break;

		default:
		BoxLSty = 1;
		BoxLWid = 1;
		break;
	}

	if( A<=18 ){
		BoxGPos[0] = (double)N;
		BoxGPos[1] = (double)Z;
	}else if(
		   (A>=19)
		&& (A<=28)
	){
		BoxGPos[0] = (double)N+0.5;
		BoxGPos[1] = (double)Z-13.+0.5;
	}else if(
		   (A>=29)
		&& (A<=33)
	){
		BoxGPos[0] = (double)N;
		BoxGPos[1] = (double)Z-13.;
	}else if(
		   (A>=40)
		&& (A< 41)
	){
		BoxGPos[0] = (double)N-2.5;
		BoxGPos[1] = (double)Z-13.-2.5;
	}else if(
		   (A>=41)
		&& (A< 49)
	){
		BoxGPos[0] = (double)N-9.0;
		BoxGPos[1] = (double)Z-13.-2.5;
	}else if(
		   (A>=50)
		&& (A< 53)
	){
		BoxGPos[0] = (double)N-7.5;
		BoxGPos[1] = (double)Z-13.-4.5;
	}else if(
		   (A>=56)
		&& (A< 57)
	){
		BoxGPos[0] = (double)N-7.5;
		BoxGPos[1] = (double)Z-13.-6.0;
	}else if(
		   (A>=88)
		&& (A< 90)
	){
		BoxGPos[0] = (double)N-26.0;
		BoxGPos[1] = (double)Z-13.-17.5;
	}else if(
		   (A>=120)
		&& (A< 140)
	){
		BoxGPos[0] = (double)N-56.;
		BoxGPos[1] = (double)Z-13.-34.5;
	}else if(
		   (A>=200)
		&& (A< 210)
	){
		BoxGPos[0] = (double)N-95.5;
		BoxGPos[1] = (double)Z-13.-57.0;
	}else if( A==777 ){
		BoxGPos[0] = ExamplePos_X;
		BoxGPos[1] = ExamplePos_Y;
	}else{
		BoxGPos[0] = (double)N;
		BoxGPos[1] = (double)Z;
	}
	
	BoxEdgePos[0] = BoxGPos[0]-0.5;
	BoxEdgePos[1] = BoxGPos[1]-0.5;
	BoxEdgePos[2] = BoxGPos[0]+0.5;
	BoxEdgePos[3] = BoxGPos[1]+0.5;
	if( A==777 ){
		BoxEdgePos[0] = BoxGPos[0]-1.;
		BoxEdgePos[1] = BoxGPos[1]-1.;
		BoxEdgePos[2] = BoxGPos[0]+1.;
		BoxEdgePos[3] = BoxGPos[1]+1.;
	}else;

	double AngEdgeX[4][4]= { {BoxEdgePos[0], BoxEdgePos[0]+0.25, BoxEdgePos[0]     , BoxEdgePos[0]}, 
	           	             {BoxEdgePos[0], BoxEdgePos[0]+0.25, BoxEdgePos[0]     , BoxEdgePos[0]},
	                         {BoxEdgePos[2], BoxEdgePos[2]-0.25, BoxEdgePos[2]     , BoxEdgePos[2]},
	                         {BoxEdgePos[2], BoxEdgePos[2]-0.25, BoxEdgePos[2]     , BoxEdgePos[2]} };
	double AngEdgeY[4][4]= { {BoxEdgePos[1], BoxEdgePos[1]     , BoxEdgePos[1]+0.25, BoxEdgePos[1]}, 
	                         {BoxEdgePos[3], BoxEdgePos[3]     , BoxEdgePos[3]-0.25, BoxEdgePos[3]},
	                         {BoxEdgePos[3], BoxEdgePos[3]     , BoxEdgePos[3]-0.25, BoxEdgePos[3]},
	                         {BoxEdgePos[1], BoxEdgePos[1]     , BoxEdgePos[1]+0.25, BoxEdgePos[1]} };
	if( A==777 ){
		AngEdgeX[0][1] = BoxEdgePos[0]+0.50;
		AngEdgeX[1][1] = BoxEdgePos[0]+0.50;
		AngEdgeX[2][1] = BoxEdgePos[2]-0.50;
		AngEdgeX[3][1] = BoxEdgePos[2]-0.50;
		
		AngEdgeY[0][2] = BoxEdgePos[1]+0.50;
		AngEdgeY[1][2] = BoxEdgePos[3]-0.50;
		AngEdgeY[2][2] = BoxEdgePos[3]-0.50;
		AngEdgeY[3][2] = BoxEdgePos[1]+0.50;
	}else;

	GamlayPosX[0] = ( AngEdgeX[0][0]+AngEdgeX[0][1]+AngEdgeX[0][2] )/3.;
	GamlayPosX[1] = ( AngEdgeX[1][0]+AngEdgeX[1][1]+AngEdgeX[1][2] )/3.;
	GamlayPosY[0] = ( AngEdgeY[0][0]+AngEdgeY[0][1]+AngEdgeY[0][2] )/3.;
	GamlayPosY[1] = ( AngEdgeY[1][0]+AngEdgeY[1][1]+AngEdgeY[1][2] )/3.;

	HypBox = new TPave( BoxEdgePos[0], BoxEdgePos[1],BoxEdgePos[2], BoxEdgePos[3], 1, "NB" );
	Set -> Setting_Pave( HypBox, BoxLCol, BoxLWid, BoxLSty, BoxFCol );
//	HypBox -> SetBorderSize(1);

	HypSymbol = new TLatex();
	Set -> Setting_Latex( HypSymbol, 42, 22, 1, 0.0300 );
	if( A==777 ){
		HypSymbol->SetTextSize(0.060);
	}else;

	for(int i=0; i<NofReactionType; i++){
		HypReaction[i] = new TPolyLine(4, AngEdgeX[i], AngEdgeY[i]);
		Set -> Setting_Line( HypReaction[i], AngCol[i], 1, 1, AngCol[i], 1001);
	}

	for(int i=0; i<2; i++){
		GamLaySpec[i] = new TGraph();
		Set -> Setting_Graph( GamLaySpec[i], "", "", "", 1, 1, 42, GamLayCol, 20 );
		GamLaySpec[i] -> SetPoint( 0, GamlayPosX[i], GamlayPosY[i] );
		if( A==777 ){
			GamLaySpec[i] -> SetMarkerSize(1.50*GamLaySize);
		}else{
			GamLaySpec[i] -> SetMarkerSize(GamLaySize);
		}
	}

	EmulsionMearsure = new TPolyLine();
	EmulsionMearsure -> SetLineColor(AngCol[3]);
	EmulsionMearsure -> SetLineStyle(1);
	EmulsionMearsure -> SetLineWidth(1);
	EmulsionMearsure -> SetFillColor(AngCol[3]);
	EmulsionMearsure -> SetFillStyle(1001);
	for(int i=0; i<4; i++){EmulsionMearsure->SetPoint(i, AngEdgeX[3][i], AngEdgeY[3][i]);}
	
}

/*--------------------------------------------------------------------------------------------------------*/
LamHyp::~LamHyp(){
	cout<<"aaaaaaa"<<endl;
}

/*--------------------------------------------------------------------------------------------------------*/
void LamHyp::SetNameTable(){
	string str="";
	int i=0;

	NameTable = "../dat/NameTable00.dat";
	fin_NameTable.open( NameTable.c_str(), ios_base::in );
	
	if(!fin_NameTable){
		cout<<"!!!WARNING!!!   LamHyp::SetNameTable   !!!WARNING!!!"<<endl;
		cout<<"There is no such file: "<<NameTable<<endl;
		return;
	}else;

	while(fin_NameTable>>str){
		NuclName[i] = str;
		i++;
	}	
}

/*--------------------------------------------------------------------------------------------------------*/
void LamHyp::GetNuclNameFromNT(){
	if(Z>0){
		HypName = Form("{}^{%d}_{#Lambda}",A)+NuclName[Z-1];
	}else if(Z==0){
		HypName = "#Lambda";
	}else{
		HypName = "#Lambda";
	}
}

/*--------------------------------------------------------------------------------------------------------*/
void LamHyp::HypDraw(){
	HypBox->Draw();

	HypSymbol -> DrawLatex( BoxGPos[0], BoxGPos[1], HypName.c_str() );

	for(int i=0; i<NofReactionType; i++){
		if(ReactFlag[i]==1){
			HypReaction[i]->Draw("f");
			HypReaction[i]->Draw();
		}else;
	}

	if(EmulsionFlag==1){
		EmulsionMearsure->Draw("f");
		EmulsionMearsure->Draw();
	}else;

	switch(GamlayFlag){
		case 0:
		break;

		case 1:
		GamLaySpec[0]->Draw("sameP");
		break;

		case 2:
		GamLaySpec[0]->Draw("sameP");
		GamLaySpec[1]->Draw("sameP");
		break;

		default:
		break;
	}
}

/*--------------------------------------------------------------------------------------------------------*/
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////

int main( int argc, char** argv){
	const int NofHypNcl = 43;

	int CanvH = 1600;
//	int CanvV = 1000;
	int CanvV = 1250;

//	double Hmin = -0.75;
//	double Hmax = 32.50;
//	double Vmin = -0.75;
//	double Vmax = 16.50;
	double Hmin = -0.75;
	double Hmax = 14.50;
	double Vmin = -0.75;
//	double Vmax = 9.50;
	double Vmax = 10.00;

	ifstream fin_ParamTable;
	string ParamTableName = "../dat/ParamTable00.dat";
	int I=0;
	int val[9];
	int TableA[NofHypNcl];
	int TableZ[NofHypNcl];
	int TableN[NofHypNcl];
	int TableKpi[NofHypNcl];
	int TablepiK[NofHypNcl];
	int TableeeK[NofHypNcl];
	int TableGamLay[NofHypNcl];
	int TableEmulsion[NofHypNcl];
	int TableMeasure[NofHypNcl];

	TApplication *theApp;
	TCanvas *Ca;
	TH1D *h_fr;
	TArrow *Ax[2];
	//0: Z (Proton num.)	
	//1: N (Neutron num.)	
	TF1 *f_Stable;
	TGraph *g_Stable;
	TLine *L_Connect;
	TF1 *f_StableHeavier;
	TGraph *g_StableHeavier;
	TLatex *Lat;
	TLatex *Lat_tit;

	LamHyp *Lambda;
	LamHyp *nnL;
	LamHyp *HypNcl[NofHypNcl];
	LamHyp *HypExample;
	TLatex *Lat_ind;
	TLine *Ln_ind[5];

	theApp = new TApplication( "App", &argc, argv );

	fin_ParamTable.open( ParamTableName.c_str(), ios_base::in );
	
	while( fin_ParamTable>>val[0]>>val[1]>>val[2]>>val[3]>>val[4]>>val[5]>>val[6]>>val[7]>>val[8] ){
		TableA[I]        = val[0];	cout<<TableA[I]<<" ";
		TableZ[I]        = val[1];	cout<<TableZ[I]<<" ";
		TableN[I]        = val[2];	cout<<TableN[I]<<" ";
		TableKpi[I]      = val[3];	cout<<TableKpi[I]<<" ";
		TablepiK[I]      = val[4];	cout<<TablepiK[I]<<" ";
		TableeeK[I]      = val[5];	cout<<TableeeK[I]<<" ";
		TableGamLay[I]   = val[6];	cout<<TableGamLay[I]<<" ";
		TableMeasure[I]  = val[7];	cout<<TableMeasure[I]<<" ";
		TableEmulsion[I] = val[8];	cout<<TableEmulsion[I]<<endl;
		I++;
	}

	gROOT->SetBatch(1);
	gStyle->SetFrameLineWidth(0);
	
	Ca = new TCanvas( "Ca", "Ca", CanvH, CanvV );
	Ca -> SetMargin( .075, .025, .075, .100 );

	h_fr = new TH1D( "h_fr", "h_fr", 1, Hmin, Hmax );
	h_fr -> SetStats(0);
	h_fr -> SetTitle("");
	h_fr -> SetLineColor(0);
	h_fr -> GetXaxis()->SetNdivisions(0);
	h_fr -> GetYaxis()->SetNdivisions(0);
	h_fr -> GetYaxis()->SetRangeUser( Vmin, Vmax);

	Ax[0] = new TArrow( Hmin, Vmin, Hmin, Vmax );
	Ax[1] = new TArrow( Hmin, Vmin, Hmax, Vmin );
	for(int i=0; i<2; i++){
		Ax[i] -> SetArrowSize(0.025);
		Ax[i] -> SetAngle(30.);
	}
	
	Lat = new TLatex();
	Lat -> SetTextFont(62);
	Lat -> SetTextAlign(22);

	Lat_tit = new TLatex();
	Lat_tit -> SetTextFont(72);
	Lat_tit -> SetTextAlign(22);
	Lat_tit -> SetTextColor(602);
	Lat_tit -> SetTextSize(.045);

	f_Stable = new TF1( "f_Stable", "x", Hmin, 10. );
	g_Stable = new TGraph(f_Stable);
	g_Stable -> SetLineColor(602);
	g_Stable -> SetLineWidth(1);

	f_StableHeavier = new TF1( "f_StableHeavier", "x+[0]", 12.5, Hmax-2. );
	f_StableHeavier -> SetParameter(0, -13.);
	g_StableHeavier = new TGraph(f_StableHeavier);
	g_StableHeavier -> SetLineColor(602);
	g_StableHeavier -> SetLineWidth(1);

	L_Connect = new TLine( 12., 12., 12.5, -0.5 );
	L_Connect -> SetLineColor(602);
	L_Connect -> SetLineWidth(1);
	L_Connect -> SetLineStyle(4);

	h_fr            -> Draw("");
	g_Stable        -> Draw("same, L");
//	L_Connect       -> Draw();
//	g_StableHeavier -> Draw("same, L");


	HypExample = new LamHyp("", 777, 0, 0, 1, 1, 1, 1, 1, 1 );
	HypExample -> HypDraw();
	Ln_ind[0] = new TLine( ExamplePos_X-1.0 , ExamplePos_Y-1.0  , ExamplePos_X-1.5 , ExamplePos_Y-1.5  );
	Ln_ind[1] = new TLine( ExamplePos_X-1.0 , ExamplePos_Y+1.0  , ExamplePos_X-1.5 , ExamplePos_Y+1.5  );
	Ln_ind[2] = new TLine( ExamplePos_X+1.0 , ExamplePos_Y+1.0  , ExamplePos_X+1.5 , ExamplePos_Y+1.5  );
	Ln_ind[3] = new TLine( ExamplePos_X+1.0 , ExamplePos_Y-1.0  , ExamplePos_X+1.5 , ExamplePos_Y-1.5  );
	Ln_ind[4] = new TLine( ExamplePos_X-0.83, ExamplePos_Y-0.925, ExamplePos_X-0.83, ExamplePos_Y-2.25 );
	for(int i=0; i<5; i++){Ln_ind[i] -> Draw();}
	Lat_ind = new TLatex();
	Lat_ind -> SetTextAlign(22);
	Lat_ind -> SetTextFont(42);
	Lat_ind -> SetTextSize(.040);
	Lat_ind -> DrawLatex( ExamplePos_X, ExamplePos_Y+2.50, "Studied by");
	Lat_ind -> SetTextSize(.030);
	Lat_ind -> DrawLatex( ExamplePos_X-1.75, ExamplePos_Y-1.75, "(#it{K^{-}}, #pi^{-})");
	Lat_ind -> DrawLatex( ExamplePos_X-1.75, ExamplePos_Y+1.75, "(#pi^{+}, #it{K^{+}})");
	Lat_ind -> DrawLatex( ExamplePos_X+1.75, ExamplePos_Y+1.75, "(#it{e}, #it{e'K^{+}})");
	Lat_ind -> DrawLatex( ExamplePos_X+1.75, ExamplePos_Y-1.75, "Emulsion");
	Lat_ind -> DrawLatex( ExamplePos_X-0.83, ExamplePos_Y-2.50, "#gamma-ray spectroscopy");

	Lambda = new LamHyp( "#Lambda", 1, 0, 0, 0, 0, 0, 0, 1, 0 );
	Lambda -> HypDraw();

	nnL = new LamHyp( "n", 3, 0, 2, 0, 0, 0, 0, 0, 0);
	nnL -> HypDraw();
	
//	for(int i=0; i<NofHypNcl; i++){
	for(int i=0; i<28; i++){
		HypNcl[i] = new LamHyp( "", 
		                        TableA[i], TableZ[i], TableN[i], 
		                        TableKpi[i], TablepiK[i], TableeeK[i], TableGamLay[i], 
		                        TableMeasure[i], TableEmulsion[i]
		                      );
		HypNcl[i] -> SetNameTable();
		HypNcl[i] -> GetNuclNameFromNT();
		HypNcl[i] -> HypDraw();
	}

	for(int i=0; i<2; i++){Ax[i] -> Draw("same, |>");}
//	Lat -> SetTextFont(72);
	Lat -> DrawLatexNDC(.80, .04, "N #scale[0.75]{(Neutron Number)}");
	Lat -> SetTextAngle(90.);
	Lat -> DrawLatexNDC(.035, .80, "Z #scale[0.75]{(Proton Number)}");
	Lat -> SetTextAngle(45.);
	Lat -> SetTextColor(602);
	Lat -> SetTextSize(.040);
	Lat -> DrawLatex(9., 9.5, "Z=N");

	Lat_tit -> DrawLatexNDC( .350, .950, "#Lambda hypernuclear chart (2020)" );
	Lat_tit -> SetTextSize(.025);
	Lat_tit -> SetTextFont(42);
	Lat_tit -> DrawLatexNDC( .350, .900, "(Only single #Lambda and light mass region)");
	

	Ca->Print("../fig/test_05.pdf", "pdf");	

	//Test for name table
//	ifstream ifs( "../dat/NameTable00.dat", ios_base::in );
//	string val;
//	int i=1;
//	while(ifs>>val){
//		cout<<Form("%03d", i)<<"  "<<val<<endl;
//		i++;
//	}

	gSystem->Exit(-1);

	theApp->Run();
	return 0;
}
