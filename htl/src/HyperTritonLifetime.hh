/* * * * * * * * * * * * * * * 
 *  HyperTritonLifitime.hh   *
 *                           *
 *  T. Fujiwara              *
 *  2020. 11. 27 (Fri)       *
 * * * * * * * * * * * * * * */
#ifndef HyperTritonLifitime_h
#define HyperTritonLifitime_h 1

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
#include "TPolyLine.h"
#include "TBox.h"
#include "TPave.h"
#include "TPaveStats.h"
#include "TPaveText.h"

#include "./Setting.h"

const int NofData = 99;

class HyperTritonLifetime
{
	public:
		 HyperTritonLifetime( int ID=0 );
		~HyperTritonLifetime();
		void LTDraw();
		void SetDT();
		void GetDT();
		int SetMCol(int Type=0);
		int SetMSty(int ID=0, int Type=0);
		void SetDataPointFromDT();
	
	private:
		TGraphAsymmErrors* gr_lt;
		TBox* SystErr;
		TLatex* Lat;
		TPaveText* Pt;
		Setting* Set;
		
		int DataID;
		int DataType;
		double Xpos;
		double BoxEdgePos[4];
			//0: X1
			//1: X2
			//2: Y1
			//3: Y2
		double Lifetime_Val;
		double Lifetime_StatErr[2];	//0: Upper, 1: Lower
		double Lifetime_SystErr[2];	//0: Upper, 1: Lower
		double RefPos[2];			//0: Xpos,  1: Ypos
		string RefText;
		string RefPaveText[3];

		int Col;
		int Sty;

		string datamname;
		ifstream ifs;
		int LT_ID[NofData];
		int LT_DataType[NofData];
		double LT_Val[NofData];
		double LT_StatErr[NofData][2];
		double LT_SystErr[NofData][2];
		string LT_Ref[NofData][6];

		double BoxWidth =0.20;
		double BoxHeight[2];
		
};

#endif
