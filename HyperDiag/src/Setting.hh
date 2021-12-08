#ifndef Setting_h
#define Setting_h 1

/*
	2020. 12. 11 (Fri)
	For FToF Optical sim.
 */

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
#include "TCurlyLine.h"
#include "TArrow.h"
#include "TEllipse.h"
#include "TArc.h"
#include "TBox.h"
#include "TPave.h"
#include "TPaveStats.h"
#include "TPaveText.h"


/////////////////////////////////////////////////////////////////////////////////////
class Setting {
	public:
		Setting();
		~Setting();
		void Setting_Gene  ( int BatchFlag=1 );

		void Setting_Hist1D( 
			TH1D *hist, 
			TString HTitle="", 
			TString XTitle="", 
			TString YTitle="", 
			int LCol=602, 
			int LSty=1, 
			int Font=42, 
			int FCol=422, 
			int FSty=1001 
		);

		void Setting_Hist2D( 
			TH2D *hist, 
			TString HTitle="", 
			TString XTitle="", 
			TString YTitle="", 
			TString ZTitle="", 
			double Min=0.95 
		);


		void Setting_Graph ( 
			TGraph *gr, 
			TString GTitle="", 
			TString XTitle="", 
			TString YTitle="", 
			int LCol=602, 
			int LSty=1, 
			int Font=42, 
			int MCol=6, 
			int MSty=29, 
			double MSiz=2. 
		);

		void Setting_GError( 
			TGraphErrors *gr, 
			TString GTitle="", 
			TString XTitle="", 
			TString YTitle="", 
			int LCol=602, 
			int LSty=1, 
			int Font=42, 
			int MCol=6, 
			int MSty=29, 
			double MSiz=2. 
		);

		void Setting_GError( 
			TGraphAsymmErrors *gr, 
			TString GTitle="", 
			TString XTitle="", 
			TString YTitle="", 
			int LCol=602, 
			int LSty=1, 
			int Font=42, 
			int MCol=6, 
			int MSty=29, 
			double MSiz=2. 
		);

		void Setting_Func  ( TF1 *func , int LCol=4, int LSty=1 );

		void Setting_Legend( TLegend* leg, int Font=42, int Al=22, int Col=602, double Size=0.05, int TranspFlag=0 );
		void Setting_Latex ( TLatex *lat , int Font=42, int Al=22, int Col=602, double Size=.05 );

		void Setting_Line  ( TLine *lin     , int LCol=602, int LWid=1, int LSty=1 );
		void Setting_Line  ( TPolyLine *Poll, int LCol=602, int LWid=1, int LSty=1, int FCol=0, int FSty=0 );
		void Setting_Line  ( TCurlyLine* lin, int LCol=602, int LWid=1, int LSty=1, double CulAmp=0.0125, double CulWL=0.01 );
		void Setting_Arrow ( TArrow* arr    , int LCol=602, int LWid=1, int LSty=1 );
		void Setting_Arc   ( TArc* arc      , int LCol=602, int LWid=1, int LSty=1, int FCol=0, int FSty=0 );
		void Setting_Box   ( TBox *box      , int LCol=602, int Wid=1, int Sty=1, int FCol=0, int FSty=0 );

		void Setting_Pave( 
			TPaveText* Pt,
			int Font=42, int TCol=602, int TAli=22, double TSize=0.05,	//Text Attributes
			int LCol=1, int LWid=1, int LSty=1,							//Line Attributes
			int FCol=0, int FSty=0,										//Fill Attributes
			int BSize=4
		);

		void RedrawFrame(TFrame *frame);
	
		TString DrawOption( int Num=0, int Type=0 );
};

#endif
