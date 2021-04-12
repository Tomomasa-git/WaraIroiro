#ifndef Setting_h
#define Setting_h

/*
	2020. 07. 08 (Wed)
	Analysis for ToF cosmic
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
#include "TAxis.h"
#include "TGaxis.h"
#include "TMath.h"
#include "TSpectrum.h"
#include "TString.h"
#include "TPaveStats.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TPolyLine.h"
#include "TBox.h"
#include "TArc.h"
#include "TArrow.h"

class Setting
{
	public:
		Setting();
		~Setting();
		void Setting_Gene  ( int BatchFlag );
		void Setting_Hist1D( TH1D *hist      , TString HTitle, TString XTitle, TString YTitle, int LCol, int LSty, int Font, int FCol, int FSty );
		void Setting_Hist2D( TH2D *hist      , TString HTitle, TString XTitle, TString YTitle, TString ZTitle, double Min );
		void Setting_Graph ( TGraph *gr      , TString GTitle, TString XTitle, TString YTitle, int LCol, int LSty, int Font, int MCol, int MSty );
		void Setting_GError( TGraphErrors *gr, TString GTitle, TString XTitle, TString YTitle, int LCol, int LSty, int Font, int MCol, int MSty );
		void Setting_Func  ( TF1 *func , int LCol, int LSty );
		void Setting_Legend( TLegend *leg   , int Font, int Al, int Col, double Size );
		void Setting_Latex ( TLatex *lat    , int Font, int Al, int Col, double Size );
		void Setting_Line  ( TLine *lin     , int LCol, int LWid, int LSty );
		void Setting_Line  ( TPolyLine *Poll, int LCol, int LWid, int LSty, int FCol, int FSty );
		void Setting_Pave  ( TPave *pave    , int LCol, int LWid, int LSty, int FCol );
		void Setting_Box   ( TBox *box      , int LCol, int LWid, int LSty, int FCol );
};

#endif
