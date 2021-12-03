/* * * * * * * * * * * * 
 *  MaruiChohokei.h    *
 *                     *
 *  T. Fujiwara        *
 *  2021. 01. 10 (Sun) *
 * * * * * * * * * * * */
#ifndef MaruiChohokei_h
#define MaruiChohokei_h 1

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

#include "./Setting.h"

const int NofAngle = 4;


class MaruiChohokei
{
	public:
		 MaruiChohokei( double gx, double gy, int N, int pid );
		~MaruiChohokei();
		void SetAttribute( int LCol=602, int LSty=1, int LWid=1 );
		void SetQuarks();
		void Draw();
		void DrawQuarks();
		double  GetQPosX();
		double  GetGPosX();
		double  GetGPosY();
		double* GetGPos();
		double  GetQuarkPosX(int i=0);
		double  GetQuarkPosY(int i=0);
		TString GetParName(int i=999);
		int GetParCol(int i=999);
		void SetTextSize(double TSize=0.05);
		double GetTextSize();
		

	private:
		TLine* Ln[NofAngle];
		TEllipse* Ell[NofAngle];
		TArc* Ang[NofAngle];
		TLatex* Lat;

		Setting* MySetting;

		int ParticleID;
		int NofQuarks;

		int LColor = 602;
		int LStyle = 1;
		int LWidth = 4;

		double gPos[2];	//0:X, 1:Y
		double Hsize=0.0;
		double Vsize;

		double LnPos[NofAngle][2][2];
		double AngCenter[4][2];
		double Rad;
	//	double Rad_V;
		double PhiRange[NofAngle][2];
		double Theta[NofAngle];

		double VerticalRange;
		double CharPos[3][2];
		TString Quark[3];
		TString ParName;
};
#endif
