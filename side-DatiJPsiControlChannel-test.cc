//
//g++ -O3 -o side-DatiJPsiControlChannel-test side-DatiJPsiControlChannel-test.cc `root-config --cflags --libs`  -lRooFit -lGenVector
//
// 
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <sys/stat.h>
#include "Riostream.h"
#include <map>
#include <string>
//#include <boost/algorithm/string.hpp>
//#include <boost/algorithm/string/trim.hpp>
#include <vector>
#include <math.h>
//#include <TCint.h>
//#include <TGenericClassInfo.h> 
#include <TMath.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TF1.h>
#include <TF2.h>
#include <TF3.h>
#include <TROOT.h>
#include <TEnv.h>
#include <TSystem.h>
#include <TTree.h>
#include "TBranch.h"
#include <TApplication.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TStyle.h> 
#include <TF1.h>  
#include <TF2.h> 
#include <TLorentzVector.h>
#include <TVector3.h>
#include "TDSet.h"
#include "TChain.h"
#include <time.h> 
#include <TSystemDirectory.h>
#include <TEfficiency.h>
#include <TLegend.h>
#include <TMinuit.h>
#include "Math/WrappedMultiTF1.h"
#include "TRandom.h" 
#include "TRandom3.h" 
#include  <TStopwatch.h>
#include "TH1F.h"
#include "TH2F.h"			// unused?
#include "TStyle.h"
#include "TCanvas.h"
#include <TGraphAsymmErrors.h>
#include <TFrame.h>
#include <TFitResult.h>
#include <TFitter.h>
#include "Fit/Fitter.h"
#include <TMatrixDSym.h>
#include <TBinomialEfficiencyFitter.h>
#include <TKDTreeBinning.h>
#include <TH2Poly.h>
#include <Math/Vector4Dfwd.h>
#include <Math/GenVector/LorentzVector.h>
#include <Math/GenVector/PtEtaPhiM4D.h>//#include <RooMath/GenVector/PtEtaPhiM4D.h>
//#if !defined(__CINT__) || defined(__MAKECINT__)
#include <RooRandom.h>
#include <RooFit.h>
#include <RooMinuit.h>
#include <RooRealVar.h>
#include <RooConstVar.h>
#include <RooPlot.h>
#include <RooDataSet.h>
#include <RooDataHist.h>
#include <RooHistPdf.h>
#include <RooRealVar.h>
#include <RooAbsReal.h>
#include <RooArgSet.h>
#include <RooRealProxy.h>
#include <RooCategoryProxy.h>
#include <RooAbsCategory.h>
#include <RooProdPdf.h>
#include <RooAddPdf.h>
#include <RooProduct.h>
#include <RooGaussian.h>
#include <RooCBShape.h>
#include <RooPolynomial.h>
#include <RooChebychev.h>
#include <RooWorkspace.h>
#include <RooExponential.h>
#include <RooErrorVar.h>
#include <RooFitResult.h>
#include <RooRangeBinning.h>
#include <RooBinning.h>
#include <RooNumGenConfig.h>
#include <RooBernstein.h>
#include <RooPolynomial.h>
#include <RooExtendPdf.h>
#include "RooFoamGenerator.h"
#include "RooAcceptReject.h"
#include "GBRMath.h"
#include "RooDoubleCBFast.h"
#include "RooBernsteinSideband.h"
#include "RooMCStudy.h"
#include "TFoam.h"
#include "TRatioPlot.h"

#include <sys/time.h>
#include <sys/times.h>
#include <iostream>
// System stuff
#include <fstream> 
#include <sys/time.h>
#include <sys/times.h>

timeval startTime, stopTime, totalTime;
clock_t startCPU, stopCPU; 
tms startProc, stopProc; 

using namespace std; 
using namespace ROOT;
using namespace RooFit;
void AnalizeDatiJPsi_test();
std::string grep(std::string& filename, std::string keyword);
//
TFile*OutFileInputHisto;

//
char * RunEra;
char OutFileNameInputHisto[300]        = "/gwpool/users/dini/p5prime/2018/skims/newphi/2018Data_All_finalSelection.root";
char OutFileNameInputHistoMC[300]      = "/gwpool/users/dini/p5prime/2018/skims/newphi/2018MC_JPSI.root";
char OutFileNameInputHistoMCJPsi[300] = "/gwpool/users/dini/p5prime/sidebands/2018MC_BToPsi2SK.root";

char OutputRecoB0TreeName[10]     = "ntuple";
char PNGName[300]="AnalizeDatiJPsiControlChannel_test2018";
char NameList[300]="xcut-namelist.lis";

std::stringstream sss;
double XMinSign = 5.0;
double XMaxSign = 5.6;
//
double XMinSBL  = XMinSign;
double XMaxSBL  = 5.15;
double XMinSBR  = 5.4;
double XMaxSBR  = XMaxSign;
//
double XStepSign = 0.0025*2;
float xMassHBin = (XMaxSign -XMinSign)/XStepSign;
float xMassHBin2   =  xMassHBin ; // plot only!
//
bool MC = false;
bool CorreTag = false;
bool WrongTag = false;
bool Cutkst2 =  false;
bool Cutmmpi2=  false ;
bool CutCosK =  false ;
//
  double x0Cut=-0.4;
  double y0Cut= 0.3;
  double x1Cut= 0.6;
  double y1Cut=-0.10;
//  
  double x0Cut2=-0.4;
  double y0Cut2= 0.3;
  double x1Cut2= 0.6;
  double y1Cut2=-0.10;
//  
  double x_0Cut=3;
  double y_0Cut=3.8;
  double x_1Cut=3.6;
  double y_1Cut=4.8;
//  
  double CutX1=3.2;
  double CutX2=3.6;
  double CutY1=4.7;
  double CutY2=4.9;
//  

//
void replaceAll(std::string& str, const std::string& from, const std::string& to) ;
void replaceChar(char * txt, const char * txt1, const char * txt2) ;
bool Background( double m12, double m13, double DmMass, double m1 , double m2, double m3);
double  KineLimit(  double *x, double *par);
double  FunCutkst2( double *x, double *par);
double  FunCutmmpi2( double *x, double *par);
double  BackgroundPlot( double m12, double m13, double DmMass, double m1 , double m2, double m3);
std::map<std::string, std::string>  ReadNamelist(int argc, char** argv);
//
int main (int argc, char** argv) {
  if (argc<=1 ){
   std::cout<<"Error: please set the ERA [2016 2017 2018]\n"<<std::endl;
   exit(1);
  }else{ 
     RunEra = (char *) malloc(strlen(argv[1])+1);
     strcpy(RunEra,argv[1]);
     if ((strcmp(argv[1],"2016") == 0)){
       std::cout<<Form("==>Setting Era = %s\n",RunEra)<<std::endl;
       replaceChar(PNGName,"2018",RunEra);
       replaceChar(OutFileNameInputHisto  ,"2018",RunEra);
       replaceChar(OutFileNameInputHistoMC,"2018",RunEra);
       replaceChar(OutFileNameInputHistoMCJPsi,"2018",RunEra);
     }else if ((strcmp(argv[1],"2017") == 0)){
       std::cout<<Form("==>Setting Era = %s\n",RunEra)<<std::endl;
       replaceChar(PNGName,"2018",RunEra);
       replaceChar(OutFileNameInputHisto  ,"2018",RunEra);
       replaceChar(OutFileNameInputHistoMC,"2018",RunEra);
       replaceChar(OutFileNameInputHistoMCJPsi,"2018",RunEra);
     }else if ((strcmp(argv[1],"2018") == 0)){
       std::cout<<Form("==>Setting Era = %s\n",RunEra)<<std::endl;
     }else{
      std::cout<<Form("Error: Era %s not recognize => set the ERA [2016 2017 2018]\n",RunEra)<<std::endl;
      exit(1);
     }
     if (argc>2 && ((strcmp(argv[2],"m") == 0)||(strcmp(argv[2],"M") == 0)) ){
      cout<<"Read MC  !!!" <<endl;
      sprintf(OutFileNameInputHisto,"%s",OutFileNameInputHistoMC);
      replaceChar(PNGName,"Dati","MC");
      MC = true;
      if (argc>3 && ((strcmp(argv[3],"tag") == 0)||(strcmp(argv[3],"TAG") == 0)) ){
     	replaceChar(PNGName,"MC","MCCorreTag");
     	CorreTag=true;
        std::cout<<Form("==>Setting MC Correctly Tagged RunEra=\n",RunEra)<<std::endl;
      } else if (argc>3 && ((strcmp(argv[3],"mis") == 0)||(strcmp(argv[3],"MIS") == 0)) ){
     	replaceChar(PNGName,"MC","MCWrongTag");
     	WrongTag=true;
        std::cout<<Form("==>Setting MC Wrongly   Tagged RunEra=\n",RunEra)<<std::endl;
      }else{
        std::cout<<Form("==>Setting MC All Tags RunEra=\n",RunEra)<<std::endl;
      }
 
     }else if (argc>2 && ((strcmp(argv[2],"m2") == 0)||(strcmp(argv[2],"M2") == 0)) ){
 
      cout<<"Read MC  Psi !!!" <<endl;
//   	if (argc>2 && ((strcmp(argv[2],"tag") == 0)||(strcmp(argv[2],"TAG") == 0)) ) CorreTag=true;
//   	if (argc>2 && ((strcmp(argv[2],"mis") == 0)||(strcmp(argv[2],"MIS") == 0)) ) WrongTag=true;
      sprintf(OutFileNameInputHisto,"%s",OutFileNameInputHistoMCJPsi);
      replaceChar(PNGName,"Dati","MCJPsi");
      XStepSign=XStepSign*5;
      xMassHBin = xMassHBin/5;
      xMassHBin2 = xMassHBin;
      MC = true;
 
     }else{
      cout<<"Read Data!!!" <<endl;
      MC = false;
     }
  }
  char*argn[]={NameList};
  
  std::map<std::string, std::string> mappa = ReadNamelist(1,argn );
//
  x0Cut  =  atof (mappa["x0Cut"].c_str() ) ;
  y0Cut  =  atof (mappa["y0Cut"].c_str() ) ;
  x1Cut  =  atof (mappa["x1Cut"].c_str() ) ;
  y1Cut  =  atof (mappa["y1Cut"].c_str() ) ;

  x0Cut2 =  atof (mappa["x0Cut2"].c_str() ) ;
  y0Cut2 =  atof (mappa["y0Cut2"].c_str() ) ;
  x1Cut2 =  atof (mappa["x1Cut2"].c_str() ) ;
  y1Cut2 =  atof (mappa["y1Cut2"].c_str() ) ;

  x_0Cut =  atof (mappa["x_0Cut"].c_str() ) ;
  y_0Cut =  atof (mappa["y_0Cut"].c_str() ) ;
  x_1Cut =  atof (mappa["x_1Cut"].c_str() ) ;
  y_1Cut =  atof (mappa["y_1Cut"].c_str() ) ;

  CutX1  =  atof (mappa["CutX1"].c_str() ) ;
  CutX2  =  atof (mappa["CutX2"].c_str() ) ;
  CutY1  =  atof (mappa["CutY1"].c_str() ) ;
  CutY2  =  atof (mappa["CutY2"].c_str() ) ;

  Cutkst2 =  ( atoi(mappa["Cutkst2"].c_str())==1) ;
  Cutmmpi2=  ( atoi(mappa["Cutmmpi2"].c_str())==1) ;
  CutCosK=  ( atoi(mappa["CutCosK"].c_str())==1) ;
  
  map<string,string>::iterator  it= mappa.find("XMaxSBL");
  if(it != mappa.end()) {
   XMaxSBL  =  atof (mappa["XMaxSBL"].c_str() ) ;
   std::cout<<Form(" Setting XMaxSBL = %f \n",XMaxSBL)<<std::endl;
  }
  map<string,string>::iterator  ip= mappa.find("XMinSBR");
  if(ip != mappa.end()) {
   XMinSBR  =  atof (mappa["XMinSBR"].c_str() ) ;
   std::cout<<Form(" Setting XMinSBR = %f \n",XMinSBR)<<std::endl;
  }
  
  if(Cutkst2 ) {
   std::cout<<" Setting Cut kst2 \n"<<std::endl;
   sprintf(PNGName,"%s-Cutkst2",PNGName);
  } 
  if(Cutmmpi2) {
   std::cout<<" Setting Cut mmpi2\n"<<std::endl;
   sprintf(PNGName,"%s-Cutmmpi2",PNGName);
  } 
  if(CutCosK) {
   std::cout<<" Setting Cut for plot CutCosK\n"<<std::endl;
   sprintf(PNGName,"%s-ScatterPlotCutCosK",PNGName);
  } 
  
  
  
  
  TStopwatch TimeWatch;
  TimeWatch.Start();

  AnalizeDatiJPsi_test(); 
//  app.Run() ;
  cout<<"esco..." <<endl;
  TimeWatch.Stop();
  TimeWatch.Print();
  
  return 0 ;
}
void  AnalizeDatiJPsi_test(){


  double cos_theta_l	  ;
  double cos_theta_k	  ;
  double phi_kst_mumu	  ;
  double tagged_mass	  ;
  double mumuMass	  ;
  double mumuMassE	  ;
  double mumuPt   	  ;
  double mumuPhi   	  ;
  double mumuEta   	  ;
  int genBro3PdgId	  ;
  int genBro2PdgId	  ;
  int genBro1PdgId	  ;
  int genJpsiAncestorPdgId;
  double muonMass_ = 0.1056583745;
  double piMass = 0.13957039;
  double kMass = 0.493677;
  double BpMass = 5.2791;
  double B0Mass = 5.27962;
  double PMass = 0.939;
  double KstarMass = 0.896;
//  double KstarMass = 0.892;
  double gen_jpsimass_12 ;
  double mumTrkp ;
  double mupTrkm ;
  double mmk1    ;
  double mmk2    ;
  double kstarmass ;
  double kstBarMass ;
  double kstMass ;
  double MassPsi2S=3.696;
  double MassJPsi =3.0969;
  double mmpip;
  double mmpim;
  double mmpipkm;
  double mmpimkp;
  double mmpiKaon;
  double mmpiKaon2;
  double mmpi1;
  double mmpi2;
  double mmpipi;
  double tagB0;
  double mumPt;
  double mumPhi;
  double mumEta;
  double mupPt;
  double mupPhi;
  double mupEta;
  double kstTrkmPt;
  double kstTrkmPhi;
  double kstTrkmEta;
  double kstTrkpPt;
  double kstTrkpPhi;
  double kstTrkpEta;
  bool   passB0Psi_jpsi ;
  bool   passB0Psi_psip ;
  double  genSignal;
  
  TCanvas* cStudies = new TCanvas("cStudies","Cutmmpi2&Cutkst2",200,10,900,780);
  cStudies->Divide(4,2);
  TCanvas* cDalitz      = new TCanvas("cDalitz","#mu#mu#pi^2%#mu#muK^2",200,10,900,780);
  TCanvas* cDalitzSBL   = new TCanvas("cDalitzSBL","#mu#mu#pi^2%#mu#muK^2",200,10,900,780);
  TCanvas* cDaliCut     = new TCanvas("cDaliCut","#mu#mu#pi^2%#mu#muK^2 Cutted",200,10,900,780);

  if (!TFile::Open(OutFileNameInputHisto,"READ"))
  {
    cout<<"File:"<<OutFileNameInputHisto<<" not found!!! create..."<<endl;
    exit(1);
  }else{
   OutFileInputHisto = TFile::Open(OutFileNameInputHisto,"READ");
   cout<<"File:"<<OutFileNameInputHisto<<" FOUND !!!"<<endl;
  } 
  TTree *RecoB0TreeOut     = (TTree*)OutFileInputHisto->Get(OutputRecoB0TreeName);
   if(!RecoB0TreeOut ){
     cout<<"TTree Reco Data: "<< OutputRecoB0TreeName <<" not found!!! Suggestion: remove this file e try again..."<<endl;
     exit(1);
   }else{
     cout<<"TTree Reco Data: "<< OutputRecoB0TreeName <<" OK FOUND!!!"<<endl;
   }  
     TH1D* Hxpsi2s_Sign	     = new TH1D( "Hxpsi2s_Sign"	    , "psi2 mass"	                             ,120, 3.5, 3.8);
     TH1D* Hxpsi2s_Sign_c    = new TH1D( "Hxpsi2s_Sign_c"   , "psi2 mass"	                             ,120, 3.5, 3.8);
     TH1D* HxLambdab_Sign	     = new TH1D( "HxLambdab_Sign"	    , "Lambdab Mass"	                             ,50, 5.45, 5.9);
     TH1D* HxLambdab_Sign_c	     = new TH1D( "HxLambdab_Sign_c"	    , "Lambdab Mass"	                             ,50, 5.45, 5.9);
     TH1D* HxLambdab_SBL	     = new TH1D( "HxLambdab_SBL"	    , "Lambdab Mass"	                             ,50, 5.45, 5.9);
     TH1D* HxLambdab_SBR	     = new TH1D( "HxLambdab_SBR"	    , "Lambdab Mass"	                             ,50, 5.45, 5.9);
     TH1D* HxLambdab_Signc	     = new TH1D( "HxLambdab_Signc"	    , "Lambdab Mass c"	                             ,50, 5.45, 5.9);
     TH1D* HxBsMass	     = new TH1D( "HxBsMass"	    , "Bs Mass"	                             ,xMassHBin2, XMinSign, XMaxSign);
     TH1D* HxBsMass_Sign     = new TH1D( "HxBsMass_Sign"    , "Bs Mass Sign"	                     ,xMassHBin2, XMinSign, XMaxSign);
     TH1D* HxBsMass_SignC    = new TH1D( "HxBsMass_SignC"   , "Bs Mass Sign C"	                     ,xMassHBin2, XMinSign, XMaxSign);
     TH1D* HxBsMass_SBL	     = new TH1D( "HxBsMass_SBL"	    , "Bs Mass SBL"	                     ,xMassHBin2, XMinSign, XMaxSign);
     TH1D* HxBsMass_SBR	     = new TH1D( "HxBsMass_SBR"	    , "Bs Mass SBR"	                     ,xMassHBin2, XMinSign, XMaxSign);
     TH1D* HxMass	     = new TH1D( "HxMass"           , "B^{0} Mass"			     ,xMassHBin2, XMinSign, XMaxSign);
     TH1D* HxMass_Sign	     = new TH1D( "HxMass_Sign"      , "B^{0} Mass Sign"			     ,xMassHBin2, XMinSign, XMaxSign);
     TH1D* HxMass_SBL 	     = new TH1D( "HxMass_SBL"       , "B^{0} Mass SBL"			     ,xMassHBin2, XMinSign, XMaxSign);
     TH1D* HxMass_SBR 	     = new TH1D( "HxMass_SBR"       , "B^{0} Mass SBR"			     ,xMassHBin2, XMinSign, XMaxSign);
//   TH1D* HxMass_Cutted	     = new TH1D( "HxMass_Cutted"    , "B^{0} Mass rejected cut1"	     ,xMassHBin2, XMinSign, XMaxSign);
//   TH1D* HxMass_Cutted2	     = new TH1D( "HxMass_Cutted2"   , "B^{0} Mass rejected cut2"	     ,xMassHBin2, XMinSign, XMaxSign);
//   TH1D* HxMass_Cutted3	     = new TH1D( "HxMass_Cutted3"   , "B^{0} Mass rejected cut3"	     ,xMassHBin2, XMinSign, XMaxSign);
//   TH1D* HxMass_CutApp	     = new TH1D( "HxMass_CutApp"    , "B^{0} Mass cut1 applied" 	     ,xMassHBin2, XMinSign, XMaxSign);
//   TH1D* HxMass_CutApp2	     = new TH1D( "HxMass_CutApp2"   , "B^{0} Mass cut2 applied" 	     ,xMassHBin2, XMinSign, XMaxSign);
//   TH1D* HxMass_CutApp3	     = new TH1D( "HxMass_CutApp3"   , "B^{0} Mass cut2 applied" 	     ,xMassHBin2, XMinSign, XMaxSign);
//   TH1D* HxMass_Sign 	     = new TH1D( "HxMass_Sign"      , "B^{0} Mass (signal)"		     ,xMassHBin2, XMinSign, XMaxSign);
//   TH1D* HxMass_SBL_Cutted    = new TH1D( "HxMass_SBL_Cutted","B^{0} Mass (signal cut1 applied)"      ,xMassHBin2, XMinSign, XMaxSign);
//   TH1D* HxCosL_SBL	     = new TH1D( "HxCosL_SBL"       , "B^{0} CosL SBL"			     ,    80, -1., 1.);
//   TH1D* HxCosL_Sign	     = new TH1D( "HxCosL_Sign"      , "B^{0} CosL (signal)"		     ,    80, -1., 1.);
//   TH1D* HxCosL_Sign_CutApp   = new TH1D( "HxCosL_Sign_CutApp"    , "B^{0} CosL (signal cut1 applied)",	  80, -1., 1.);
//   TH1D* HxCosL_SBL_Cutted    = new TH1D( "HxCosL_SBL_Cutted"     , "B^{0} CosL (SBL rejected cut1)"  ,	  80, -1., 1.);
//   TH1D* HxCosL_SBL_Cutted2   = new TH1D( "HxCosL_SBL_Cutted2"    , "B^{0} CosL (SBL rejected cut2)"  ,	  80, -1., 1.);
//   TH1D* HxCosL_SBL_Cutted3   = new TH1D( "HxCosL_SBL_Cutted3"    , "B^{0} CosL (SBL rejected cut3)"  ,	  80, -1., 1.);
//   TH1D* HxCosL_SBL_CutApp    = new TH1D( "HxCosL_SBL_CutApp"     , "B^{0} CosL (SBL cut1 applied)"   ,	  80, -1., 1.);
//   TH1D* HxCosL_SBL_CutApp2   = new TH1D( "HxCosL_SBL_CutApp2"    , "B^{0} CosL (SBL cut2 applied)"   ,	  80, -1., 1.);
//   TH1D* HxCosL_SBL_CutApp3   = new TH1D( "HxCosL_SBL_CutApp3"    , "B^{0} CosL (SBL cut3 applied)"   ,	  80, -1., 1.);
     TH1D* HxCosK	     = new TH1D( "HxCosK"           , "B^{0} CosK"		    ,	 40, -1., 1.);
     TH1D* HxCosK_LK	     = new TH1D( "HxCosK_Lpi"       , "B^{0} CosK leading pi"		    ,	 40, -1., 1.);
     TH1D* HxCosK_Lpi	     = new TH1D( "HxCosK_LK"        , "B^{0} CosK leading K"		    ,	 40, -1., 1.);
     TH1D* HxCosK_SBR	     = new TH1D( "HxCosK_SBR"	    , "B^{0} CosK SBR"		    ,	 40, -1., 1.);
     TH1D* HxCosK_SBL	     = new TH1D( "HxCosK_SBL"       , "B^{0} CosK SBL"		    ,	 40, -1., 1.);
     TH1D* HxCosK_Sign	     = new TH1D( "HxCosK_Sign"      , "B^{0} CosK (signal)"	    ,	 40, -1., 1.);
     TH1D* HxCosK_SBR_Lpi    = new TH1D( "HxCosK_SBR_Lpi"   , "B^{0} CosK SBR leading pi"	   ,	40, -1., 1.);
     TH1D* HxCosK_SBL_Lpi    = new TH1D( "HxCosK_SBL_Lpi"   , "B^{0} CosK SBL leading pi"	   ,	40, -1., 1.);
     TH1D* HxCosK_Sign_Lpi   = new TH1D( "HxCosK_Sign_Lpi"  , "B^{0} CosK (signal) leading pi"	   ,	40, -1., 1.);
     TH1D* HxCosK_SBR_LK     = new TH1D( "HxCosK_SBR_LK"    , "B^{0} CosK SBR leading K"	   ,	40, -1., 1.);
     TH1D* HxCosK_SBL_LK     = new TH1D( "HxCosK_SBL_LK"    , "B^{0} CosK SBL leading K"	   ,	40, -1., 1.);
     TH1D* HxCosK_Sign_LK    = new TH1D( "HxCosK_Sign_LK"   , "B^{0} CosK (signal) leading K"      ,	40, -1., 1.);
//   TH1D* HxCosK_Sign_CutApp   = new TH1D( "HxCosK_Sign_CutApp"    , "B^{0} CosK (sign cut)",		  80, -1., 1.);
//   TH1D* HxCosK_Sign_SBL_Subtr= new TH1D( "HxCosK_Sign_SBL_Subtr" , "B^{0} CosK (sign cut SB subtr)",	  80, -1., 1.);
//   TH1D* HxCosK_Sign_SBL_SubNC= new TH1D( "HxCosK_Sign_SBL_SubNC" , "B^{0} CosK (sign SB subtr no cut)",	  80, -1., 1.);
//   TH1D* HxCosK_SBL_Cutted    = new TH1D( "HxCosK_SBL_Cutted"     , "B^{0} CosK (SBL rejected cut1)",	  80, -1., 1.);
//   TH1D* HxCosK_SBL_Cutted2   = new TH1D( "HxCosK_SBL_Cutted2"    , "B^{0} CosK (SBL rejected cut2)",	  80, -1., 1.);
//   TH1D* HxCosK_SBL_Cutted3   = new TH1D( "HxCosK_SBL_Cutted3"    , "B^{0} CosK (SBL rejected cut3)",	  80, -1., 1.);
//   TH1D* HxCosK_SBL_CutApp    = new TH1D( "HxCosK_SBL_CutApp"     , "B^{0} CosK (SBL cut1 applied)" ,	  80, -1., 1.);
//   TH1D* HxCosK_SBL_CutApp2   = new TH1D( "HxCosK_SBL_CutApp2"    , "B^{0} CosK (SBL cut2 applied)" ,	  80, -1., 1.);
//   TH1D* HxCosK_SBL_CutApp3   = new TH1D( "HxCosK_SBL_CutApp3"    , "B^{0} CosK (SBL cut3 applied)" ,      80, -1., 1.);
//   TH1D* HxCosK_SBR           = new TH1D( "HxCosK_SBR"            , "B^{0} CosK SBR"                   , 80, -1., 1.);
//   TH1D* HxCosK_SBR_Scaled    = new TH1D( "HxCosK_SBR_Scaled"     , "B^{0} CosK (SBR scaled)"          , 80, -1., 1.);
//   TH1D* HxCosK_SBR_Scaled_C  = new TH1D( "HxCosK_SBR_Scaled_C"   , "B^{0} CosK" 		     , 80, -1., 1.);
//   TH1D* HxPhi_SBL	     = new TH1D( "HxPhi_SBL_"	   	 , "B^{0} Phi SBL"    		      , 90, -TMath::Pi(), TMath::Pi());
//   TH1D* HxPhi_Sign	     = new TH1D( "HxPhi_Sign"	   	 , "B^{0} Phi (signal)" 	      , 90, -TMath::Pi(), TMath::Pi());
//   TH1D* HxPhi_Sign_CutApp    = new TH1D( "HxPhi_Sign_CutApp"     , "B^{0} Phi (signal cut1)"	      , 90, -TMath::Pi(), TMath::Pi());
//   TH1D* HxPhi_SBL_Cutted     = new TH1D( "HxPhi_SBL_Cutted"      , "B^{0} Phi (SBL rejected cut1)"    , 90, -TMath::Pi(), TMath::Pi());
//   TH1D* HxPhi_SBL_Cutted2    = new TH1D( "HxPhi_SBL_Cutted2"     , "B^{0} Phi (SBL rejected cut2)"    , 90, -TMath::Pi(), TMath::Pi());
//   TH1D* HxPhi_SBL_Cutted3    = new TH1D( "HxPhi_SBL_Cutted3"     , "B^{0} Phi (SBL rejected cut2)"    , 90, -TMath::Pi(), TMath::Pi());
//   TH1D* HxPhi_SBL_CutApp     = new TH1D( "HxPhi_SBL_CutApp"      , "B^{0} Phi (SBL cut1 applied)"     , 90, -TMath::Pi(), TMath::Pi());
//   TH1D* HxPhi_SBL_CutApp2    = new TH1D( "HxPhi_SBL_CutApp2"     , "B^{0} Phi (SBL cut2 applied)"     , 90, -TMath::Pi(), TMath::Pi());
//   TH1D* HxPhi_SBL_CutApp3    = new TH1D( "HxPhi_SBL_CutApp3"     , "B^{0} Phi (SBL cut3 applied)"     , 90, -TMath::Pi(), TMath::Pi());
//   TH1D* HxPhi_SBR	     = new TH1D( "HxPhi_SBR"	  , "B^{0} Phi SBR" 	    ,		90, -TMath::Pi(), TMath::Pi());
//   TH2D* Hxmmk2mmk1	     = new TH2D( "Hxmmk2mmk1"  , "B^{0} mmk2%mmk1"	,	    150,  3.,6., 150,  3.,6.);
//   TH2D* Hxmmk2mmk1_SBR	     = new TH2D( "Hxmmk2mmk1_SBR" , "B^{0} mmk2%mmk1 SBR"	,	    150,  3.,6., 150,  3.,6.);
//   TH2D* Hxmmk2mmk1S	     = new TH2D( "Hxmmk2mmk1S" , "B^{0} mmk2%mmk1"	,	    150,  3.,6., 150,  3.,6.);
// //
//   TH2D* Hxmmpi2mmka1_CutApp      = new TH2D( "Hxmmpi2mmka1_CutApp", "B^{0} mmpi2%mmka1 (cut1 applied)"  	    ,    	150,  3.,6., 150,  3.5,6.);
//   TH2D* Hxmmpi2mmka1_Cutted	 = new TH2D( "Hxmmpi2mmka1_Cutted", "B^{0} mmpi2%mmka1  rejected"   		    ,		150,  3.,6., 150,  3.5,6.);
//   TH2D* Hxmmpi2mmka1_Sign        = new TH2D( "Hxmmpi2mmka1_Sign"  , "B^{0} mmpi2%mmka1 (Sign)"			    ,		150,  3.,6., 150,  3.5,6.);
//   TH2D* Hxmmpi2mmka1_Sign_CutApp = new TH2D( "Hxmmpi2mmka1_Sign_CutApp", "B^{0} mmpi2%mmka1 (Sign cut1 applied)"    ,		150,  3.,6., 150,  3.5,6.);
//   TH2D* Hxmmpi2mmka1_Sign_Cutted = new TH2D( "Hxmmpi2mmka1_Sign_Cutted", "B^{0} mmpi2%mmka1 (Sign rejected cut1)"   ,		150,  3.,6., 150,  3.5,6.);
//   TH2D* Hxmmpi2mmka1_SBL 	 = new TH2D( "Hxmmpi2mmka1_SBL"        , "B^{0} mmpi2%mmka1 (SBL)"		    ,		150,  3.,6., 150,  3.5,6.);
//   TH2D* Hxmmpi2mmka1_SBL_CutApp	 = new TH2D( "Hxmmpi2mmka1_SBL_CutApp" , "B^{0} mmpi2%mmka1 (SBL cut1 applied)"     ,		150,  3.,6., 150,  3.5,6.);
//   TH2D* Hxmmpi2mmka1_SBL_Cutted	 = new TH2D( "Hxmmpi2mmka1_SBL_Cutted" , "B^{0} mmpi2%mmka1 (SBL rejected cut1)"    ,		150,  3.,6., 150,  3.5,6.);
// //
//   TH2D* HxdKstdmmpi2		 = new TH2D( "HxdKstdmmpi2", "d-kst2%d-d-mmpi2" 	  ,	       150,  -0.4,1.3, 150,  -1.3,0.8);
//   TH2D* HxdKstdmmpi2_CutApp	 = new TH2D( "HxdKstdmmpi2_CutApp", "d-kst2%d-d-mmpi2 Cut",	       150,  -0.4,1.3, 150,  -1.3,0.8);
// //
//   TH2D* HxdKst2dBp		 = new TH2D( "HxdKst2dBp"	      , "B^{0} (kst2-K0*Mass)%(B+ Mass-mmpiKaon)"	           ,	100,  -0.4,0.6, 100,  -0.4,0.6);
//   TH2D* HxdKst2dBp_CutApp	 = new TH2D( "HxdKst2dBp_CutApp"      , "B^{0} (kst2-K0*Mass)%(B+ Mass-mmpiKaon) (Cut1 applied)" ,	100,  -0.4,0.6, 100,  -0.4,0.6);
//   TH2D* HxdKst2dBp_Cutted	 = new TH2D( "HxdKst2dBp_Cutted"      , "B^{0} (kst2-K0*Mass)%(B+ Mass-mmpiKaon) (Cut1 rejected)" ,	100,  -0.4,0.6, 100,  -0.4,0.6);
//   TH2D* HxdKst2dBp_SBL		 = new TH2D( "HxdKst2dBp_SBL"	      , "B^{0} (kst2-K0*Mass)%(B+ Mass-mmpiKaon) (SBL)"		   ,	100,  -0.4,0.6, 100,  -0.4,0.6);
//   TH2D* HxdKst2dBp_SBL_Cutted	 = new TH2D( "HxdKst2dBp_SBL_Cutted"  , "B^{0} (kst2-K0*Mass)%(B+ Mass-mmpiKaon) (SBL Cut1 rejected)",	100,  -0.4,0.6, 100,  -0.4,0.6);
//   TH2D* HxdKst2dBp_SBL_CutApp	 = new TH2D( "HxdKst2dBp_SBL_CutApp"  , "B^{0} (kst2-K0*Mass)%(B+ Mass-mmpiKaon) (SBL Cut1 applied)" ,	100,  -0.4,0.6, 100,  -0.4,0.6);
//   TH2D* HxdKst2dBp_Sign 	 = new TH2D( "HxdKst2dBp_Sign"        , "B^{0} (kst2-K0*Mass)%(B+ Mass-mmpiKaon) (Sign)"  	   ,	100,  -0.4,0.6, 100,  -0.4,0.6);
//   TH2D* HxdKst2dBp_Sign_Cutted	 = new TH2D( "HxdKst2dBp_Sign_Cutted" , "B^{0} (kst2-K0*Mass)%(B+ Mass-mmpiKaon) (Sign Cut1 rejected)",	100,  -0.4,0.6, 100,  -0.4,0.6);
//   TH2D* HxdKst2dBp_Sign_CutApp	 = new TH2D( "HxdKst2dBp_Sign_CutApp" , "B^{0} (kst2-K0*Mass)%(B+ Mass-mmpiKaon) (Sign Cut1 applied)",	100,  -0.4,0.6, 100,  -0.4,0.6);
// //
//   TH2D* HxdKst1dBp	       = new TH2D( "HxdKst1dBp"              , "B^{0} (kst1-K0*Mass)%(B+ Mass-mmpiKaon2)"             ,	       100,  -0.4,0.6, 100,  -0.4,0.6);
//   TH2D* HxdKst1dBp_CutApp      = new TH2D( "HxdKst1dBp_CutApp"       , "B^{0} (kst1-K0*Mass)%(B+ Mass-mmpiKaon2) Cut1 applied",        100,  -0.4,0.6, 100,  -0.4,0.6);
//   TH2D* HxdKst1dBp_Cutted      = new TH2D( "HxdKst1dBp_Cutted"       , "B^{0} (kst1-K0*Mass)%(B+ Mass-mmpiKaon2) Cut1 applied",        100,  -0.4,0.6, 100,  -0.4,0.6);
//   TH2D* HxdKst1dBp_SBL	       = new TH2D( "HxdKst1dBp_SBL"          , "B^{0} (kst1-K0*Mass)%(B+ Mass-mmpiKaon2) SBL"         ,	       100,  -0.4,0.6, 100,  -0.4,0.6);
//   TH2D* HxdKst1dBp_SBL_CutApp  = new TH2D( "HxdKst1dBp_SBL_CutApp"   , "B^{0} (kst1-K0*Mass)%(B+ Mass-mmpiKaon2) SBL Cut1 applied",    100,  -0.4,0.6, 100,  -0.4,0.6);
//   TH2D* HxdKst1dBp_SBL_Cutted  = new TH2D( "HxdKst1dBp_SBL_Cutted"   , "B^{0} (kst1-K0*Mass)%(B+ Mass-mmpiKaon2) SBL Cut1 applied",    100,  -0.4,0.6, 100,  -0.4,0.6);
//   TH2D* HxdKst1dBp_Sign	       = new TH2D( "HxdKst1dBp_Sign"         , "B^{0} (kst1-K0*Mass)%(B+ Mass-mmpiKaon2) Sign"            ,    100,  -0.4,0.6, 100,  -0.4,0.6);
//   TH2D* HxdKst1dBp_Sign_CutApp = new TH2D( "HxdKst1dBp_Sign_CutApp"  , "B^{0} (kst1-K0*Mass)%(B+ Mass-mmpiKaon2) Sign Cut1 applied",   100,  -0.4,0.6, 100,  -0.4,0.6);
//   TH2D* HxdKst1dBp_Sign_Cutted = new TH2D( "HxdKst1dBp_Sign_Cutted"  , "B^{0} (kst1-K0*Mass)%(B+ Mass-mmpiKaon2) Sign Cut1 applied",   100,  -0.4,0.6, 100,  -0.4,0.6);
// //
//   TH1D* HxMuTrkD	     = new TH1D( "HxMuTrkD"   , "B^{0} mmk1-mmk2"      ,	    200,  -2.,2.);
//   TH1D* Hxmmk1_SBL	     = new TH1D( "Hxmmk1_SBL" , "B^{0} mmk1 SBL" 	  ,	       100,  5.,6.);
//   TH1D* Hxmmk2_SBL	     = new TH1D( "Hxmmk2_SBL" , "B^{0} mmk2 SBL" 	  ,	       100,  5.,6.);
//   TH1D* Hxmmk1		     = new TH1D( "Hxmmk1"     , "B^{0} mmk1"	       ,	    100,  5.,6.);
//   TH1D* Hxmmk2		     = new TH1D( "Hxmmk2"     , "B^{0} mmk2"	       ,	    100,  5.,6.);
//   TH1D* HxmumuMass	     = new TH1D( "HxmumuMass"  , "mumuMass"	       ,	    100,  2.,4.);
//   TH1D* HxmumuMass_SBL	     = new TH1D( "HxmumuMass_SBL"  , "mumuMass SBL"	       ,	    100,  2.,4.);
//   TH1D* HxmumuMass_SBL_Cutted= new TH1D( "HxmumuMass_SBL_Cutted" , "mumuMass"	       ,	    100,  2.,4.);
//   TH1D* HxmumuMass_SBL_CutApp= new TH1D( "HxmumuMass_SBL_CutApp" , "mumuMass"	       ,	    100,  2.,4.);
//   TH1D* Hxkstarmass_SBL	     = new TH1D( "Hxkstarmass_SBL"   	 , "K^{*0} SBL"        ,	    100,  0.,2.);
//   TH1D* Hxkstarmass_SBL_Cutted = new TH1D( "Hxkstarmass_SBL_Cutted"  , "K^{*0} SBL rejected cut1 ",	    100,  0.,2.);
//   TH1D* HxkstarmassC	     = new TH1D( "HxkstarmassC"  , "K^{*0}"	       ,	    100,  0.,2.);
//   TH1D* HxMin		     = new TH1D( "HxMin"      , "xmin"  	       ,	    100,  -0.5,0.5);
//   TH1D* Hxmmpipi	     = new TH1D( "Hxmmpipi"   , "min mumupipi"  	 ,	    100,  3.,4.);
//   TH1D* HxMin2		     = new TH1D( "HxMin2"     , "min ps2s-mumupi"      ,	    100,  0.,1.);
//   TH1D* HxmmpiKaon	     = new TH1D( "HxmmpiKaon" , "min mumupiKaon"       ,	    100,  4.,7.);
//   TH1D* Hxkst2		     = new TH1D( "Hxkst2"      , "K^{*0} swap"         ,	    100,  0.6,2.);
//   TH1D* Hxkst2kstarmass        = new TH1D( "Hxkst2kstarmass", "K^{*0} comb 2- K^{*0} PDG",     100,  -0.4,0.6);
//      TH2D* Hxmmpikstarmass_Sign       = new TH2D( "Hxmmpikstarmass_Sign"      , "mmpi vs kstar mass",	          100, 14,26,100,0.3,1.2);
//      TH2D* Hxmmpikstarmass_Sign_c     = new TH2D( "Hxmmpikstarmass_Sign_c"    , "mmpi vs kstar mass",	          100, 14,26,100,0.3,1.2);
//      TH2D* Hxmmpikstarmass_Sign_LK    = new TH2D( "Hxmmpikstarmass_Sign_LK"   , "mmpi vs kstar mass lead K",          100, 14,26,100,0.3,2.2);
//      TH2D* Hxmmpikstarmass_Sign_LK_c  = new TH2D( "Hxmmpikstarmass_Sign_LK_c" , "mmpi vs kstar mass lead K comp",     100, 14,26,100,0.3,2.2);
//      TH2D* Hxmmpikstarmass_Sign_Lpi   = new TH2D( "Hxmmpikstarmass_Sign_Lpi"  , "mmpi vs kstar mass lead pi",         100, 14,26,100,0.3,1.2);
//      TH2D* Hxmmpikstarmass_Sign_Lpi_c = new TH2D( "Hxmmpikstarmass_Sign_Lpi_c", "mmpi vs kstar mass lead pi comp",    100, 14,26,100,0.3,1.2);

     TH2D* Hxmmpikstarmass_Sign_LK    = new TH2D( "Hxmmpikstarmass_Sign_LK"   , "mmpi vs kstar mass lead K",          100, 10,24,100,0.3,2.2);
     TH2D* Hxmmpikstarmass_Sign_LK_c  = new TH2D( "Hxmmpikstarmass_Sign_LK_c" , "mmpi vs kstar mass lead K comp",     100, 10,24,100,0.3,2.2);
     TH2D* Hxmmpikstarmass_Sign_Lpi   = new TH2D( "Hxmmpikstarmass_Sign_Lpi"  , "mmpi vs kstar mass lead pi",         100, 10,24,100,0.3,1.2);
     TH2D* Hxmmpikstarmass_Sign_Lpi_c = new TH2D( "Hxmmpikstarmass_Sign_Lpi_c", "mmpi vs kstar mass lead pi comp",    100, 10,24,100,0.3,1.2);
//   TH1D* HxmmpiKaonBpMass     = new TH1D( "HxmmpiKaonBpMass", "B^{+} PDG - mumupiK comb",   100,  -0.4,0.6);
//   double xbin1 = -0.4;
//   double xbin2 = 0.4;
//   double ybin1 = -1;
//   double ybin2 = 1;
//   
//   TH2D* Hxmmpi2BptMass_CutApp2 = new TH2D( "Hxmmpi2BptMass_CutApp2", "B^{0} mmpi2%BpMass"   ,	      150,  xbin1,xbin2, 150,  ybin1,ybin2);
//   TH2D* Hxmmpi2BptMass_SBL_CutApp2= new TH2D( "Hxmmpi2BptMass_SBL_CutApp2", "B^{0} mmpi2%BpMass Left"   ,   150,  xbin1,xbin2, 150,  ybin1,ybin2);
//   TH2D* Hxmmpi2BptMass_SBR_CutApp2= new TH2D( "Hxmmpi2BptMass_SBR_CutApp2", "B^{0} mmpi2%BpMass Right"   ,  150,  xbin1,xbin2, 150,  ybin1,ybin2);
//   TH2D* Hxmmpi2mmpi1   = new TH2D( "Hxmmpi2mmpi1", "B^{0} mmpi2%mmpi1"   ,	      150,  3.,6., 150,  3.,6.);

  TH2D* HxMassCosK      = new TH2D( "HxMassCosK"     , "B^{0} Mass%Cosk"       ,	      60,  5.,5.6, 100,  -1.,1.);
  TH2D* HxMassCosK_Sign = new TH2D( "HxMassCosK_Sign", "B^{0} Mass%Cosk Sign"  ,	      60,  5.,5.6, 100,  -1.,1.);
  TH2D* HxMassCosK_SBL  = new TH2D( "HxMassCosK_SBL" , "B^{0} Mass%Cosk SBL"   ,	      60,  5.,5.6, 100,  -1.,1.);
  TH2D* HxMassCosK_SBR  = new TH2D( "HxMassCosK_SBR" , "B^{0} Mass%Cosk SBR"   ,	      60,  5.,5.6, 100,  -1.,1.);

  int nentries = (int)RecoB0TreeOut->GetEntries(); 
  cout<<"nentries: "<< nentries<<endl;
  RecoB0TreeOut->SetBranchAddress("tagB0"                ,&tagB0);
  RecoB0TreeOut->SetBranchAddress("passB0Psi_jpsi" 	 ,&passB0Psi_jpsi );
  RecoB0TreeOut->SetBranchAddress("passB0Psi_psip" 	 ,&passB0Psi_psip );
  RecoB0TreeOut->SetBranchAddress("cos_theta_l" 	 ,&cos_theta_l );
  RecoB0TreeOut->SetBranchAddress("cos_theta_k" 	 ,&cos_theta_k );
  RecoB0TreeOut->SetBranchAddress("phi_kst_mumu"	 ,&phi_kst_mumu);
  RecoB0TreeOut->SetBranchAddress("tagged_mass" 	 ,&tagged_mass );
  RecoB0TreeOut->SetBranchAddress("mumuMass"		 ,&mumuMass );
  RecoB0TreeOut->SetBranchAddress("mumuPt"		 ,&mumuPt );
  RecoB0TreeOut->SetBranchAddress("mumuPhi"		 ,&mumuPhi );
  RecoB0TreeOut->SetBranchAddress("mumuEta"		 ,&mumuEta );
  RecoB0TreeOut->SetBranchAddress("tagged_mass" 	 ,&tagged_mass );
  RecoB0TreeOut->SetBranchAddress("mumuMassE"		 ,&mumuMassE );
  RecoB0TreeOut->SetBranchAddress("mumTrkp"		 ,&mumTrkp );
  RecoB0TreeOut->SetBranchAddress("mupTrkm"		 ,&mupTrkm );
  RecoB0TreeOut->SetBranchAddress("mmk1"		 ,&mmk1 );
  RecoB0TreeOut->SetBranchAddress("mmk2"		 ,&mmk2 );
  RecoB0TreeOut->SetBranchAddress("kstarmass"		 ,&kstarmass );
  RecoB0TreeOut->SetBranchAddress("kstMass"		 ,&kstMass );
  RecoB0TreeOut->SetBranchAddress("kstBarMass"		 ,&kstBarMass );
  RecoB0TreeOut->SetBranchAddress("mumPt"		 ,&mumPt);
  RecoB0TreeOut->SetBranchAddress("mumPhi"		 ,&mumPhi);
  RecoB0TreeOut->SetBranchAddress("mumEta"		 ,&mumEta);
  RecoB0TreeOut->SetBranchAddress("mupPt"		 ,&mupPt);
  RecoB0TreeOut->SetBranchAddress("mupPhi"		 ,&mupPhi);
  RecoB0TreeOut->SetBranchAddress("mupEta"		 ,&mupEta);
  RecoB0TreeOut->SetBranchAddress("kstTrkmPt"		 ,&kstTrkmPt);
  RecoB0TreeOut->SetBranchAddress("kstTrkmPhi"		 ,&kstTrkmPhi);
  RecoB0TreeOut->SetBranchAddress("kstTrkmEta"		 ,&kstTrkmEta);
  RecoB0TreeOut->SetBranchAddress("kstTrkpPt"		 ,&kstTrkpPt);
  RecoB0TreeOut->SetBranchAddress("kstTrkpPhi"		 ,&kstTrkpPhi);
  RecoB0TreeOut->SetBranchAddress("kstTrkpEta"		 ,&kstTrkpEta);
  if(MC){
   RecoB0TreeOut->SetBranchAddress("genSignal"            ,&genSignal);
  }
  double pi1Pt=0;
  double pi2Pt=0;
  double pimPt=0;
  double pipPt=0;
  double piKaon =0;
  double pipk =0;
  double pimk =0;
  double mmkp =0;
  double mmkm =0;
  double mmka1 =0;
  double mmka2 =0;
  double kstKp =0;
  double kstKm =0;
  double kst1  =0;
  double kst2  =0;
  double lambdab1=0;
  double lambdab2=0;
  double lambdab=0;
  double nolambdab=0;
  double kstarmassx;
  double kstarmassAll;
  double mmpix;
  double mmkk;
  double masskk;
  double masspipi;
    
  
  double DmMass = B0Mass;
  double m1	= piMass;
  double m2	= kMass;
  double m3	= MassJPsi;
  float m12max=(DmMass-m3)*(DmMass-m3);
  float m13max=(DmMass-m2)*(DmMass-m2);
  float m12min=(m1+m2)*(m1+m2);
  float m13min=(m1+m3)*(m1+m3);
  cout<<Form("m12min=%f m12max=%f m13min=%f m13max=%f",m12min,m12max,m13min,m13max)<<endl;
  TH2F *HxDalitz = new TH2F("HxDalitz","",150,m12min,m12max,150,m13min,m13max); 
  TH2D* Hxmmpikstarmass_SBL	   = new TH2D( "Hxmmpikstarmass_SBL"	   , "mmpi vs kstar mass SBL",         100,m12min,m12max,100,m13min,m13max );
  TH2D* Hxmmpikstarmass_Sign	   = new TH2D( "Hxmmpikstarmass_Sign"	   , "mmpi vs kstar mass",	       100,m12min,m12max,100,m13min,m13max );
  TH2D* Hxmmpikstarmass_Sign_b     = new TH2D( "Hxmmpikstarmass_Sign_b"    , "mmpi vs kstar mass",	       100,m12min,m12max,100,m13min,m13max );
  TH2D* Hxmmpikstarmass_Sign_c     = new TH2D( "Hxmmpikstarmass_Sign_c"    , "mmpi vs kstar mass",	       100,m12min,m12max,100,m13min,m13max );
  TH1D* Hxmmpi_Sign	  = new TH1D( "Hxmmpi_Sign"	 , "JPsi pi mass"	      ,48, m13min,m13max);
  TH1D* Hxmmpiq_SBL	  = new TH1D( "Hxmmpiq_SBL"	 , "(JPsi pi mass)^2"	      ,48, m13min,m13max);
  TH1D* Hxmmpiq_Sign	  = new TH1D( "Hxmmpiq_Sign"	 , "(JPsi pi mass)^2"	      ,48, m13min,m13max);
  TH1D* Hxkstarmass2_SBL   = new TH1D( "Hxkstarmass2_SBL"      , "(K^{*0} mass)^2 SBL"        , 100,  m12min,m12max);
  TH1D* Hxkstarmass2_Sign  = new TH1D( "Hxkstarmass2_Sign"     , "(K^{*0} mass)^2 Sign"       , 100,  m12min,m12max);
  TH1D* Hxkstarmass2_Sign_b= new TH1D( "Hxkstarmass2_Sign_b"   , "(K^{*0} mass)^2 Sign"       , 100,  m12min,m12max);
  TF2  *DKinematic = new TF2("DKinematic",KineLimit,m12min,m12max,m13min,m13max,4);
//  TF2  *DKinematic = new TF2("DKinematic",KineLimit,round(m12min),round(m12max),round(m13min),round(m13max),4);
//  TF2  *DKinematic = new TF2("DKinematic",KineLimit,0.,28.,0.,28.,4);
  TF2  *DKinematic2 = new TF2("DKinematic2",KineLimit,0.,28.,0.,28.,4);
  Double_t Dalipar[4];
  Dalipar[0]=DmMass;
  Dalipar[1]=m1;
  Dalipar[2]=m2;
  Dalipar[3]=m3;
  DKinematic->SetNpy(1000);
  DKinematic->SetNpx(1000);
  DKinematic->SetParameters(Dalipar);
  Dalipar[0]=DmMass;
  Dalipar[1]=m1;
  Dalipar[2]=m2+piMass;
  Dalipar[3]=m3+piMass;
  DKinematic2->SetNpy(1000);
  DKinematic2->SetNpx(1000);
  DKinematic2->SetParameters(Dalipar);
  TF2  *FCutkst2 = new TF2("FCutkst2 ",FunCutkst2,  -0.4,0.6,  -0.4,0.6,4);
  Double_t kst2par[4];
  kst2par[0]=x0Cut;
  kst2par[1]=x1Cut;
  kst2par[2]=y0Cut;
  kst2par[3]=y1Cut;
  FCutkst2->SetNpx(2000);
  FCutkst2->SetNpy(2000);
  FCutkst2->SetParameters(kst2par);
  TF2  *FCutmmpi2 =  new TF2("FCutmmpi2",FunCutmmpi2,3.,6.,3.5,6.,8);
  Double_t mmpi2par[8];
  mmpi2par[0]=CutX1;
  mmpi2par[1]=CutX2;
  mmpi2par[2]=CutY1;
  mmpi2par[3]=CutY2;
  mmpi2par[4]=x_0Cut;
  mmpi2par[5]=x_1Cut;
  mmpi2par[6]=y_0Cut;
  mmpi2par[7]=y_1Cut;
  FCutmmpi2->SetNpx(2000);
  FCutmmpi2->SetNpy(2000);
  FCutmmpi2->SetParameters(mmpi2par);
  
//  TF2  *DKinematic = new TF2("DKinematic",KineLimit, m12min,m12max,m13min,m13max,1);
//
//  double y1Cut2=-0.25;
//nentries=300000;
  for (Int_t i=0;i<nentries;i++) {
         RecoB0TreeOut->GetEntry(i);
//	 
        if( WrongTag&&MC&& (tagB0 == 1 && genSignal == 1) || ( tagB0 == 0 && genSignal == 2)) continue;  //(cut correctly tagged)
        if( CorreTag&&MC&& (tagB0 == 1 && genSignal == 2) || ( tagB0 == 0 && genSignal == 1)) continue;  //(cut wrongly   tagged)
//	 
//	 if(!passB0Psi_psip) continue;
	 if(!passB0Psi_jpsi) continue;
         if (tagged_mass>XMaxSign || tagged_mass<XMinSign) continue;
 	 ROOT::Math::PtEtaPhiMVector MupVec (mupPt,mupEta,mupPhi,muonMass_);
 	 ROOT::Math::PtEtaPhiMVector MumVec (mumPt,mumEta,mumPhi,muonMass_);
 	 ROOT::Math::PtEtaPhiMVector JpsiVec (mumuPt,mumuEta,mumuPhi,mumuMass);
// 	 ROOT::Math::PtEtaPhiMVector JpsiVec (MupVec+MumVec);

//  	 ROOT::Math::PtEtaPhiMVector MKpVec (mupPt,mupEta,mupPhi,muonMass_);
//  	 ROOT::Math::PtEtaPhiMVector MKmVec (mumPt,mumEta,mumPhi,muonMass_);
//  	 ROOT::Math::PtEtaPhiMVector KKVec (MKpVec+MKmVec);

 	 ROOT::Math::PtEtaPhiMVector PipVec (kstTrkpPt,kstTrkpEta,kstTrkpPhi,piMass);
 	 ROOT::Math::PtEtaPhiMVector PimVec (kstTrkmPt,kstTrkmEta,kstTrkmPhi,piMass);
 
 	 ROOT::Math::PtEtaPhiMVector KpVec (kstTrkpPt,kstTrkpEta,kstTrkpPhi,kMass);
 	 ROOT::Math::PtEtaPhiMVector KmVec (kstTrkmPt,kstTrkmEta,kstTrkmPhi,kMass);

 	 ROOT::Math::PtEtaPhiMVector PpVec (kstTrkpPt,kstTrkpEta,kstTrkpPhi,PMass);
 	 ROOT::Math::PtEtaPhiMVector PmVec (kstTrkmPt,kstTrkmEta,kstTrkmPhi,PMass);
	 bool leadK =false;
	 mmpipi = (JpsiVec+PipVec+PimVec).mass();
	 mmpip  = (JpsiVec+PipVec).mass();
	 mmpim  = (JpsiVec+PimVec).mass();
	 mmkp   = (JpsiVec+KpVec).mass();
	 mmkm   = (JpsiVec+KmVec).mass();
	 kstKp  = (KpVec+PimVec).mass();
	 kstKm  = (KmVec+PipVec).mass();
	 pipPt  = (PipVec).pt();
	 pimPt  = (PimVec).pt();
	 mmpipkm = (JpsiVec+PipVec+KmVec).mass();
	 mmpimkp = (JpsiVec+PimVec+KpVec).mass();
	 mmpipi  = (JpsiVec+PimVec+PipVec).mass();
	 mmkk    = (JpsiVec+KmVec+KpVec).mass();
	 masskk = (KpVec+KmVec).mass();
	 masspipi = (PipVec+PimVec).mass();
	 
	 pipk = (PipVec+KmVec).mass();
	 pimk = (PimVec+KpVec).mass();
         if (tagB0>0) { //B0
	  if((KpVec.Pt()>PimVec.Pt())) leadK =true;
//	  lambdab1 =  (JpsiVec+KmVec+PpVec).mass();
//	  lambdab2 =  (JpsiVec+KpVec+PmVec).mass();
	  mmpi1=mmpip;
	  mmpi2=mmpim;
	  kst1 = kstKp;
	  kst2 = kstKm;
	  mmpiKaon=mmpipkm;
	  mmpiKaon2=mmpimkp;
	  pi1Pt=pipPt;
	  pi2Pt=pimPt;
	  piKaon = pipk;
	  mmka1=mmkp;
	  mmka2=mmkm;
	  kstarmassAll=kstMass;
	 }else{ //B0bar
	  if((KmVec.Pt()>PipVec.Pt())) leadK =true;
//	  lambdab1  = (JpsiVec+KpVec+PmVec).mass();
//	  lambdab2 =  (JpsiVec+KmVec+PpVec).mass();
	  
	  mmpi1=mmpim;
	  mmpi2=mmpip;
	  kst1 = kstKm;
	  kst2 = kstKp;
	  mmpiKaon=mmpimkp;
	  mmpiKaon2=mmpipkm;
	  pi1Pt=pimPt;
	  pi2Pt=pipPt;
	  piKaon = pimk;
	  mmka1=mmkm;
	  mmka2=mmkp;
	  kstarmassAll=kstBarMass;
	 }
	 if(PmVec.Pt()>PpVec.Pt()){
	  lambdab  =(JpsiVec+KpVec+PmVec).mass();
	  nolambdab=(JpsiVec+KmVec+PpVec).mass();
	 }else{
	  lambdab  =(JpsiVec+KmVec+PpVec).mass();
	  nolambdab=(JpsiVec+KpVec+PmVec).mass();
	 } 
	 
//	 kstarmass=kstarmassAll;
         mmpix=mmpi2;
	 kstarmassx=kstarmass;
	 
//	 if(kstarmass*kstarmass<0.5) continue;
//	 if(kstarmass>0.8) continue;
//	 if(kstarmass>0.8&&kstarmass<0.98) continue;
//	 if(kstarmass<0.8||kstarmass>1.05) continue;
	 
//	 if(mumuMass>3.685) continue;
	 
//	 if((lambdab<5.63&&lambdab>5.61)) continue;
	 
//	 if ((masskk> 0.96&&masskk<1.08)) continue;
//	 if(((lambdab1<5.64&&lambdab1>5.6) || (lambdab2<5.64&&lambdab2>5.6)) ) continue; 

//         bool KLimit0 = Background( kstarmass*kstarmass, mmpix*mmpix,  tagged_mass , piMass, kMass, MassPsi2S);
//         bool KLimit1 = Background( kst2*kst2, mmpi1*mmpi1,  mmpiKaon , piMass, kMass, MassPsi2S);
         bool KLimit1 = Background( kst2*kst2          , mmpi1*mmpi1,  B0Mass , piMass, kMass, MassPsi2S);
         bool KLimit0 = Background( kstarmassx*kstarmassx, mmpix*mmpix,  B0Mass , piMass, kMass, MassPsi2S);
	 HxMass->Fill(tagged_mass);
//	 if(!KLimit0) continue;
//	 if(!KLimit1) continue;
//	 if(!KLimit0||kstarmass*kstarmass>1.1) continue;
//	 if((!KLimit0&&!KLimit1)) continue;
//	 if( kstarmass*kstarmass<0.7||!KLimit0||kstarmass*kstarmass>1.1) continue;

	 bool XCutCosK = true;
	 if(CutCosK) XCutCosK = cos_theta_k>-0.6&&cos_theta_k<-0.2;
//	 XCutCosK = cos_theta_k>-0.8&&cos_theta_k<-0.2;
//	 XCutCosK = cos_theta_k<-0.8||cos_theta_k>-0.2;
//	 HxMassCosK->Fill(tagged_mass,cos_theta_k);
//	 HxCosK->Fill(cos_theta_k);
	 if(XCutCosK) HxCosK->Fill(cos_theta_k);
	 if(XCutCosK&&!leadK) HxCosK_Lpi->Fill(cos_theta_k);
	 if(XCutCosK&& leadK) HxCosK_LK ->Fill(cos_theta_k);
//
// LEFT SIDEBAND
//
         if(tagged_mass>XMinSBL && tagged_mass<XMaxSBL ){
//  	  HxMassCosK_SBL->Fill(tagged_mass,cos_theta_k);
  	  if(XCutCosK) HxCosK_SBL->Fill(cos_theta_k);
	  if(XCutCosK&&!leadK) HxCosK_SBL_Lpi->Fill(cos_theta_k);
	  if(XCutCosK&& leadK) HxCosK_SBL_LK ->Fill(cos_theta_k);
	  if(XCutCosK) HxBsMass_SBL->Fill(mmpiKaon);
	  if(XCutCosK) HxLambdab_SBL->Fill(lambdab);
	  HxMass_SBL->Fill(tagged_mass);
//	  if(cos_theta_k>-0.6&&cos_theta_k<-0.2){
           Hxmmpikstarmass_SBL->Fill(kstarmassx*kstarmassx,mmpix*mmpix);
           Hxmmpiq_SBL->Fill(mmpix*mmpix);
           Hxkstarmass2_SBL->Fill(kstarmassx*kstarmassx);
//	  }
	 } 
//
// RIGHT SIDEBAND
//
         if(tagged_mass>XMinSBR && tagged_mass<XMaxSBR ){
//  	  HxMassCosK_SBR->Fill(tagged_mass,cos_theta_k);
  	  if(XCutCosK) HxCosK_SBR->Fill(cos_theta_k);
	  if(XCutCosK&&!leadK) HxCosK_SBR_Lpi->Fill(cos_theta_k);
	  if(XCutCosK&& leadK) HxCosK_SBR_LK ->Fill(cos_theta_k);
	  if(XCutCosK) HxLambdab_SBR->Fill(lambdab);
	  if(XCutCosK) HxBsMass_SBR->Fill(mmpiKaon);
	  HxMass_SBR->Fill(tagged_mass);
	 }     
//
//SIGNAL REGION
//
         if(tagged_mass>XMaxSBL && tagged_mass<XMinSBR ){
//  	  HxMassCosK_Sign->Fill(tagged_mass,cos_theta_k);
  	  if(XCutCosK) HxCosK_Sign->Fill(cos_theta_k);
	  if(XCutCosK&&!leadK) HxCosK_Sign_Lpi->Fill(cos_theta_k);
	  if(XCutCosK&& leadK) HxCosK_Sign_LK ->Fill(cos_theta_k);
	  if(XCutCosK) HxBsMass_Sign->Fill(mmpiKaon);
	  if(!XCutCosK) HxBsMass_SignC->Fill(mmpiKaon);
	  HxMass_Sign->Fill(tagged_mass);
//	  if(XCutCosK) HxLambdab_Sign->Fill(lambdab1);
//	  if(XCutCosK) HxLambdab_Sign->Fill(lambdab2);
//          Hxmmpi_Sign->Fill(mmpi1);

//          Hxmmpikstarmass_Sign->Fill(kstarmass*kstarmass,mmpix*mmpix);

//          if(mmpix*mmpix>22) Hxkstarmass2_Sign->Fill(kstarmass*kstarmass);
          Hxmmpikstarmass_Sign->Fill(kstarmassx*kstarmassx,mmpix*mmpix);
          Hxkstarmass2_Sign->Fill(kstarmassx*kstarmassx);
          Hxmmpiq_Sign->Fill(mmpix*mmpix);
	  if(cos_theta_k>-0.6&&cos_theta_k<-0.2){
           Hxmmpikstarmass_Sign_b->Fill(kstarmassx*kstarmassx,mmpix*mmpix);
//           if(!leadK) Hxmmpi_Sign->Fill(mmpi1);
           if( leadK) Hxmmpikstarmass_Sign_LK->Fill(mmpix*mmpix,kstarmassx*kstarmassx);
           if(!leadK) Hxmmpikstarmass_Sign_Lpi->Fill(mmpix*mmpix,kstarmassx*kstarmassx);
           Hxpsi2s_Sign->Fill(mumuMass);
	   HxLambdab_Sign_c->Fill(lambdab);
	  }else{ 
           if( leadK) Hxmmpikstarmass_Sign_LK_c->Fill(mmpix*mmpix,kstarmassx*kstarmassx);
           if(!leadK) Hxmmpikstarmass_Sign_Lpi_c->Fill(mmpix*mmpix,kstarmassx*kstarmassx);
	   HxLambdab_Sign->Fill(lambdab);
           Hxpsi2s_Sign_c->Fill(mumuMass);
	  } 
	 } 
  }
  gStyle->SetPadBorderMode(0);
  gStyle->SetOptStat("e");
//  gStyle->SetOptStat(000000);
  gStyle->SetOptFit(000000);
  HxCosK->SetMinimum(0);
  HxCosK_Sign->SetMinimum(0);
  HxCosK_SBL->SetMinimum(0);
  HxCosK_SBR->SetMinimum(0);
  
//=======================================  
  cStudies->cd(1);
  HxCosK_LK->SetLineColor(kRed);
  HxCosK_Lpi->SetLineColor(kBlue);
  HxCosK->Draw();
  HxCosK_LK->Draw("same");
  HxCosK_Lpi->Draw("same");
  cStudies->cd(2);
//  HxMassCosK_SBL->Draw("contz");
  HxCosK_SBL_LK->SetLineColor(kRed);
  HxCosK_SBL_Lpi->SetLineColor(kBlue);
  HxCosK_SBL->Draw();
  HxCosK_SBL_LK->Draw("same");
  HxCosK_SBL_Lpi->Draw("same");
  cStudies->cd(3);
//  HxMassCosK_Sign->Draw("contz");
  HxCosK_Sign_LK->SetLineColor(kRed);
  HxCosK_Sign_Lpi->SetLineColor(kBlue);
  HxCosK_Sign->Draw();
  HxCosK_Sign_LK->Draw("same");
  HxCosK_Sign_Lpi->Draw("same");
  cStudies->cd(4);
//  HxMassCosK_SBR->Draw("contz");
  HxCosK_SBR_LK->SetLineColor(kRed);
  HxCosK_SBR_Lpi->SetLineColor(kBlue);
  HxCosK_SBR->Draw();
  HxCosK_SBR_LK->Draw("same");
  HxCosK_SBR_Lpi->Draw("same");
  cStudies->cd(5);
  Hxmmpikstarmass_Sign->Draw();
  Hxmmpikstarmass_Sign_b->SetMarkerColor(kRed);
  Hxmmpikstarmass_Sign_b->Draw("same");
  DKinematic->Draw("same");
//  Hxmmpikstarmass_Sign_LK->Draw();
//   HxMass_SBL->Draw();
//   HxBsMass_SBL->Draw("same");
  cStudies->cd(6);
//  Hxmmpikstarmass_Sign_LK_c->Draw();
  Hxmmpiq_Sign->Draw();
  Hxmmpiq_SBL->SetLineColor(kRed);
  Hxmmpiq_SBL->Draw("same");
//   HxMass_Sign->Draw();
//   HxBsMass_Sign->Draw("same");
//   HxBsMass_SignC->Draw("same");
  cStudies->cd(7);
//   HxMass_SBR->Draw();
//   HxBsMass_SBR->Draw("same");
//  Hxmmpikstarmass_Sign_Lpi->Draw();
  Hxkstarmass2_Sign->Draw();
  Hxkstarmass2_SBL->SetLineColor(kRed);
  Hxkstarmass2_SBL->Draw("same");
  cStudies->cd(8);
  Hxmmpikstarmass_SBL->Draw();
  DKinematic->Draw("same");
//   HxLambdab_Sign->Draw();
//   HxLambdab_Sign_c->Draw("same");
//   HxLambdab_Sign_c->SetLineColor(kMagenta);
//   HxLambdab_SBR->SetLineColor(kRed);
//   HxLambdab_SBL->SetLineColor(kBlue);
//   HxLambdab_SBL->Draw("same");
//   HxLambdab_SBR->Draw("same");
  
//  Hxmmpi_Sign->Draw();
//   HxMass->Draw();
//   HxMass_Sign->Draw("same");
//   HxMass_SBL->Draw("same");
//   HxMass_SBR->Draw("same");
//  Hxmmpikstarmass_Sign_Lpi_c->Draw();
//  Hxpsi2s_Sign->Draw();
// 
  sss.clear();
  sss.str("");
  sss<<PNGName<<"-studies"<<".png";
  gSystem->Exec(Form("mv %s %s.tmp",sss.str().c_str(),sss.str().c_str()));
  cStudies->Print(sss.str().c_str());
  sss.clear();
  sss.str("");
//   sss<<PNGName<<"-studies"<<".pdf";
//   gSystem->Exec(Form("mv %s %s.tmp",sss.str().c_str(),sss.str().c_str()));
//   cStudies->Print(sss.str().c_str());
//=======================================  
}  

//=========================================================================================
//
//=========================================================================================
//
// Namelist Routine
//
//=========================================================================================
//
//=========================================================================================
std::map<std::string, std::string> ReadNamelist(int argc, char** argv){
   if ( argc>=1 && (strcmp(argv[0],"namelist")>=0) ){
     std::cout<<"Defined namelist: "<<argv[0]<<std::endl;
   }else{
     std::cout<<"Namelist:"<<argv[0]<<"  should be named/renamed namelist*.list "<<argc<<std::endl;
     exit(1);
   }
   std::vector<std::string> split( char *str, char c = ' ');
   ifstream indata;
   std::map<std::string, std::string> mappa;
   std::string line;
   std::vector<std::string>vstring ;
//
    indata.open(argv[0]);
   if(!indata) { // file couldn't be opened
   	std::cout <<"Line: "<<__LINE__ <<" "<<argv[0]<< " Error: fileList can not be opened" << std::endl;
   	exit(1);
   }
   while(std::getline(indata, line)) {
	 line.erase(std::remove(line.begin(), line.end(), '\t'), line.end());
	 line.erase(std::remove(line.begin(), line.end(), '\n'), line.end());
	 line.erase(std::remove(line.begin(), line.end(), ' ' ), line.end());

 	 char *cstr = new char [line.size()+1];


 	 strcpy (cstr, line.c_str());
//	 cout <<"stringa->"<< cstr << endl;
	 vstring = split(cstr,'=');
	 mappa.insert( std::pair<string,string>(vstring[0],vstring[1]) );
    }
    std::cout<<"//////////////////////////////////////////////////////////////////////////////////////////////////////////////"<<std::endl;
    for (map<string,string>::iterator imap = mappa.begin();
    			       imap != mappa.end();
    			       ++imap)
    {
   	std::cout <<"mappa->"<< (*imap).first<<" = "<<(*imap).second << std::endl;
    }
    std::cout<<"/////////////////////////////////////////////////////////////////////////////////////////////////////////////"<<std::endl;
    indata.close();	
  return mappa ;
}
//===============================================================================================================
std::vector<std::string> split( char *str, char c = ' ')
{
    std::vector<std::string> result;

    while(1)
    {
         char *begin = str;

        while(*str != c && *str)
                str++;

        result.push_back(string(begin, str));

        if(0 == *str++)
                break;
    }

    return result;
}
//===============================================================================================================

std::string grep(std::string& filename, std::string keyword)
{
    std::string out("GREP: NOT FOUND!!!");
    std::string line;
    std::size_t found;
    ifstream in(filename.c_str());
//    std::cout<<"keyword = "<<keyword<<"\n"<<std::endl;
//    std::cout<<"filename= "<<filename<<"\n"<<std::endl;
    if( in.is_open())
    {
	  while( getline(in,line) )
	  {
	      found=line.find(keyword);
//	           std::cout<<"line "<<line<<"\n"<<std::endl;
	      if( found!=std::string::npos)
		   return line ;
	  }
    }else{
     std::cout<<"GREP: can't open filename= "<<filename<<"\n"<<std::endl;
    }
    return out;
}
//===============================================================================================================

bool Background( double m12, double m13, double DmMass, double m1 , double m2, double m3){
//
//     
//
	bool Background = true;
//
	double m12_min     = ( m1 + m2 )*( m1 + m2 );	
	double m12_max     = ( DmMass - m3 )*( DmMass - m3 );	
	double m13_tempmin = ( m1 + m3 )*( m1 + m3 );	
	double m13_tempmax = ( DmMass - m2 )*( DmMass - m2 );	
//
	if( m12 < m12_min )     Background = false;
	if( m12 > m12_max )     Background = false;
	if( m13 < m13_tempmin ) Background = false;
	if( m13 > m13_tempmax ) Background = false;
//	
	double E1 = ( m12                   +m1*m1-m2*m2)/(2.* sqrt( m12 ));
	double E3 = ( DmMass*DmMass-m12  -m3*m3)/(2.* sqrt( m12 )) ;
	double A  = E1*E1 - m1*m1 ;
	double B  = E3*E3 - m3*m3 ;
	if( A < 0 ) A = 0.;
	if( B < 0 ) B = 0.;
	double m13_min = ( E1 + E3 )*( E1 + E3 ) - ( sqrt( A ) + sqrt( B ))*( sqrt( A ) + sqrt( B ));
	double m13_max = ( E1 + E3 )*( E1 + E3 ) - ( sqrt( A ) - sqrt( B ))*( sqrt( A ) - sqrt( B ));
	if( m13 < m13_min ) Background = false;
	if( m13 > m13_max ) Background = false;
	
//
	return Background;
}
//---------------------------------------------------------------------------------------	
//---------------------------------------------------------------------------------------	
  double KineLimit( double* x,  double *par){
  double BackgroundPlot( double m12, double m13, double DmMass, double m1 , double m2, double m3);
  double m12 =x[0];
  double m13 =x[1];
  double DMass = par[0];
  double m1 = par[1];
  double m2 = par[2];
  double m3 = par[3];
//   double piMass = 0.13957039;
//   double kMass = 0.493677;
//   double BpMass = 5.2791;
//   double B0Mass = 5.27962;
//  double KstarMass = 0.896;
//  double KstarMass = 0.892;
//   double MassPsi2S=3.696;
//   double MassJPsi =3.0969;
  double ret = BackgroundPlot(  m12,  m13,  DMass,  m1 ,  m2, m3);
//  if ((m12>=m13)) ret=1.;
//  cout<<Form("%f vs %f = %f",m12,  m13, ret)<<endl;
  return ret;
}
//===============================================================================================================
double FunCutkst2(double* x, double *par){

double  BpMass = 5.2791;
double  KstarMass = 0.896;
//double  KstarMass = 0.892;
double  Deltakst2     = x[0];
double  DeltammpiKaon = x[1];
double  Fx0Cut =par[0];
double  Fx1Cut =par[1];
double  Fy0Cut =par[2];
double  Fy1Cut =par[3];

double ret =0.;

//diagonal
if ( fabs(((DeltammpiKaon)-Fy0Cut)/(Fy1Cut-Fy0Cut)-((Deltakst2)-Fx0Cut)/(Fx1Cut-Fx0Cut))<0.001 && Deltakst2>0.001 ) ret=1.;

//vertical line
if ( (((DeltammpiKaon)-Fy0Cut)/(Fy1Cut-Fy0Cut)<((Deltakst2)-Fx0Cut)/(Fx1Cut-Fx0Cut)) && fabs(Deltakst2)<0.001 ) ret=1.;
//if (fabs(Deltakst2)<0.0001&) ret=1.;

return ret; 

}	
//===============================================================================================================
double FunCutmmpi2(double* x, double *par){
double  BpMass = 5.2791;
double  KstarMass = 0.896;
//double  KstarMass = 0.892;
double  Fmmpi2 = x[0];
double  Fmmka1 = x[1];
double  FCutX1 =par[0];
double  FCutX2 =par[1];
double  FCutY1 =par[2];
double  FCutY2 =par[3];

double  Fx_0Cut=par[4];
double  Fx_1Cut=par[5];
double  Fy_0Cut=par[6];
double  Fy_1Cut=par[7];
double  ret =0.;

double xc =((FCutY1-Fy_0Cut)/(Fy_1Cut-Fy_0Cut))*(Fx_1Cut-Fx_0Cut)+Fx_0Cut;
//cout<<"xc="<<xc<<endl;
double yc =((FCutX2-Fx_0Cut)/(Fx_1Cut-Fx_0Cut))*(Fy_1Cut-Fy_0Cut)+Fy_0Cut;
//cout<<"yc="<<yc<<endl;
//exit(1);

if ( (Fmmpi2>FCutX1&&Fmmpi2<xc)&&fabs(Fmmka1-FCutY1)<0.001 &&((Fmmka1-Fy_0Cut)/(Fy_1Cut-Fy_0Cut))>((Fmmpi2-Fx_0Cut)/(Fx_1Cut-Fx_0Cut))) ret =1.;
if ( (Fmmpi2>FCutX1&&Fmmpi2<FCutX2)&&fabs(Fmmka1-FCutY2)<0.001 &&((Fmmka1-Fy_0Cut)/(Fy_1Cut-Fy_0Cut))>((Fmmpi2-Fx_0Cut)/(Fx_1Cut-Fx_0Cut))) ret =1.;

if (fabs(Fmmpi2-FCutX1)<0.001&&(Fmmka1>FCutY1&&Fmmka1<FCutY2)  &&((Fmmka1-Fy_0Cut)/(Fy_1Cut-Fy_0Cut))>((Fmmpi2-Fx_0Cut)/(Fx_1Cut-Fx_0Cut)))  ret=1.;
if (fabs(Fmmpi2-FCutX2)<0.001&&(Fmmka1>yc&&Fmmka1<FCutY2)  &&((Fmmka1-Fy_0Cut)/(Fy_1Cut-Fy_0Cut))>((Fmmpi2-Fx_0Cut)/(Fx_1Cut-Fx_0Cut)))  ret=1.;



//if (fabs((Fmmka1-Fy_0Cut)/(Fy_1Cut-Fy_0Cut))-((Fmmpi2-Fx_0Cut)/(Fx_1Cut-Fx_0Cut))<0.00001&&Fmmpi2>xc&&Fmmka1>FCutY1)  ret=1.;

if (((Fmmka1-Fy_0Cut)/(Fy_1Cut-Fy_0Cut))>((Fmmpi2-Fx_0Cut)/(Fx_1Cut-Fx_0Cut))&&((Fmmka1-Fy_0Cut)/(Fy_1Cut-Fy_0Cut))-((Fmmpi2-Fx_0Cut)/(Fx_1Cut-Fx_0Cut))<0.01&&(Fmmpi2>xc&&Fmmpi2<FCutX2)&&(Fmmka1>FCutY1&&Fmmka1<yc))  ret=1.;
//if (fabs((Fmmka1-Fy_0Cut)/(Fy_1Cut-Fy_0Cut))-((Fmmpi2-Fx_0Cut)/(Fx_1Cut-Fx_0Cut))<0.01)  ret=1.;
return ret; 
}	
//---------------------------------------------------------------------------------------	
double BackgroundPlot( double m12, double m13, double DmMass, double m1 , double m2, double m3){
//
//     
//
	double Background = 1.;
//
	double m12_min     = ( m1 + m2 )*( m1 + m2 );	
	double m12_max     = ( DmMass - m3 )*( DmMass - m3 );	
	double m13_tempmin = ( m1 + m3 )*( m1 + m3 );	
	double m13_tempmax = ( DmMass - m2 )*( DmMass - m2 );	
//
	if( m12 < m12_min )     Background = 0.;
	if( m12 > m12_max )     Background = 0.;
	if( m13 < m13_tempmin ) Background = 0.;
	if( m13 > m13_tempmax ) Background = 0.;
//	
	double E1 = ( m12                   +m1*m1-m2*m2)/(2.* sqrt( m12 ));
	double E3 = ( DmMass*DmMass-m12  -m3*m3)/(2.* sqrt( m12 )) ;
	double A  = E1*E1 - m1*m1 ;
	double B  = E3*E3 - m3*m3 ;
	if( A < 0 ) A = 0.;
	if( B < 0 ) B = 0.;
	double m13_min = ( E1 + E3 )*( E1 + E3 ) - ( sqrt( A ) + sqrt( B ))*( sqrt( A ) + sqrt( B ));
	double m13_max = ( E1 + E3 )*( E1 + E3 ) - ( sqrt( A ) - sqrt( B ))*( sqrt( A ) - sqrt( B ));
	if( m13 < m13_min ) Background = 0.;
	if( m13 > m13_max ) Background = 0.;
	
//
	return Background;
}
//===============================================================================================================
void replaceAll(std::string& str, const std::string& from, const std::string& to) {
    if(from.empty())
        return;
    size_t start_pos = 0;
    while((start_pos = str.find(from, start_pos)) != std::string::npos) {
        str.replace(start_pos, from.length(), to);
        start_pos += to.length(); // In case 'to' contains 'from', like replacing 'x' with 'yx'
    }
}

//===============================================================================================================
void replaceChar(char * txt,const  char * txt1,const  char * txt2) {

  std::stringstream sss,sss1,sss2;
  sss<<txt;
  sss1<<txt1;
  sss2<<txt2;
  std::string ss=sss.str();
  replaceAll( ss,  sss1.str(), sss2.str());
  strcpy(txt,ss.c_str());
  sss.str("");
  sss.clear();
  sss1.str("");
  sss1.clear();
  sss2.str("");
  sss2.clear();
  printf ("replaceChar output=>%s\n",txt);
}  

