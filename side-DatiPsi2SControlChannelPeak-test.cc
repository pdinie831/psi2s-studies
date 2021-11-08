//
//g++ -O3 -o side-DatiPsi2SControlChannelPeak-test side-DatiPsi2SControlChannelPeak-test.cc `root-config --cflags --libs`  -lRooFit -lGenVector -lMathCore
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
#include <cmath>
#include <TMath.h>
#include <TLine.h>
#include <TArrow.h>
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
#include <Math/Math.h>
#include <Math/GenVector/GenVector_exception.h>
#include <Math/Vector4Dfwd.h>
#include <Math/GenVector/LorentzVector.h>
#include <Math/GenVector/PtEtaPhiM4D.h>//#include <RooMath/GenVector/PtEtaPhiM4D.h>
#include <Math/GenVector/Boost.h>

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
#include "TPaveLabel.h"

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
void AnalizeDatiPsi2S_test();
double computeCosine (double Vx,double Vy,double Vz,double Wx, double Wy, double  Wz);
void replaceAll(std::string& str, const std::string& from, const std::string& to) ;
void replaceChar(char * txt, const char * txt1, const char * txt2) ;
bool Background( double m12, double m13, double DmMass, double m1 , double m2, double m3);
double  KineLimit(  double *x, double *par);
double  FunCutkst2( double *x, double *par);
double  FunCutmmpi2( double *x, double *par);
double  BackgroundPlot( double m12, double m13, double DmMass, double m1 , double m2, double m3);
double  P_OnShell_Compute( double XM12, double XM1,double  XM2 );
double  vettore(double XDmMass,double Am12,double Am13, double Am1,double Am2, double Am3, double Vemas, double GamVe);
double  DalitzModel(  double *x, double *par);
double  DalitzModel2(  double *x, double *par);
double lambda(double x,double y,double z);
double cos12( double DmMass,double m1,double m2,double m3,double m12,double m13);
double BreitWigner( double* x,  double *par);

std::string grep(std::string& filename, std::string keyword);
//
TFile*OutFileInputHisto;

//
char TYPEFILE[10]=".png";
char TXTMASSW[100];

char * RunEra;
char OutFileNameInputHistoNoKstar[350] = "/gwpool/users/dini/NeuralNetwork/Keras/sub_samples/all_2018_data_Charmonium_addDNN.root.scaled";
//char OutFileNameInputHistoNoKstar[350] = "/gwpool/users/dini/NeuralNetwork/Keras/sub_samples/all_2018_data_Charmonium_addDNN_test2.root.scaled";
//char OutFileNameInputHistoNoKstar[350] = "/gwpool/users/dini/NeuralNetwork/Keras/sub_samples/all_2018_data_Charmonium_addDNN.root.final";
//char OutFileNameInputHistoNoKstar[350] = "/gwpool/users/dini/NeuralNetwork/Keras/sub_samples/all_2018_data_Charmonium_addDNN.root.ok8";
																       //char OutFileNameInputHistoNoKstar[350] = "/gwpool/users/dini/NeuralNetwork/Keras/sub_samples/all_2018_data_Charmonium_addDNN.root.ok10";
//char OutFileNameInputHistoNoKstar[350]        = "/gwteray/users/fiorendi/p5prime/data2018/flat_ntuples/Oct2020/2018D_noKstar.root";
char OutFileNameInputHisto[300]        = "/gwpool/users/dini/p5prime/2018/skims/newphi/2018Data_All_finalSelection.root";
char OutFileNameInputHistoMC[300]      = "/gwpool/users/dini/p5prime/2018/skims/newphi/2018MC_JPSI.root";
char OutFileNameInputHistoMCNoKstar[350] = "/gwpool/users/dini/NeuralNetwork/Keras/sub_samples/sample_2018_MC_Charmonium_all_addDNN.root";
//char OutFileNameInputHistoNoKstar[300]      = "/gwpool/users/dini/p5prime/2018/skims/newphi/2018MC_PSI.root";
char OutFileNameInputHistoMCPsi2S[300] = "/gwpool/users/dini/p5prime/sidebands/2018MC_BToPsi2SK.root";

char OutputRecoB0TreeName[10]     = "ntuple";
char PDFName[300]="AnalizeDatiPsi2SPeak_test2018";
char NameList[300]="xcut-namelist.lis";

std::stringstream sss;
double XMinSign = 5.0;
double XMaxSign = 5.6;
//
double XMinSBL  = XMinSign;
double XMaxSBL  = 5.15;
double XMinSBR  = 5.4;
double XMaxSBR  = XMaxSign;

double WMinSign= XMaxSBL;
double WMaxSign= XMinSBR;
//
double dnn_cut = 0.85;
//
double XStepSign = 0.0025*2;
float xMassHBin = (XMaxSign -XMinSign)/XStepSign;
float xMassHBin2   =  xMassHBin ; // plot only!
//
bool MC = false;
bool NK = false;
bool INV = false;
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
       replaceChar(PDFName,"2018",RunEra);
       replaceChar(OutFileNameInputHisto  ,"2018",RunEra);
       replaceChar(OutFileNameInputHistoMC,"2018",RunEra);
       replaceChar(OutFileNameInputHistoMCPsi2S,"2018",RunEra);
     }else if ((strcmp(argv[1],"2017") == 0)){
       std::cout<<Form("==>Setting Era = %s\n",RunEra)<<std::endl;
       replaceChar(PDFName,"2018",RunEra);
       replaceChar(OutFileNameInputHisto  ,"2018",RunEra);
       replaceChar(OutFileNameInputHistoMC,"2018",RunEra);
       replaceChar(OutFileNameInputHistoMCPsi2S,"2018",RunEra);
     }else if ((strcmp(argv[1],"2018") == 0)){
       std::cout<<Form("==>Setting Era = %s\n",RunEra)<<std::endl;
     }else{
      std::cout<<Form("Error: Era %s not recognize => set the ERA [2016 2017 2018]\n",RunEra)<<std::endl;
      exit(1);
     }
     if (argc>2 && ((strcmp(argv[2],"m") == 0)||(strcmp(argv[2],"M") == 0)) ){
      cout<<"Read MC  !!!" <<endl;
      sprintf(OutFileNameInputHisto,"%s",OutFileNameInputHistoMC);
      replaceChar(PDFName,"Dati","MC");
      MC = true;
      if (argc>3 && ((strcmp(argv[3],"tag") == 0)||(strcmp(argv[3],"TAG") == 0)) ){
     	replaceChar(PDFName,"MC","MCCorreTag");
     	CorreTag=true;
        std::cout<<Form("==>Setting MC Correctly Tagged RunEra=%s\n",RunEra)<<std::endl;
      } else if (argc>3 && ((strcmp(argv[3],"mis") == 0)||(strcmp(argv[3],"MIS") == 0)) ){
     	replaceChar(PDFName,"MC","MCWrongTag");
     	WrongTag=true;
        std::cout<<Form("==>Setting MC Wrongly   Tagged RunEra=%s\n",RunEra)<<std::endl;
      }else{
        std::cout<<Form("==>Setting MC All Tags RunEra=%s\n",RunEra)<<std::endl;
      }
 
     }else if (argc>2 && ((strcmp(argv[2],"nk") == 0)||(strcmp(argv[2],"NK") == 0)) ){
 
      cout<<"Read No Kstar Data Ntupla !!!" <<endl;
//   	if (argc>2 && ((strcmp(argv[2],"tag") == 0)||(strcmp(argv[2],"TAG") == 0)) ) CorreTag=true;
//   	if (argc>2 && ((strcmp(argv[2],"mis") == 0)||(strcmp(argv[2],"MIS") == 0)) ) WrongTag=true;
      sprintf(OutFileNameInputHisto,"%s",OutFileNameInputHistoNoKstar);
      replaceChar(PDFName,"Dati","DatiNoKstar");
      XStepSign=XStepSign*2;
      xMassHBin = xMassHBin/2;
      xMassHBin2 = xMassHBin;
      NK = true;
      if (argc>3 && ((strcmp(argv[3],"inv") == 0)||(strcmp(argv[3],"INV") == 0)) ){
     	replaceChar(PDFName,"Dati","DatiInverted");
     	INV=true;
        std::cout<<Form("==>Setting Data swap pi<->K RunEra=\n",RunEra)<<std::endl;
      } 
     }else if (argc>2 && ((strcmp(argv[2],"mnk") == 0)||(strcmp(argv[2],"MNK") == 0)) ){
 
      cout<<"===>>>>> MC!!! Read No Kstar MC Ntupla !!!" <<endl;
      sprintf(OutFileNameInputHisto,"%s",OutFileNameInputHistoMC);
      replaceChar(PDFName,"Dati","MC");
      NK = true;
      MC = true;
      if (argc>3 && ((strcmp(argv[3],"tag") == 0)||(strcmp(argv[3],"TAG") == 0)) ){
     	replaceChar(PDFName,"MC","MCCorreTag");
     	CorreTag=true;
        std::cout<<Form("==>Setting MC Correctly Tagged RunEra=%s\n",RunEra)<<std::endl;
      } else if (argc>3 && ((strcmp(argv[3],"mis") == 0)||(strcmp(argv[3],"MIS") == 0)) ){
     	replaceChar(PDFName,"MC","MCWrongTag");
     	WrongTag=true;
        std::cout<<Form("==>Setting MC Wrongly   Tagged RunEra=%s\n",RunEra)<<std::endl;
      }else{
        std::cout<<Form("==>Setting MC All Tags RunEra=%s\n",RunEra)<<std::endl;
      }
//   	if (argc>2 && ((strcmp(argv[2],"tag") == 0)||(strcmp(argv[2],"TAG") == 0)) ) CorreTag=true;
//   	if (argc>2 && ((strcmp(argv[2],"mis") == 0)||(strcmp(argv[2],"MIS") == 0)) ) WrongTag=true;
      sprintf(OutFileNameInputHisto,"%s",OutFileNameInputHistoMCNoKstar);
      replaceChar(PDFName,"MC","MCNoKstar");
      XStepSign=XStepSign*2;
      xMassHBin = xMassHBin/2;
      xMassHBin2 = xMassHBin;
      if (argc>3 && ((strcmp(argv[3],"inv") == 0)||(strcmp(argv[3],"INV") == 0)) ){
     	replaceChar(PDFName,"MC","MCInverted");
     	INV=true;
        std::cout<<Form("==>Setting MC swap pi<->K RunEra=\n",RunEra)<<std::endl;
      } 
 
     }else{
      cout<<"Read Data!!!" <<endl;
      MC = false;
      NK = false;
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
  
  map<string,string>::iterator  iw= mappa.find("WMinSign");
  if(iw != mappa.end()) {
   WMinSign  =  atof (mappa["WMinSign"].c_str() ) ;
   std::cout<<Form(" Setting WMinSign = %f \n",WMinSign)<<std::endl;
  }
  map<string,string>::iterator  is= mappa.find("WMaxSign");
  if(is != mappa.end()) {
   WMaxSign  =  atof (mappa["WMaxSign"].c_str() ) ;
   std::cout<<Form(" Setting WMaxSign = %f \n",WMaxSign)<<std::endl;
  }
  map<string,string>::iterator  id= mappa.find("WMaxSign");
  if(id != mappa.end()) {
   dnn_cut  =  atof (mappa["dnn_cut"].c_str() ) ;
   std::cout<<Form(" Setting dnn_cut = %f \n",dnn_cut)<<std::endl;
  }
  if(Cutkst2 ) {
   std::cout<<" Setting Cut kst2 \n"<<std::endl;
   sprintf(PDFName,"%s-Cutkst2",PDFName);
  } 
  if(Cutmmpi2) {
   std::cout<<" Setting Cut mmpi2\n"<<std::endl;
   sprintf(PDFName,"%s-Cutmmpi2",PDFName);
  } 
  if(CutCosK) {
   std::cout<<" Setting Cut for plot CutCosK\n"<<std::endl;
   sprintf(PDFName,"%s-ScatterPlotCutCosK",PDFName);
  } 
  
  if(WMinSign<XMaxSBL) {
   std::cout<<"Error!!! WMinSign<XMinSBR\n"<<std::endl;
   exit(1);
  } 
  if(WMaxSign>XMinSBR) {
   std::cout<<"Error!!! WMaxSign>XMinSBR\n"<<std::endl;
   exit(1);
  } 
  
  if(NK){  
      char *DNNTXT = (Form("dnncut%3.2f",dnn_cut));
      replaceChar(DNNTXT,".","p");
      cout<< DNNTXT<<endl;
      sprintf(PDFName,Form("%s-%s",PDFName,DNNTXT));
  }    
  
  sprintf(TXTMASSW,"-SBL-%3.2f-%3.2f-SGN-%3.2f-%3.2f-SBR-%3.2f-%3.2f",XMinSBL,XMaxSBL,WMinSign,WMaxSign,XMinSBR,XMaxSBR);
  replaceChar(TXTMASSW,".","p");
  TStopwatch TimeWatch;
  TimeWatch.Start();

  AnalizeDatiPsi2S_test(); 
//  app.Run() ;
  cout<<"esco..." <<endl;
  TimeWatch.Stop();
  TimeWatch.Print();
  
  return 0 ;
}
void  AnalizeDatiPsi2S_test(){


  double cos_theta_l	  ;
  double cos_theta_k	  ;
  double phi_kst_mumu	  ;
  double tagged_mass	  ;
  double bMass		  ;
  double bBarMass	  ;
  double bPt		  ;
  double bEta		  ;
  double bPhi		  ;
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
  double KstarMassWT = 1.;
  double KstarGamm = 0.047;  
  double ZMass = 4.478;
  double ZGamm = 0.181;
  double RhoMass=0.775;
  double RhoGam=0.148;
//  double KstarMass = 0.892;
  double gen_jpsimass_12 ;
  double mumTrkp ;
  double mupTrkm ;
  double mmk1    ;
  double mmk2    ;
  double kstarmass ;
  double kstPt	   ;
  double kstEta	   ;
  double kstPhi	   ;
  double kstMass    ;
  double kstBarMass ;
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
  double bLBS,bLBSE,bVtxCL;
  double kstTrkmDCABS,kstTrkpDCABS,kstTrkmDCABSE,kstTrkpDCABSE;
  double kstTrkmMinIP2D, kstTrkpMinIP2D;
  double kkMass;
  double bCosAlphaBS;
  float  dnn_prob;
  bool   passB0Psi_jpsi ;
  bool   passB0Psi_psip ;
  bool   passB0Psi_lmnr;
  double  genSignal;
  TCanvas* cStudies = new TCanvas("cStudies","Cutmmpi2&Cutkst2",200,10,900,780);
  cStudies->Divide(4,2);
  TCanvas* cDalitz      = new TCanvas("cDalitz","#psi(2S)#pi^{2}%#psi(2S)K^{2}",200,10,900,780);
  TCanvas* cDalitzSBL   = new TCanvas("cDalitzSBL","#psi(2S)#pi^{2}%#psi(2S)K^{2}",200,10,900,780);
  TCanvas* cDaliCut     = new TCanvas("cDaliCut","#psi(2S)#pi^{2}%#psi(2S)K^{2} Cutted",200,10,900,780);
  TCanvas* cHely    = new TCanvas("cHely","Helicity",200,10,900,780);
  TCanvas* cHelyLeft= new TCanvas("cHelyLeft","Helicity Left",200,10,900,780);
  TCanvas* ckstarproj= new TCanvas("ckstarproj","Helicity Left",200,10,900,780);
  TCanvas* cMass = new TCanvas("cMass","Mass spectrum",200,10,900,780);
  TCanvas* cPeak = new TCanvas("cPeak","peak studies",200,10,1200,900);
  TCanvas* cLambda = new TCanvas("cLambda","Lambda studies",200,10,1200,780);
  TCanvas* cBelle = new TCanvas("cBelle","Dalitz slices",200,10,1400,600);
  cPeak->Divide(3,2);
  cBelle->Divide(4,2);

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
//      TH1D* Hxmmpi_Sign	     = new TH1D( "Hxmmpi_Sign"	    , "psi2spi mass"	                             ,48, 14, 26);
//      TH1D* Hxmmpiq_SBL	     = new TH1D( "Hxmmpiq_SBL"	    , "(psi2spi mass)^{2}"	                             ,48, 14, 26);
//      TH1D* Hxmmpiq_Sign	     = new TH1D( "Hxmmpiq_Sign"	    , "(psi2spi mass)^{2}"	                             ,48, 14, 26);
     TH1D* HxLambdab_Sign	     = new TH1D( "HxLambdab_Sign"	    , "Lambdab Mass"	                             ,80, 5.45, 5.9);
     TH1D* HxLambdab_Sign_c	     = new TH1D( "HxLambdab_Sign_c"	    , "Lambdab Mass"	                             ,80, 5.45, 5.9);
     TH1D* HxLambdab_SBL	     = new TH1D( "HxLambdab_SBL"	    , "Lambdab Mass"	                             ,80, 5.45, 5.9);
     TH1D* HxLambdab_SBR	     = new TH1D( "HxLambdab_SBR"	    , "Lambdab Mass"	                             ,80, 5.45, 5.9);
     TH1D* HxLambdab_Signc	     = new TH1D( "HxLambdab_Signc"	    , "Lambdab Mass c"	                             ,80, 5.45, 5.9);
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
     TH1D* HxHely_Sign	     = new TH1D( "HxHely_Sign"      , "B^{0} K* Helicity (signal)"  ,	 120, -1.2, 1.2);
     TH1D* HxCosK_Sign_H     = new TH1D( "HxCosK_Sign_H"    , "B^{0} CosK (signal)"	    ,	 120, -1.2, 1.2);
     TH1D* HxCosK_Sign_R     = new TH1D( "HxCosK_Sign_R"    , "B^{0} CosK (signal)"	    ,	 120, -1.2, 1.2);
     TH1D* HxCosK_Sign_W     = new TH1D( "HxCosK_Sign_W"    , "B^{0} CosK (signal)"	    ,	 120, -1.2, 1.2);
     TH1D* HxHely_Left	     = new TH1D( "HxHely_Left"      , "B^{0} K* Helicity (Left)"  ,	 120, -1.2, 1.2);
     TH1D* HxCosK_Left_H     = new TH1D( "HxCosK_Left_H"    , "B^{0} CosK (Left)"	    ,	 120, -1.2, 1.2);
     TH1D* HxCosK_Left_R     = new TH1D( "HxCosK_Left_R"    , "B^{0} CosK (Left)"	    ,	 120, -1.2, 1.2);
//     
     int scala = 3 ;  
//     
     TH1D* HxCosK	     = new TH1D( "HxCosK"           , "B^{0} CosK"		    ,	 10*scala, -1., 1.);
     TH1D* HxCosK_LK	     = new TH1D( "HxCosK_Lpi"       , "B^{0} CosK leading pi"	    ,	 10*scala, -1., 1.);
     TH1D* HxCosK_Lpi	     = new TH1D( "HxCosK_LK"        , "B^{0} CosK leading K"	    ,	 10*scala, -1., 1.);
     TH1D* HxCosK_SBR	     = new TH1D( "HxCosK_SBR"	    , "B^{0} cos(#theta_{K}) SB Right"		    ,	 10*scala, -1., 1.);
     TH1D* HxCosK_SBL	     = new TH1D( "HxCosK_SBL"       , "B^{0} cos(#theta_{K}) SB Left"		    ,	 10*scala, -1., 1.);
     TH1D* HxCosK_Sign	     = new TH1D( "HxCosK_Sign"      , Form("cos(#theta_{K}) (B^{0} signal: %3.2f<mass<%3.2f GeV/c^{2})",WMinSign,WMaxSign)	    ,	 10*scala, -1., 1.);
     TH1D* HxCosK_SBR_Lpi    = new TH1D( "HxCosK_SBR_Lpi"   , "B^{0} cos(#theta_{K}) SBR leading pi"	   ,	10*scala, -1., 1.);
     TH1D* HxCosK_SBL_Lpi    = new TH1D( "HxCosK_SBL_Lpi"   , "B^{0} cos(#theta_{K}) SBL leading pi"	   ,	10*scala, -1., 1.);
     TH1D* HxCosK_Sign_Lpi   = new TH1D( "HxCosK_Sign_Lpi"  , "B^{0} CosK (signal) leading pi"	   ,	10*scala, -1., 1.);
     TH1D* HxCosK_SBR_LK     = new TH1D( "HxCosK_SBR_LK"    , "B^{0} cos(#theta_{K}) SBR leading K"	   ,	10*scala, -1., 1.);
     TH1D* HxCosK_SBL_LK     = new TH1D( "HxCosK_SBL_LK"    , "B^{0} cos(#theta_{K}) SBL leading K"	   ,	10*scala, -1., 1.);
     TH1D* HxCosK_Sign_LK    = new TH1D( "HxCosK_Sign_LK"   , "B^{0} CosK (signal) leading K"      ,	10*scala, -1., 1.);
     TH1D* HxCosK_Sign_Ks    = new TH1D( "HxCosK_Sign_Ks"      , "B^{0} CosK (signal) Kstar band"  ,	 10*scala, -1., 1.);
     TH1D* HxCosK_Sign_NK    = new TH1D( "HxCosK_Sign_NK"      , "B^{0} CosK (signal) Ktsar veto"	   ,	 10*scala, -1., 1.);
     TH1D* HxCosK_Sign_Md    = new TH1D( "HxCosK_Sign_Md"      , "B^{0} CosK (signal) 1<K#pi mass^{2}<1.75",	 10*scala, -1., 1.);
     TH1D* HxCosK_Sign_Tl    = new TH1D( "HxCosK_Sign_Tl"      , "B^{0} CosK (signal) K#pi mass^{2}>1.75",	 5*scala, -1., 1.);
     TH2D* Hxmmpikstarmass_Sign_LK    = new TH2D( "Hxmmpikstarmass_Sign_LK"   , "#psi(2S)#pi mass vs kstar mass lead K",          100, 14,26,100,0.3,2.2);
     TH2D* Hxmmpikstarmass_Sign_LK_c  = new TH2D( "Hxmmpikstarmass_Sign_LK_c" , "#psi(2S)#pi mass vs kstar mass lead K comp",     100, 14,26,100,0.3,2.2);
     TH2D* Hxmmpikstarmass_Sign_Lpi   = new TH2D( "Hxmmpikstarmass_Sign_Lpi"  , "#psi(2S)#pi mass vs kstar mass lead pi",         100, 14,26,100,0.3,1.2);
     TH2D* Hxmmpikstarmass_Sign_Lpi_c = new TH2D( "Hxmmpikstarmass_Sign_Lpi_c", "#psi(2S)#pi mass vs kstar mass lead pi comp",    100, 14,26,100,0.3,1.2);
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

  TH2D* HxMassCosK      = new TH2D( "HxMassCosK"     , "B^{0} Mass%cos(#theta_{K})"       ,	      60,  5.,5.6, 100,  -1.,1.);
  TH2D* HxMassCosK_Sign = new TH2D( "HxMassCosK_Sign", "B^{0} Mass%cos(#theta_{K}) Sign"  ,	      60,  5.,5.6, 100,  -1.,1.);
  TH2D* HxMassCosK_SBL  = new TH2D( "HxMassCosK_SBL" , "B^{0} Mass%cos(#theta_{K}) SBL"   ,	      60,  5.,5.6, 100,  -1.,1.);
  TH2D* HxMassCosK_SBR  = new TH2D( "HxMassCosK_SBR" , "B^{0} Mass%cos(#theta_{K}) SBR"   ,	      60,  5.,5.6, 100,  -1.,1.);

  int nentries = (int)RecoB0TreeOut->GetEntries(); 
  cout<<"nentries: "<< nentries<<endl;
  if (!NK){
  RecoB0TreeOut->SetBranchAddress("passB0Psi_lmnr" 	 ,&passB0Psi_lmnr );
  RecoB0TreeOut->SetBranchAddress("passB0Psi_jpsi" 	 ,&passB0Psi_jpsi );
  RecoB0TreeOut->SetBranchAddress("passB0Psi_psip" 	 ,&passB0Psi_psip );
  RecoB0TreeOut->SetBranchAddress("mumTrkp"		 ,&mumTrkp );
  RecoB0TreeOut->SetBranchAddress("mupTrkm"		 ,&mupTrkm );
  }else{
   RecoB0TreeOut->SetBranchAddress("dnn_prob"             ,&dnn_prob);
  } 
  RecoB0TreeOut->SetBranchAddress("kstarmass"		 ,&kstarmass );
  RecoB0TreeOut->SetBranchAddress("tagged_mass" 	 ,&tagged_mass );
  RecoB0TreeOut->SetBranchAddress("bMass"                ,&bMass);
  RecoB0TreeOut->SetBranchAddress("bBarMass"             ,&bBarMass);
  RecoB0TreeOut->SetBranchAddress("tagB0"                ,&tagB0);
  RecoB0TreeOut->SetBranchAddress("cos_theta_l" 	 ,&cos_theta_l );
  RecoB0TreeOut->SetBranchAddress("cos_theta_k" 	 ,&cos_theta_k );
  RecoB0TreeOut->SetBranchAddress("phi_kst_mumu"	 ,&phi_kst_mumu);
  RecoB0TreeOut->SetBranchAddress("bPt" 	 	 ,&bPt );
  RecoB0TreeOut->SetBranchAddress("bEta" 	 	 ,&bEta );
  RecoB0TreeOut->SetBranchAddress("bPhi" 	 	 ,&bPhi );
  RecoB0TreeOut->SetBranchAddress("mumuMass"		 ,&mumuMass );
  RecoB0TreeOut->SetBranchAddress("mumuPt"		 ,&mumuPt );
  RecoB0TreeOut->SetBranchAddress("mumuPhi"		 ,&mumuPhi );
  RecoB0TreeOut->SetBranchAddress("mumuEta"		 ,&mumuEta );
  RecoB0TreeOut->SetBranchAddress("mumuMassE"		 ,&mumuMassE );
  RecoB0TreeOut->SetBranchAddress("mmk1"		 ,&mmk1 );
  RecoB0TreeOut->SetBranchAddress("mmk2"		 ,&mmk2 );
  RecoB0TreeOut->SetBranchAddress("kstPt" 	 	 ,&kstPt );
  RecoB0TreeOut->SetBranchAddress("kstEta" 	 	 ,&kstEta );
  RecoB0TreeOut->SetBranchAddress("kstPhi" 	 	 ,&kstPhi );
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
  RecoB0TreeOut->SetBranchAddress("bLBS"		 ,&bLBS);
  RecoB0TreeOut->SetBranchAddress("bLBSE"		 ,&bLBSE);
  RecoB0TreeOut->SetBranchAddress("kkMass"		 ,&kkMass);
  RecoB0TreeOut->SetBranchAddress("bCosAlphaBS"		 ,&bCosAlphaBS);
  RecoB0TreeOut->SetBranchAddress("bVtxCL"		 ,&bVtxCL);
  RecoB0TreeOut->SetBranchAddress("kstTrkmDCABS"	 ,&kstTrkmDCABS);
  RecoB0TreeOut->SetBranchAddress("kstTrkmDCABSE"	 ,&kstTrkmDCABSE);
  RecoB0TreeOut->SetBranchAddress("kstTrkpDCABS"	 ,&kstTrkpDCABS);
  RecoB0TreeOut->SetBranchAddress("kstTrkpDCABSE"        ,&kstTrkpDCABSE);
  RecoB0TreeOut->SetBranchAddress("kstTrkmMinIP2D"       ,&kstTrkmMinIP2D);
  RecoB0TreeOut->SetBranchAddress("kstTrkpMinIP2D"       ,&kstTrkpMinIP2D);
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
  double kstarmassInvAll;
  double tagged_massInvAll   ;
  double tagged_massAll      ;
  double mmpix;
  double mmkk;
  double masskk;
  double masspipi;
  double Kstarmum;
  double Kstarmup;
  double hely_ang;
  double kstarmass_rnd,mmpi2_rnd,hely_rnd,hely_wnd,hely_zed,hely_zmd;
  
  double DmMass = B0Mass;
//   double m1	= piMass;
//   double m2	= kMass;
//   double m3	= MassPsi2S;
//   double m1	= KstarMass;
//   double m2	= muonMass_;
//   double m3	= muonMass_;
  double m3	= MassPsi2S;
  double m2	= kMass;
  double m1	= piMass;
  float m12max=(DmMass-m3)*(DmMass-m3);
  float m13min=(m1+m3)*(m1+m3);
  float m13max=(DmMass-m2)*(DmMass-m2);
  float m12min=(m1+m2)*(m1+m2);
  double ang_ZMass = -cos12(  DmMass, m1, m2, m3, KstarMass, ZMass);
  double ang_ZMassWT= -cos12(  DmMass, m2, m1, m3, KstarMassWT, ZMass+kMass-piMass);
  cout<<Form("m12min=%f m12max=%f m13min=%f m13max=%f",m12min,m12max,m13min,m13max)<<endl;
  cout<<Form("Z4430 peak   -> cosk = %f", ang_ZMass)<<endl;
  cout<<Form("Z4430 peak wt-> cosk = %f", ang_ZMassWT)<<endl;
  TH2F *HxDalitz = new TH2F("HxDalitz","",150,m12min,m12max,150,m13min,m13max); 
  TH2D* Hxmmpikstarmass_SBL	   = new TH2D( "Hxmmpikstarmass_SBL"	   , "(K#pi mass)^{2} vs (#psi(2S)#pi mass)^{2} SB Left",         100,m12min,m12max,100,m13min,m13max );
  TH2D* Hxmmpikstarmass_SBR	   = new TH2D( "Hxmmpikstarmass_SBR"	   , "(K#pi mass)^{2} vs (#psi(2S)#pi mass)^{2} SB Right",         100,m12min,m12max,100,m13min,m13max );
  TH2D* Hxmmpikstarmass_Sign	   = new TH2D( "Hxmmpikstarmass_Sign"	   , "(K#pi mass)^{2} vs (#psi(2S)#pi mass)^{2} Sign",	  100,m12min,m12max,100,m13min,m13max );
  TH2D* Hxmmpikstarmass_Sign_b     = new TH2D( "Hxmmpikstarmass_Sign_b"    , "(K#pi mass)^{2} vs (#psi(2S)#pi mass)^{2} ",	  100,m12min,m12max,100,m13min,m13max );
  TH2D* Hxmmpikstarmass_Sign_c     = new TH2D( "Hxmmpikstarmass_Sign_c"    , "(K#pi mass)^{2} vs (#psi(2S)#pi mass)^{2} ",	  100,m12min,m12max,100,m13min,m13max );
  float scale1= 3;
  TH2D* HxmmpiCosK_Sign	          = new TH2D( "HxmmpiCosK_Sign"	          , "cos(#theta_{K}) vs #psi(2S)#pi mass (Signal, all events)",	       scale1*20,-1,1,scale1*40,sqrt(m13min),sqrt(m13max) );
  TH2D* HxmmpiCosK_Sign_Ks	  = new TH2D( "HxmmpiCosK_Sign_Ks"	  , "cos(#theta_{K}) vs #psi(2S)#pi mass (Signal, only K^{*0})" , scale1*20,-1,1,scale1*40,sqrt(m13min),sqrt(m13max) );
  TH2D* HxmmpiCosK_Sign_Md	  = new TH2D( "HxmmpiCosK_Sign_Md"	  , "cos(#theta_{K}) vs #psi(2S)#pi mass (Signal, 1<M(K#pi)^{2}<1.75)", scale1*20,-1,1,scale1*40,sqrt(m13min),sqrt(m13max) );
  TH2D* HxmmpiCosK_Sign_Tl	  = new TH2D( "HxmmpiCosK_Sign_Tl"	  , "cos(#theta_{K}) vs #psi(2S)#pi mass Sign (M^{2}>1.75)",   scale1*20,-1,1,scale1*40,sqrt(m13min),sqrt(m13max) );
  TH2D* HxmmpiCosK_Sign_NK	  = new TH2D( "HxmmpiCosK_Sign_NK"	  , "cos(#theta_{K}) vs #psi(2S)#pi mass Sign (Ks veto)" , scale1*20,-1,1,scale1*40,sqrt(m13min),sqrt(m13max) );
  TH2D* HkstarmassCosK_Sign= new TH2D( "HkstarmassCosK_Sign"	  , "kstarmass vs CosK_Sign",	       scale1*10,-1,1,scale1*60,m12min,m12max );
//  TH1D* Hxmmpi_Sign	  = new TH1D( "Hxmmpi_Sign"	 , "#psi(2S)#pi mass (Sign)"	 ,scale1*10, sqrt(m13min),sqrt(m13max));
  TH1D* Hxmmpi_SBL	  = new TH1D( "Hxmmpi_SBL"	 , "#psi(2S)#pi mass SB Left" 	 ,scale1*10, sqrt(m13min),sqrt(m13max));
  TH1D* Hxmmpi_SBR	  = new TH1D( "Hxmmpi_SBR"	 , "#psi(2S)#pi mass SB Right" 	 ,scale1*10, sqrt(m13min),sqrt(m13max));
  TH1D* Hxmmpiq_SBL	  = new TH1D( "Hxmmpiq_SBL"	 , "(#psi(2S)#pi mass)^{2} SB Left"	                 ,scale1*10, m13min,m13max);
  TH1D* Hxmmpiq_Sign	  = new TH1D( "Hxmmpiq_Sign"	 , "(#psi(2S)#pi mass)^{2} Sign"			 ,scale1*10, m13min,m13max);
  TH1D* Hxmmpiq_Sign_Ks	  = new TH1D( "Hxmmpiq_Sign_Ks"	 , "(#psi(2S)#pi mass)^{2} (Sign) (Kstar)"	 ,scale1*10, m13min,m13max);
  TH1D* Hxmmpiq_Sign_Md	  = new TH1D( "Hxmmpiq_Sign_Md"	 , "(#psi(2S)#pi mass)^{2} (Sign) (mass K#pi 1-1.75)" ,scale1*10, m13min,m13max);
  TH1D* Hxmmpiq_Sign_Tl	  = new TH1D( "Hxmmpiq_Sign_Tl"	 , "(#psi(2S)#pi mass)^{2} (Sign) (mass K#pi>1.75)"	 ,scale1*10, m13min,m13max);
  TH1D* Hxmmpiq_Sign_NK	  = new TH1D( "Hxmmpiq_Sign_NK"	 , "(#psi(2S)#pi mass)^{2} (Sign) (no Ks)"		 ,scale1*10, m13min,m13max);
  TH1D* Hxmmpi_Sign	 = new TH1D( "Hxmmpi_Sign"	, Form("#psi(2S)#pi mass (B^{0} signal: %3.2f<mass<%3.2f GeV/c^{2})",WMinSign,WMaxSign)  ,scale1*10, sqrt(m13min),sqrt(m13max));
  TH1D* Hxmmpi_Sign_Ks   = new TH1D( "Hxmmpi_Sign_Ks"	, "#psi(2S)#pi mass (Sign) (Kstar)"	     ,scale1*10, sqrt(m13min),sqrt(m13max));
  TH1D* Hxmmpi_Sign_Md   = new TH1D( "Hxmmpi_Sign_Md"	, "#psi(2S)#pi mass (Sign) (mass K#pi 1-1.75)"  ,scale1*10, sqrt(m13min),sqrt(m13max));
  TH1D* Hxmmpi_Sign_Tl   = new TH1D( "Hxmmpi_Sign_Tl"	, "#psi(2S)#pi mass (Sign) (mass K#pi>1.75)"    ,scale1*5, sqrt(m13min),sqrt(m13max));
  TH1D* Hxmmpi_Sign_NK   = new TH1D( "Hxmmpi_Sign_NK"	, "#psi(2S)#pi mass (Sign) (no Ks)"	     ,scale1*10, sqrt(m13min),sqrt(m13max));

  TH1D* Hxmmpiq_Sign_B0	  = new TH1D( "Hxmmpiq_Sign_B0"	 , "M(#psi(2S)#pi)^{2}       M(K#pi)<0.796 GeV/c^{2}" ,scale1*10, m13min,m13max);
  TH1D* Hxmmpiq_Sign_B1	  = new TH1D( "Hxmmpiq_Sign_B1"	 , "M(#psi(2S)#pi)^{2} 0.796<M(K#pi)<0.996 GeV/c^{2}" ,scale1*10, m13min,m13max);
  TH1D* Hxmmpiq_Sign_B2	  = new TH1D( "Hxmmpiq_Sign_B2"	 , "M(#psi(2S)#pi)^{2} 0.996<M(K#pi)<1.332 GeV/c^{2}" ,scale1*10, m13min,m13max);
  TH1D* Hxmmpiq_Sign_B3	  = new TH1D( "Hxmmpiq_Sign_B3"	 , "M(#psi(2S)#pi)^{2} 1.332<M(K#pi)<1.532 GeV/c^{2}" ,scale1*10, m13min,m13max);
  TH1D* Hxmmpiq_Sign_B4	  = new TH1D( "Hxmmpiq_Sign_B4"	 , "M(#psi(2S)#pi)^{2} 	   M(K#pi)>1.532 GeV/c^{2}" ,scale1*10, m13min,m13max);

  TH1D* Hxkstarmass2_Sign_B1  = new TH1D( "Hxkstarmass2_Sign_B1"     , "M(K#pi)^{2}    M(#psi(2S)#pi)^{2}<19 GeV^{2}/c^{4}"  , scale1*20,  m12min,m12max);
  TH1D* Hxkstarmass2_Sign_B2  = new TH1D( "Hxkstarmass2_Sign_B2"     , "M(K#pi)^{2} 19<M(#psi(2S)#pi)^{2}<20.5 GeV^{2}/c^{4}"  , scale1*20,  m12min,m12max);
  TH1D* Hxkstarmass2_Sign_B3  = new TH1D( "Hxkstarmass2_Sign_B3"     , "M(K#pi)^{2}    M(#psi(2S)#pi)^{2}>20.5 GeV^{2}/c^{4}"  , scale1*20,  m12min,m12max);

//   TH1D* Hxkstarmass2_SBL   = new TH1D( "Hxkstarmass2_SBL"      , "(K^{*0} mass)^{2} SBL"        , 100,  m12min,m12max);
//   TH1D* Hxkstarmass2_Sign  = new TH1D( "Hxkstarmass2_Sign"     , "(K^{*0} mass)^{2} Sign"       , 100,  m12min,m12max);
//   TH1D* Hxkstarmass2_Sign_b= new TH1D( "Hxkstarmass2_Sign_b"   , "(K^{*0} mass)^{2} Sign"       , 100,  m12min,m12max);
  TH1D* Hxkstarmass2_SBL   = new TH1D( "Hxkstarmass2_SBL"      , "(K#pi mass)^{2} SB Left"    , scale1*20,  m12min,m12max);
  TH1D* Hxkstarmass2_Sign  = new TH1D( "Hxkstarmass2_Sign"     , "(K#pi mass)^{2} Sign"       , scale1*20,  m12min,m12max);
  TH1D* Hxkstarmass2_Sign_b= new TH1D( "Hxkstarmass2_Sign_b"   , "(K#pi mass)^{2} Sign"       , scale1*20,  m12min,m12max);
//
  TF2  *Dalitz2 = new TF2("Dalitz2",DalitzModel2,m12min,m12max,m13min,m13max,8);
  Double_t parDali2[8];
  parDali2[0]=DmMass;
  parDali2[1]=m1;
  parDali2[2]=m2;
  parDali2[3]=m3;
  parDali2[4]=KstarMass;
  parDali2[5]=KstarGamm;
  parDali2[6]=ZMass;
  parDali2[7]=ZGamm;
  Dalitz2->SetNpy(1000);
  Dalitz2->SetNpx(1000);
  Dalitz2->SetParameters(parDali2);
  TH2* HDalitz2 = (TH2*)Dalitz2->CreateHistogram();
  HDalitz2->Scale(1);
  TH1* HDalitz2Prox =  HDalitz2->ProjectionX("HDalitz2Prox");
  TH1* HDalitz2Proy =  HDalitz2->ProjectionY("HDalitz2Proy");
//
  TF2  *Dalitz = new TF2("Dalitz",DalitzModel,m12min,m12max,m13min,m13max,6);
  Double_t parDali[6];
  parDali[0]=DmMass;
  parDali[1]=m1;
  parDali[2]=m2;
  parDali[3]=m3;
  parDali[4]=KstarMass;
  parDali[5]=KstarGamm;
  Dalitz->SetNpy(1000);
  Dalitz->SetNpx(1000);
  Dalitz->SetParameters(parDali);
  TH2* HDalitz = (TH2*)Dalitz->CreateHistogram();
  TH1* HDalitzProx =  HDalitz->ProjectionX("HDalitzProx");
  TH1* HDalitzProy =  HDalitz->ProjectionY("HDalitzProy");
//
//
  TF2  *DKinematic = new TF2("DKinematic",KineLimit,m12min,m12max,m13min,m13max,4);
  
//  TF2  *DKinematic = new TF2("DKinematic",KineLimit,0.,28.,0.,28.,4);
  TF2  *DKinematic2 = new TF2("DKinematic2",KineLimit,0.,28.,0.,28.,4);
  Double_t Dalipar[4];
  Dalipar[0]=DmMass;
  Dalipar[1]=m1;
  Dalipar[2]=m2;
  Dalipar[3]=m3;
  DKinematic->SetNpy(4000);
  DKinematic->SetNpx(1000);
  DKinematic->SetParameters(Dalipar);
  Dalipar[0]=DmMass;
  Dalipar[1]=m1;
  Dalipar[2]=m2+piMass;
  Dalipar[3]=m3+piMass;
  DKinematic2->SetNpy(1000);
  DKinematic2->SetNpx(1000);
  DKinematic2->SetParameters(Dalipar);
//   TF2  *FCutkst2 = new TF2("FCutkst2 ",FunCutkst2,  -0.4,0.6,  -0.4,0.6,4);
//   Double_t kst2par[4];
//   kst2par[0]=x0Cut;
//   kst2par[1]=x1Cut;
//   kst2par[2]=y0Cut;
//   kst2par[3]=y1Cut;
//   FCutkst2->SetNpx(2000);
//   FCutkst2->SetNpy(2000);
//   FCutkst2->SetParameters(kst2par);
//   TF2  *FCutmmpi2 =  new TF2("FCutmmpi2",FunCutmmpi2,3.,6.,3.5,6.,8);
//   Double_t mmpi2par[8];
//   mmpi2par[0]=CutX1;
//   mmpi2par[1]=CutX2;
//   mmpi2par[2]=CutY1;
//   mmpi2par[3]=CutY2;
//   mmpi2par[4]=x_0Cut;
//   mmpi2par[5]=x_1Cut;
//   mmpi2par[6]=y_0Cut;
//   mmpi2par[7]=y_1Cut;
//   FCutmmpi2->SetNpx(2000);
//   FCutmmpi2->SetNpy(2000);
//   FCutmmpi2->SetParameters(mmpi2par);
  
//  TF2  *DKinematic = new TF2("DKinematic",KineLimit, m12min,m12max,m13min,m13max,1);
//
//  double y1Cut2=-0.25;
//nentries=300000;
//double xmax=-999;
  hely_zed = -cos12(  B0Mass, piMass, kMass, MassPsi2S, KstarMass,ZMass);
  cout<<Form("cosk(Z4430) = %f",hely_zed)<<endl;
  hely_zmd = -cos12(  B0Mass, piMass, kMass, MassPsi2S,sqrt(1.375),ZMass);
  cout<<Form("cosk(sqrt(1.375)) = %f (middle position in Dalitz region > K* on Dalitz)",hely_zmd)<<endl;
   
  bool firstmsg = true; 
  if (NK ){
   cout<<"Warning: set passB0Psi_*=true \n"<<endl;
   passB0Psi_jpsi = true;
   passB0Psi_psip = true;
   passB0Psi_lmnr = true;
  }else{
   dnn_prob = 9999;
  }
  for (Int_t i=0;i<nentries;i++) {
         RecoB0TreeOut->GetEntry(i);
//	 
        if( (WrongTag&&MC) && (tagB0 == 1 && genSignal == 1) || ( tagB0 == 0 && genSignal == 2)) continue;  //(cut correctly tagged)
        if( (CorreTag&&MC) && (tagB0 == 1 && genSignal == 2) || ( tagB0 == 0 && genSignal == 1)) continue;  //(cut wrongly   tagged)
        if ( dnn_prob<dnn_cut) continue;
//	  if(!passB0Psi_psip) continue;
//	 if(!passB0Psi_jpsi) continue;
         if(!passB0Psi_psip) continue;
//         if(!passB0Psi_lmnr) continue;
         if (!NK && (tagged_mass>XMaxSign || tagged_mass<XMinSign) ) continue;
         ROOT::Math::PtEtaPhiMVector B0Vec (bPt,bEta,bPhi,tagged_mass);
         ROOT::Math::PtEtaPhiMVector KstarVec (kstPt,kstEta,kstPhi,kstarmass);
//         ROOT::Math::PtEtaPhiMVector KstarVec (kstPt,kstEta,kstPhi,kstarmass);
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

	 TLorentzVector b0_lv_boosted (B0Vec.Px(),B0Vec.Py(),B0Vec.Pz(),B0Vec.E());
	 TLorentzVector KstarLorentz (KstarVec.Px(),KstarVec.Py(),KstarVec.Pz(),KstarVec.E());
	 TLorentzVector *kaon_lv_boosted =0;
	 ROOT::Math::PtEtaPhiMVector  *KappaVec=0;
	 
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
	 
// 	 if(NK){
// 	  if (pipPt>pimPt){
// 	    tagB0=1;
// 	  }else{
// 	    tagB0=0;
// 	  }
// 	 } 
	 
         if (tagB0>0) { //B0
	  if((KpVec.Pt()>PimVec.Pt())) leadK =true;
//	  lambdab1 =  (JpsiVec+KmVec+PpVec).mass();
//	  lambdab2 =  (JpsiVec+KpVec+PmVec).mass();  
          
          kaon_lv_boosted = new TLorentzVector(KpVec.Px(),KpVec.Py(),KpVec.Pz(),KpVec.E());
          Kstarmup = (MupVec+PimVec+KpVec).mass();
          Kstarmum = (MumVec+PimVec+KpVec).mass();
	  KappaVec = new ROOT::Math::PtEtaPhiMVector(kstTrkpPt,kstTrkpEta,kstTrkpPhi,kMass);
//	  boostKst = new BoostVector(Kstarmup);
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
	  tagged_massAll=bMass;
	  kstarmassInvAll=kstBarMass;
	  tagged_massInvAll=bBarMass;
	 }else{ //B0bar
	  if((KmVec.Pt()>PipVec.Pt())) leadK =true;
          kaon_lv_boosted = new TLorentzVector(KmVec.Px(),KmVec.Py(),KmVec.Pz(),KmVec.E());
          Kstarmup = (MupVec+PipVec+KmVec).mass();
          Kstarmum = (MumVec+PipVec+KmVec).mass();
	  KappaVec = new ROOT::Math::PtEtaPhiMVector(kstTrkmPt,kstTrkmEta,kstTrkmPhi,kMass);
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
	  tagged_massAll=bBarMass;
	  kstarmassInvAll=kstMass;
	  tagged_massInvAll=bMass;
	 }
	 if(PmVec.Pt()>PpVec.Pt()){
	  lambdab  =(JpsiVec+KpVec+PmVec).mass();
	  nolambdab=(JpsiVec+KmVec+PpVec).mass();
	 }else{
	  lambdab  =(JpsiVec+KmVec+PpVec).mass();
	  nolambdab=(JpsiVec+KpVec+PmVec).mass();
	 } 
 	 if (NK){
	   if(INV){
//
  	    kstarmass=kstarmassInvAll ;
	    tagged_mass=tagged_massInvAll;
//
 		 if( firstmsg ){
		  cout<<"Warning: set kstarmass=kstarmassInvAll \n"<<endl;
 		  cout<<"Warning: set tagged_mass=tagged_massInvAll \n"<<endl;
		 }
	   }else{

	    kstarmass=kstarmassAll ;
	    tagged_mass=tagged_massAll;
//	    
 		 if( firstmsg ){
		  cout<<"Warning: set kstarmass=kstarmassAll \n"<<endl;
 		  cout<<"Warning: set tagged_mass=tagged_massAll \n"<<endl;
		 }
	   } 
	   firstmsg=false;
 	 }
	 if(INV){
          mmpix=mmpi1;
 	  kstarmassx=kstarmass;
	 }else{ 
 	  mmpix=mmpi2;
	  kstarmassx=kstarmass;
	 }  
         int KLimit1 = BackgroundPlot( kst2*kst2            , mmpi1*mmpi1,  mmpiKaon , piMass, kMass, mumuMass);
         int KLimit0 = BackgroundPlot( kstarmassx*kstarmassx, mmpix*mmpix,  B0Mass , piMass, kMass, MassPsi2S);
	 
//	 if(leadK) continue;
//         int KLimit0 = BackgroundPlot( kstarmassx*kstarmassx, mmpix*mmpix,  tagged_mass , piMass, kMass, mumuMass);
//	 if(KLimit0==0 && NK) continue;
// 	 if( kkMass<1.035&& NK)  continue;
// 	 if( fabs(masspipi-RhoMass)<RhoGam && NK)  continue;
// 	 if( bLBS/bLBSE<20&& NK)  continue;
// 	 if( bCosAlphaBS<0.9 )   continue;
// 	 if( bVtxCL<0.2 )   continue;
// 	 if( fabs(kstTrkmEta)>1.2) continue;
// 	 if( fabs(kstTrkpEta)>1.2) continue;
// 	 if( fabs(mumuMass-MassJPsi)>3*mumuMassE) continue;
// 	 if( fabs(mumuMass-MassPsi2S)>3*mumuMassE) continue;
// 	 if( kstTrkmDCABS/kstTrkmDCABSE<1.5) continue;
// 	 if( kstTrkpDCABS/kstTrkpDCABSE<1.5) continue;
//  	 if(kstTrkmMinIP2D>0.002)  continue;
//  	 if(kstTrkpMinIP2D>0.002)  continue;
// 	 if( bLBS/bLBSE>10&& NK)  continue;
// 	 
//              if(kstarmassx*kstarmassx<1.0) continue;
//              if(kstarmassx*kstarmassx>1.9) continue;
//	 if(kstarmassx*kstarmassx>0.8) continue;
//         if((fabs(kstarmassx-KstarMass)>0.1)) continue;
//	 if(!(kstarmassx*kstarmassx>0.7&&kstarmassx*kstarmassx<1.`)) continue;
//	 if((kstarmassx*kstarmassx>0.65&&kstarmassx*kstarmassx<0.95)) continue;
	 
 	 TVector3 boostKst =  KstarLorentz.BoostVector();
	 b0_lv_boosted.Boost(-boostKst);
	 kaon_lv_boosted->Boost(-boostKst);
         double cos_theta_k_r = computeCosine(-b0_lv_boosted.Px(),
 					      -b0_lv_boosted.Py(),
 					      -b0_lv_boosted.Pz(),
 					       kaon_lv_boosted->Px(),
 					       kaon_lv_boosted->Py(),
 					       kaon_lv_boosted->Pz()
 					     );
/////////////////	if (NK) cos_theta_k =cos_theta_k_r;
   	  double xDmMass= tagged_mass;
//   	  double xDmMass= XMaxSBL;
//  	  double xm3    = MassPsi2S;
  	  double xm3    = mumuMass;
 	  double xm2    = kMass;
 	  double xm1    = piMass;
          double xm23max=(xDmMass-xm1)*(xDmMass-xm1);
          double xm23min=(xm2+xm3)*(xm2+xm3);;
          double xm12max=(xDmMass-xm3)*(xDmMass-xm3);
          double xm12min=(xm1+xm2)*(xm1+xm2);
          double xm13max=(xDmMass-xm2)*(xDmMass-xm2);
//	  double xm13max=(b0_lv_boosted*b0_lv_boosted-2*(*kaon_lv_boosted)*b0_lv_boosted+(*kaon_lv_boosted)*(*kaon_lv_boosted));
	  double xm13min=(xm1+xm3)*(xm1+xm3);
//	  double xm13max=21.68;
//          double xm13max=B0Mass*B0Mass+xm1*xm1+xm2*xm2+MassPsi2S*MassPsi2S-KstarMass*KstarMass-xm23min;
// 	  double xm13max=xDmMass*xDmMass+xm1*xm1+xm2*xm2+(MassPsi2S-piMass)*(MassPsi2S-piMass)-xm12min-xm23min;
//          double delta=(tagged_mass*tagged_mass+m1*m1+m2*m2+m3*m3)-(mmpi2*mmpi2+kstarmass*kstarmass+mmka1*mmka1);
	  double mm2max=(mmpi2-piMass)*(mmpi2-piMass);
	  double mm2min=4*muonMass_;
//	  cout<<delta<<endl;
//	  hely_ang = (xm13max+xm13min-2*mmpi2*mmpi2)/(xm13max-xm13min);
//	  double cosxy=cos12(  xDmMass, xm1, xm2, xm3, kstarmass, mmpi2);
	  HDalitz2->GetRandom2(kstarmass_rnd,mmpi2_rnd);
// 	  cout<<Form("kstarmass_rnd=%f mmpi2_rnd=%f",kstarmass_rnd,mmpi2_rnd)<<endl;
// 	  cout<<Form("mmka1=%f mmpi2=%f",kstarmass,mmpi2)<<endl;
	  hely_rnd = -cos12(  xDmMass, xm1, xm2, xm3, sqrt(kstarmass_rnd),sqrt(mmpi2_rnd));
	  hely_wnd = -cos12(  xDmMass, xm2, xm1, xm3, sqrt(kstarmass_rnd)+KstarMassWT-KstarMass,sqrt(mmpi2_rnd)+kMass-piMass);
	  
	  
	  

	  hely_ang = -cos12(  xDmMass, xm1, xm2, xm3, kstarmassx, mmpix);
//	  hely_ang = vettore(xDmMass,kstarmass,mmpi2,xm1,xm2,xm3,KstarMass,KstarGamm);
	  //
//          double delta=(xDmMass*xDmMass+xm1*xm1+xm2*xm2+xm3*xm3)-(xm23max+xm12max);
//	  hely_ang = (m13max+m13min-2*((tagged_mass*tagged_mass+m1*m1+m2*m2+m3*m3)-(kstarmass*kstarmass+mmka1*mmka1)))/(m13max-m13min);
	  
//	  double m13max_cosk = (m13min*(cos_theta_k+1)-2*mmpi2*mmpi2)/(cos_theta_k-1);
	  
//	  cout<<m13max_cosk-m13max <<endl;
//	  hely_ang = (m13max+m13min-2*mmpi2*mmpi2)/( ((tagged_mass-m2)+(m1+mumuMass))*((tagged_mass-m2)-(m1+mumuMass)));
          
	 
//          kaon_lv_boosted.Boost(-boostKst)
	 
//	 kstarmass=kstarmassAll;
// 	 if(cos_theta_k>1||cos_theta_k<-1) continue;
// 	 if(cos_theta_l>1||cos_theta_l<-1) continue;
// 	 if(phi_kst_mumu<-TMath::Pi()||phi_kst_mumu>TMath::Pi()) continue;
// 	 kstarmassx=mmka1;
// 	 kstarmassx=masspipi;
//          mmpix=Kstarmup;
// 	 kstarmassx=Kstarmum;
	 
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
	 
//	 HxCosK_Sign->Fill(cos_theta_k,-1);
  	  if(XCutCosK) HxHely_Left->Fill(hely_ang);
  	  if(XCutCosK) HxCosK_Left_H->Fill(cos_theta_k);
 	  if(XCutCosK) HxCosK_Left_R->Fill(cos_theta_k_r);
          Hxmmpi_SBL->Fill(mmpi2);
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
          Hxmmpi_SBR->Fill(mmpix);
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
         if(tagged_mass>WMinSign && tagged_mass<WMaxSign ){
// 	 
// Dalitz slice to be comape with Belle	results
//
          double mmpixq = mmpix*mmpix;
          if( mmpixq<19 	       )  Hxkstarmass2_Sign_B1 ->Fill(kstarmassx*kstarmassx);
          if( mmpixq>19 && mmpixq<20.5 )  Hxkstarmass2_Sign_B2 ->Fill(kstarmassx*kstarmassx);
          if(		   mmpixq>20.5 )  Hxkstarmass2_Sign_B3 ->Fill(kstarmassx*kstarmassx);
	  if( kstarmassx<0.796                     ) Hxmmpiq_Sign_B0->Fill(mmpix*mmpix);
	  if( kstarmassx>0.796 && kstarmassx<0.996 ) Hxmmpiq_Sign_B1->Fill(mmpix*mmpix);
	  if( kstarmassx>0.996 && kstarmassx<1.332 ) Hxmmpiq_Sign_B2->Fill(mmpix*mmpix);
	  if( kstarmassx>1.332 && kstarmassx<1.532 ) Hxmmpiq_Sign_B3->Fill(mmpix*mmpix);
	  if(                     kstarmassx>1.532 ) Hxmmpiq_Sign_B4->Fill(mmpix*mmpix);
//	  if(xmax<xm13max) xmax=xm13max; 
//	   double coeff = (B0Mass*B0Mass-kstarmass*kstarmass-kMass*kMass)/(2*B0Mass*kstarmass);
	  if(XCutCosK) HxHely_Sign->Fill(hely_ang);
  	  if(XCutCosK) HxCosK_Sign_H->Fill(cos_theta_k);
// 	  if(MC){
// 	    if(fabs(ZMass-mmpix)<0.1&&!WrongTag){
//  	     if(XCutCosK) HxCosK_Sign_R->Fill(cos_theta_k);
// 	    } 
// 	    if(fabs(ZMass-mmpi1)<0.1&&WrongTag){
//  	     if(XCutCosK) HxCosK_Sign_W->Fill(cos_theta_k);
// 	    }
// 	  }else{
//	    if(fabs(ZMass-mmpi2)<0.1||fabs(ZMass-mmpi1)<0.1){
 	     if(XCutCosK) HxCosK_Sign_R->Fill(hely_rnd);
 	     if(XCutCosK) HxCosK_Sign_W->Fill(hely_wnd);
// 	    if(fabs(ZMass-mmpi2)<0.03){
//  	     if(XCutCosK) HxCosK_Sign_R->Fill(cos_theta_k);
// 	    } 
// 	    if(fabs(ZMass-mmpi1)<0.03){
//  	     if(XCutCosK) HxCosK_Sign_W->Fill(cos_theta_k);
// 	    } 
//	  }   
// 	  if(XCutCosK) HxCosK_Sign_R->Fill(cos_theta_k_r);
// 	  if(XCutCosK) HxCosK_Sign_R->Fill(cos_theta_k_r,1+0.4*cos_theta_k_r*cos_theta_k_r);
//  	  HxMassCosK_Sign->Fill(tagged_mass,cos_theta_k);
  	  if(XCutCosK) HxCosK_Sign->Fill(cos_theta_k);
	  if(fabs(kstarmassx-KstarMass)<0.1) HxCosK_Sign_Ks->Fill(cos_theta_k);
	  if(fabs(kstarmassx-KstarMass)>0.1) HxCosK_Sign_NK->Fill(cos_theta_k);
	  if(kstarmassx*kstarmassx>1.1&&kstarmassx*kstarmassx<1.75) HxCosK_Sign_Md->Fill(cos_theta_k);
	  if(kstarmassx*kstarmassx>1.75) HxCosK_Sign_Tl->Fill(cos_theta_k);
	  if(XCutCosK&&!leadK) HxCosK_Sign_Lpi->Fill(cos_theta_k);
	  if(XCutCosK&& leadK) HxCosK_Sign_LK ->Fill(cos_theta_k);
	  if(XCutCosK) HxBsMass_Sign->Fill(mmpiKaon);
	  if(!XCutCosK) HxBsMass_SignC->Fill(mmpiKaon);
	  HxMass_Sign->Fill(tagged_mass);
 	  if(XCutCosK) HxLambdab_Sign->Fill(lambdab);
 	  if(!XCutCosK) HxLambdab_Sign_c->Fill(lambdab);
// 	  if(XCutCosK) HxLambdab_Sign->Fill(lambdab1);
// 	  if(XCutCosK) HxLambdab_Sign->Fill(lambdab2);
//          Hxmmpi_Sign->Fill(mmpi1);

//          Hxmmpikstarmass_Sign->Fill(kstarmass*kstarmass,mmpix*mmpix);
          HxmmpiCosK_Sign->Fill(cos_theta_k,mmpix);
          if(fabs(kstarmassx-KstarMass)<0.1) HxmmpiCosK_Sign_Ks->Fill(cos_theta_k,mmpix);
          if(fabs(kstarmassx-KstarMass)>0.1) HxmmpiCosK_Sign_NK->Fill(cos_theta_k,mmpix);
          if(kstarmassx*kstarmassx>1.&&kstarmassx*kstarmassx<1.75) HxmmpiCosK_Sign_Md->Fill(cos_theta_k,mmpix);
          if(kstarmassx*kstarmassx>1.75) HxmmpiCosK_Sign_Tl->Fill(cos_theta_k,mmpix);
          if(fabs(kstarmassx-KstarMass)>0.1) HkstarmassCosK_Sign->Fill(cos_theta_k,kstarmassx*kstarmassx);

//          if(mmpix*mmpix>22) Hxkstarmass2_Sign->Fill(kstarmass*kstarmass);
          Hxmmpikstarmass_Sign->Fill(kstarmassx*kstarmassx,mmpix*mmpix);
          Hxkstarmass2_Sign->Fill(kstarmassx*kstarmassx);
//
          Hxmmpi_Sign->Fill(mmpix);
          if(fabs(kstarmassx-KstarMass)<0.1) Hxmmpi_Sign_Ks->Fill(mmpix);
          if(fabs(kstarmassx-KstarMass)>0.1) Hxmmpi_Sign_NK->Fill(mmpix);
          if(kstarmassx*kstarmassx>1.&&kstarmassx*kstarmassx<1.75) Hxmmpi_Sign_Md->Fill(mmpix);
          if(kstarmassx*kstarmassx>1.75) Hxmmpi_Sign_Tl->Fill(mmpix);
//
          Hxmmpiq_Sign->Fill(mmpix*mmpix);
          if(fabs(kstarmassx-KstarMass)<0.1) Hxmmpiq_Sign_Ks->Fill(mmpix*mmpix);
          if(fabs(kstarmassx-KstarMass)>0.1) Hxmmpiq_Sign_NK->Fill(mmpix*mmpix);
          if(kstarmassx*kstarmassx>1.&&kstarmassx*kstarmassx<1.75) Hxmmpiq_Sign_Md->Fill(mmpix*mmpix);
          if(kstarmassx*kstarmassx>1.75) Hxmmpiq_Sign_Tl->Fill(mmpix*mmpix);
//
	  if(cos_theta_k>-0.6&&cos_theta_k<-0.2){
           Hxmmpikstarmass_Sign_b->Fill(kstarmassx*kstarmassx,mmpix*mmpix);
//           if(!leadK) Hxmmpi_Sign->Fill(mmpi1);
           if( leadK) Hxmmpikstarmass_Sign_LK->Fill(mmpix*mmpix,kstarmassx*kstarmassx);
           if(!leadK) Hxmmpikstarmass_Sign_Lpi->Fill(mmpix*mmpix,kstarmassx*kstarmassx);
           Hxpsi2s_Sign->Fill(mumuMass);
//	   HxLambdab_Sign_c->Fill(lambdab);
	  }else{ 
           if( leadK) Hxmmpikstarmass_Sign_LK_c ->Fill(mmpix*mmpix,kstarmassx*kstarmassx);
           if(!leadK) Hxmmpikstarmass_Sign_Lpi_c->Fill(mmpix*mmpix,kstarmassx*kstarmassx);
//	   HxLambdab_Sign->Fill(lambdab);
           Hxpsi2s_Sign_c->Fill(mumuMass);
	  } 
	 } 
  }
  
//  cout<<"xmax="<<xmax<<endl;
  gStyle->SetPadBorderMode(0);
  gStyle->SetOptStat("em");
//  gStyle->SetOptStat(000000);
  gStyle->SetOptFit(000000);
  HxCosK->SetMinimum(0);
  HxCosK_Sign->SetMinimum(0);
  HxCosK_SBL->SetMinimum(0);
  HxCosK_SBR->SetMinimum(0);
  
//=======================================  
  gStyle->SetOptStat(000000);
  cStudies->cd(1);
  HxCosK_LK->SetLineColor(kRed);
  HxCosK_Lpi->SetLineColor(kBlue);
  HxCosK->GetXaxis()->SetTitle("cos(#theta_{K})");
  HxCosK->Draw();
  TLegend* leg_CosK = new TLegend(0.52,0.72,0.70,0.85);
  leg_CosK->SetTextSize(0.025) ;
  leg_CosK->SetTextAlign(13);
  leg_CosK->SetBorderSize(0.);
  leg_CosK->SetFillStyle(0);
  leg_CosK->AddEntry(HxCosK ,"cos(#theta_{K})","l");
  leg_CosK->AddEntry(HxCosK_LK ,"cos(#theta_{K})  K leading track","l");
  leg_CosK->AddEntry(HxCosK_Lpi ,"cos(#theta_{K}) #pi leading track","l");
  leg_CosK->Draw("same");
  HxCosK_LK->Draw("same");
  HxCosK_Lpi->Draw("same");
  cStudies->cd(2);
//  HxMassCosK_SBL->Draw("contz");
  HxCosK_SBL_LK->SetLineColor(kRed);
  HxCosK_SBL_Lpi->SetLineColor(kBlue);
  HxCosK_SBL->GetXaxis()->SetTitle("cos(#theta_{K})");
  HxCosK_SBL->Draw();
  HxCosK_SBL_LK->Draw("same");
  HxCosK_SBL_Lpi->Draw("same");
  cStudies->cd(3);
//  HxMassCosK_Sign->Draw("contz");
  HxCosK_Sign_LK->SetLineColor(kRed);
  HxCosK_Sign_Lpi->SetLineColor(kBlue);
  HxCosK_Sign->GetXaxis()->SetTitle("cos(#theta_{K})");
  HxCosK_Sign->Draw();
  HxCosK_Sign_LK->Draw("same");
  HxCosK_Sign_Lpi->Draw("same");
  TArrow *arrowH=new TArrow(hely_zed,HxCosK_Sign->GetMinimum(),hely_zed,0.9*HxCosK_Sign->GetMaximum());
  TPaveLabel *labelH= new TPaveLabel(hely_zed,0.92*HxCosK_Sign->GetMaximum(),hely_zed+0.7,0.99*HxCosK_Sign->GetMaximum(),"Z(4430)");
  arrowH->SetLineColor(kRed);
  arrowH->SetArrowSize(0.01);
  arrowH->Draw("same,<");
  labelH->Draw("same");
  cStudies->cd(4);
//  HxMassCosK_SBR->Draw("contz");
  HxCosK_SBR_LK->SetLineColor(kRed);
  HxCosK_SBR_Lpi->SetLineColor(kBlue);
  HxCosK_SBR->GetXaxis()->SetTitle("cos(#theta_{K})");
  HxCosK_SBR->Draw();
  HxCosK_SBR_LK->Draw("same");
  HxCosK_SBR_Lpi->Draw("same");
  cStudies->cd(5);
//  Hxmmpikstarmass_SBL->Draw("COLZ");
  Hxmmpikstarmass_Sign->Draw("COLZ");
  Hxmmpikstarmass_Sign->GetYaxis()->SetTitleOffset(1.5);
  Hxmmpikstarmass_Sign->GetXaxis()->SetTitle("M(K#pi)^{2} GeV^{2}/c^{4}");
  Hxmmpikstarmass_Sign->GetYaxis()->SetTitle("M(#psi(2S)#pi)^{2} GeV^{2}/c^{4}");
//  Hxmmpikstarmass_Sign_b->SetMarkerColor(kRed);
//  Hxmmpikstarmass_Sign_b->Draw("same");
  DKinematic->Draw("same");
//  Hxmmpikstarmass_Sign_LK->Draw();
//   HxMass_SBL->Draw();
//   HxBsMass_SBL->Draw("same");
  cStudies->cd(6);
//  Hxmmpikstarmass_Sign_LK_c->Draw();
  Hxmmpiq_Sign->Draw();
  Hxmmpiq_SBL->SetLineColor(kRed);
  double NFact=0.3;
//   if(!NK){
//    HDalitzProy ->Scale(NFact*Hxmmpiq_Sign->Integral()*Hxmmpiq_Sign->GetXaxis()->GetBinWidth(1)/HDalitzProy->Integral()/HDalitzProy->GetXaxis()->GetBinWidth(1));
//    HDalitz2Proy->Scale(NFact*Hxmmpiq_Sign->Integral()*Hxmmpiq_Sign->GetXaxis()->GetBinWidth(1)/HDalitz2Proy->Integral()/HDalitz2Proy->GetXaxis()->GetBinWidth(1));
//    HDalitzProy ->SetLineColor(kMagenta);
//    HDalitz2Proy->SetLineColor(kBlue);
//    HDalitzProy ->Draw("HIST SAME C");
//    HDalitz2Proy->Draw("HIST SAME C");
//   }
  Hxmmpiq_Sign->GetXaxis()->SetTitle("M(#psi(2S)#pi)^2 GeV^{2}/c^{4}");
  Hxmmpiq_Sign->Draw("same");
  Hxmmpiq_SBL->Draw("same");
  TArrow *arrowZ2=new TArrow(ZMass*ZMass,Hxmmpiq_Sign->GetMinimum(),ZMass*ZMass,0.99*Hxmmpiq_Sign->GetMaximum());
  TPaveLabel *labelZ2= new TPaveLabel(ZMass*ZMass,0.72*Hxmmpiq_Sign->GetMaximum(),ZMass*ZMass+2,0.62*Hxmmpiq_Sign->GetMaximum(),"Z(4430)");
  arrowZ2->SetLineColor(kRed);
  arrowZ2->SetArrowSize(0.01);
  arrowZ2->Draw("same,<");
  labelZ2->Draw("same");
//   HxMass_Sign->Draw();
//   HxBsMass_Sign->Draw("same");
//   HxBsMass_SignC->Draw("same");
  cStudies->cd(7);
////   HxMass_SBR->Draw();
////   HxBsMass_SBR->Draw("same");
////  Hxmmpikstarmass_Sign_Lpi->Draw();

  Hxkstarmass2_Sign->GetXaxis()->SetTitle("M(K#pi)^{2} GeV^{2}/c^{4}");
  Hxkstarmass2_Sign->Draw();
  Hxkstarmass2_SBL->SetLineColor(kRed);

//   if(!NK){
//    HDalitzProx ->SetLineColor(kMagenta);
//    HDalitz2Prox->SetLineColor(kBlue);
//    HDalitzProx ->Scale(NFact*Hxkstarmass2_Sign->Integral()*Hxkstarmass2_Sign->GetXaxis()->GetBinWidth(1)/HDalitzProx->Integral()/HDalitzProx->GetXaxis()->GetBinWidth(1));
//    HDalitz2Prox->Scale(NFact*Hxkstarmass2_Sign->Integral()*Hxkstarmass2_Sign->GetXaxis()->GetBinWidth(1)/HDalitz2Prox->Integral()/HDalitz2Prox->GetXaxis()->GetBinWidth(1));
//    HDalitzProx ->Draw("HIST SAME C");
//    HDalitz2Prox->Draw("HIST SAME C");
//   }
  Hxkstarmass2_Sign->Draw("same");
  Hxkstarmass2_SBL->Draw("same");
  cStudies->cd(8);
//   Hxmmpi_SBR->SetLineColor(kRed);
//   Hxmmpi_SBL->Draw();
//   Hxmmpi_SBR->Draw("same");
  
   Hxmmpikstarmass_SBL->Draw("COLZ");
   Hxmmpikstarmass_SBL->GetYaxis()->SetTitleOffset(1.5);
   Hxmmpikstarmass_SBL->GetXaxis()->SetTitle("(M(K#pi)^{2} GeV^{2}/c^{4}");
   Hxmmpikstarmass_SBL->GetYaxis()->SetTitle("(M(#psi(2S)#pi)^{2} GeV^{2}/c^{4}");
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
  sss<<PDFName<<"-studies"<<TXTMASSW<<TYPEFILE;
  gSystem->Exec(Form("mv %s %s.tmp",sss.str().c_str(),sss.str().c_str()));
  cStudies->Print(sss.str().c_str());
  sss.clear();
  sss.str("");
//=============================================================
  cBelle->cd(1);
  gPad->SetLogy();
  Hxkstarmass2_Sign_B1->GetYaxis()->SetTitleOffset(1.5);
  Hxkstarmass2_Sign_B1->GetYaxis()->SetTitle(Form("Events/(%3.2f GeV/c^{2})",(m12max-m12min)/(scale1*20)));
  Hxkstarmass2_Sign_B1->GetXaxis()->SetTitle("M(K#pi)^{2} GeV^{2}/c^{4}");
  Hxkstarmass2_Sign_B1->Draw();
  cBelle->cd(2);
  gPad->SetLogy();
  Hxkstarmass2_Sign_B2->GetYaxis()->SetTitleOffset(1.5);
  Hxkstarmass2_Sign_B2->GetYaxis()->SetTitle(Form("Events/(%3.2f GeV/c^{2})",(m12max-m12min)/(scale1*20)));
  Hxkstarmass2_Sign_B2->GetXaxis()->SetTitle("M(K#pi)^{2} GeV^{2}/c^{4}");
  Hxkstarmass2_Sign_B2->Draw();
  cBelle->cd(3);
  gPad->SetLogy();
  Hxkstarmass2_Sign_B3->GetYaxis()->SetTitleOffset(1.5);
  Hxkstarmass2_Sign_B3->GetYaxis()->SetTitle(Form("Events/(%3.2f GeV/c^{2})",(m12max-m12min)/(scale1*20)));
  Hxkstarmass2_Sign_B3->GetXaxis()->SetTitle("M(K#pi)^{2} GeV^{2}/c^{4}");
  Hxkstarmass2_Sign_B3->Draw();
  cBelle->cd(4);
  Hxmmpiq_Sign_B0->GetYaxis()->SetTitleOffset(1.5);
  Hxmmpiq_Sign_B0->GetYaxis()->SetTitle(Form("Events/(%3.2f GeV/c^{2})",(m13max-m13min)/(scale1*10)));
  Hxmmpiq_Sign_B0->GetXaxis()->SetTitle("M(#psi(2S)#pi)^{2} GeV^{2}/c^{4}");
  Hxmmpiq_Sign_B0->Draw();
  cBelle->cd(5);
  Hxmmpiq_Sign_B1->GetYaxis()->SetTitleOffset(1.5);
  Hxmmpiq_Sign_B1->GetYaxis()->SetTitle(Form("Events/(%3.2f GeV/c^{2})",(m13max-m13min)/(scale1*10)));
  Hxmmpiq_Sign_B1->GetXaxis()->SetTitle("M(#psi(2S)#pi)^{2} GeV^{2}/c^{4}");
  Hxmmpiq_Sign_B1->Draw();
  cBelle->cd(6);
  Hxmmpiq_Sign_B2->GetYaxis()->SetTitleOffset(1.5);
  Hxmmpiq_Sign_B2->GetYaxis()->SetTitle(Form("Events/(%3.2f GeV/c^{2})",(m13max-m13min)/(scale1*10)));
  Hxmmpiq_Sign_B2->GetXaxis()->SetTitle("M(#psi(2S)#pi)^{2} GeV^{2}/c^{4}");
  Hxmmpiq_Sign_B2->Draw();
  cBelle->cd(7);
  Hxmmpiq_Sign_B3->GetYaxis()->SetTitleOffset(1.5);
  Hxmmpiq_Sign_B3->GetYaxis()->SetTitle(Form("Events/(%3.2f GeV/c^{2})",(m13max-m13min)/(scale1*10)));
  Hxmmpiq_Sign_B3->GetXaxis()->SetTitle("M(#psi(2S)#pi)^{2} GeV^{2}/c^{4}");
  Hxmmpiq_Sign_B3->Draw();
  cBelle->cd(8);
  Hxmmpiq_Sign_B4->GetYaxis()->SetTitleOffset(1.5);
  Hxmmpiq_Sign_B4->GetYaxis()->SetTitle(Form("Events/(%3.2f GeV/c^{2})",(m13max-m13min)/(scale1*10)));
  Hxmmpiq_Sign_B4->GetXaxis()->SetTitle("M(#psi(2S)#pi)^{2} GeV^{2}/c^{4}");
  Hxmmpiq_Sign_B4->Draw();
  sss<<PDFName<<"-cbelle"<<TXTMASSW<<TYPEFILE;
  gSystem->Exec(Form("mv %s %s.tmp",sss.str().c_str(),sss.str().c_str()));
  cBelle->Print(sss.str().c_str());
  sss.clear();
  sss.str("");
//=============================================================
  gStyle->SetOptStat("em");
  cHely->cd();
  HxCosK_Sign_H->SetLineColor(kBlue);
  HxHely_Sign->SetLineColor(kRed);
  HxCosK_Sign_R->SetLineColor(kMagenta);
  HxCosK_Sign_W->SetLineColor(kGreen);
  HxCosK_Sign_R->Scale(HxHely_Sign->Integral()/HxCosK_Sign_R->Integral()*0.05);
  HxCosK_Sign_W->Scale(HxHely_Sign->Integral()/HxCosK_Sign_W->Integral()*0.05);
  HxCosK_Sign_H->Draw();
//  TF1* FitBreitWigner = new TF1("FitBreitWigner",BreitWigner,-1.,1.,3);
  auto f2 = new TF1("f2","[3]*(TMath::Gaus(x,[0],[1])+[4]*TMath::Gaus(x,[0],[2]))",-1,1);
//  auto f2 = new TF1("f2","[5]*(ROOT::Math::crystalball_function(x,[0],[1],[2],[3])+[6]*TMath::Gaus(x,[3],[4]))",-1,1);
//  auto f2 = new TF1("f2","[7]*(ROOT::Math::crystalball_function(x,[0],[1],[2],[3])+[8]*ROOT::Math::crystalball_function(x,[4],[5],[6],[3]))",-1,1);
  Double_t parFit[2];
  parFit[0]=ang_ZMass;
  parFit[1]=0.2;
  parFit[2]=0.4;
  parFit[3]=1;
  parFit[4]=1;
  parFit[5]=1;
  parFit[6]=1;
  parFit[7]=1;
  parFit[8]=1;
//   parFit[0]=2;
//   parFit[1]=1;
//   parFit[2]=1;
//   parFit[3]=-0.4;
//   parFit[4]=0.4;
//   parFit[5]=1;
//   parFit[6]=1;
//   parFit[7]=1;
//   parFit[8]=1;
//  BreitWigner->SetNpx(1000);
  f2->SetParameters(parFit);
//   f2->SetParameter(0,2);
//   f2->SetParameter(1,1);
//   f2->SetParameter(2,1);
//   f2->SetParameter(3,-0.4);
  HxCosK_Sign_R->Fit("f2");
  parFit[0]=ang_ZMassWT;
  HxCosK_Sign_W->Fit("f2");
  HxCosK_Sign_R->SetMaximum(HxCosK_Sign_H->GetMaximum()*1.2);
  HxCosK_Sign_W->SetMaximum(HxCosK_Sign_H->GetMaximum()*1.2);
  HxCosK_Sign_R->Draw("same");
  HxCosK_Sign_W->Draw("same");
  HxHely_Sign  ->Draw("same");
  HxCosK_Sign_H->Draw("same");
  sss.clear();
  sss.str("");
  sss<<PDFName<<"-helicity"<<TXTMASSW<<TYPEFILE;
  gSystem->Exec(Form("mv %s %s.tmp",sss.str().c_str(),sss.str().c_str()));
  cHely->Print(sss.str().c_str());
  sss.clear();
  sss.str("");
  cHelyLeft->cd();
  HxHely_Left->SetLineColor(kRed);
  HxCosK_Left_R->SetLineColor(kBlue);
  HxCosK_Left_H->Draw();
  HxHely_Left->Draw("same");
//  HxCosK_Left_R->Draw("same");
//  HxHely_Left->Draw();
  sss.clear();
  sss.str("");
  sss<<PDFName<<"-helicity-Left"<<TXTMASSW<<TYPEFILE;
  gSystem->Exec(Form("mv %s %s.tmp",sss.str().c_str(),sss.str().c_str()));
  cHelyLeft->Print(sss.str().c_str());
  sss.clear();
  sss.str("");
  ckstarproj->cd();
  Hxmmpiq_Sign->Draw();
  Hxmmpiq_SBL->SetLineColor(kRed);
  Hxmmpiq_SBL->Draw("same");
//   HDalitzProy->Scale(Hxmmpiq_Sign->Integral()*Hxmmpiq_Sign->GetXaxis()->GetBinWidth(1)/HDalitzProy->Integral()/HDalitzProy->GetXaxis()->GetBinWidth(1));
//   HDalitzProy->SetLineColor(kMagenta);
//   HDalitzProy->Draw("same");
  sss.clear();
  sss.str("");
  sss<<PDFName<<"-kstarproj"<<TXTMASSW<<TYPEFILE;
  gSystem->Exec(Form("mv %s %s.tmp",sss.str().c_str(),sss.str().c_str()));
  ckstarproj->Print(sss.str().c_str());
  sss.clear();
  sss.str("");
  cMass->cd();
  HxMass->SetMinimum(0);;
  HxMass->GetYaxis()->SetTitleOffset(1.5);
  HxMass->GetYaxis()->SetTitle(Form("Events/(%3.2f GeV/c^{2})",(XMaxSign-XMinSign)/xMassHBin2));
  HxMass->GetXaxis()->SetTitle("B^{0} Mass (GeV/c^{2})");
  HxMass->Draw();
  HxMass_Sign->SetLineColor(kBlue);
//  HxMass->SetLineColor(kMagenta);
  HxMass_Sign->SetLineWidth(1.5);
  HxMass_Sign->Draw("same");
  HxMass_SBL->SetLineWidth(1.5);
  HxMass_SBR->SetLineWidth(1.5);
  HxMass_SBL->SetLineColor(kRed);
  HxMass_SBR->SetLineColor(kPink);
  HxMass_SBL->Draw("same");
  HxMass_SBR->Draw("same");
  TLegend* leg_Mass = new TLegend(0.54,0.68,0.70,0.85);
  leg_Mass->SetTextSize(0.025) ;
  leg_Mass->SetTextAlign(13);
  leg_Mass->SetBorderSize(0.);
  leg_Mass->SetFillStyle(0);
  leg_Mass->AddEntry(HxMass_Sign ,Form("Signal: %3.2f<mass<%3.2f GeV/c^{2}",WMinSign,WMaxSign),"l");
  leg_Mass->AddEntry(HxMass_SBL  ,Form("SB Left  : %3.2f<mass<%3.2f GeV/c^{2}",XMinSBL,XMaxSBL),"l");
  leg_Mass->AddEntry(HxMass_SBR  ,Form("SB Right : %3.2f<mass<%3.2f GeV/c^{2}",XMinSBR,XMaxSBR),"l");
  leg_Mass->Draw("same");
  sss.clear();
  sss.str("");
  sss<<PDFName<<"-cmass"<<TXTMASSW<<TYPEFILE;
  gSystem->Exec(Form("mv %s %s.tmp",sss.str().c_str(),sss.str().c_str()));
  cMass->Print(sss.str().c_str());
  sss.clear();
  sss.str("");
//=============================================================
  gStyle->SetOptStat(000000);
  cPeak->cd(1);
  sss<<PDFName<<"-cpeak"<<TXTMASSW<<TYPEFILE;
  gSystem->Exec(Form("mv %s %s.tmp",sss.str().c_str(),sss.str().c_str()));
  HxmmpiCosK_Sign->GetYaxis()->SetTitle("M(#psi(2S)#pi) (GeV/c^{2})");
  HxmmpiCosK_Sign->GetYaxis()->SetTitleOffset(1.5);
  HxmmpiCosK_Sign->GetXaxis()->SetTitle("cos(#theta_{K})");
  HxmmpiCosK_Sign->Draw("COLZ");
  cPeak->cd(2);
//  HkstarmassCosK_Sign->Draw("COLZ");
  TArrow *arrow04=new TArrow(-1,ZMass,hely_zed,ZMass);
  arrow04->SetLineColor(kRed);
  arrow04->Draw("same");
  HxmmpiCosK_Sign_Ks->GetXaxis()->SetTitle("cos(#theta_{K})");
  HxmmpiCosK_Sign_Ks->GetYaxis()->SetTitle("M(#psi(2S)#pi) (GeV/c^{2})");
  HxmmpiCosK_Sign_Ks->GetYaxis()->SetTitleOffset(1.5);
  HxmmpiCosK_Sign_Ks->Draw("COLZ");
  TArrow *arrow2=new TArrow(hely_zed,sqrt(m13min),hely_zed,ZMass);
  arrow2->SetLineColor(kRed);
  arrow2->SetArrowSize(0.01);
  arrow2->Draw("same,<");
  arrow04->SetArrowSize(0.01);
  arrow04->Draw("same,<");
  TLegend* leg_zed = new TLegend(0.50,0.75,0.90,0.80);
  leg_zed->SetTextSize(0.025) ;
  leg_zed->SetTextAlign(13);
  leg_zed->SetBorderSize(0.);
  leg_zed->SetFillStyle(0);
  leg_zed->AddEntry(arrow04 ,"mass Z(4430)","l");
  leg_zed->AddEntry(arrow2  ,Form( "cos(#theta_{K})(Z(4430),K^{*0}(892))"),"l");
  leg_zed->Draw("same");
  cPeak->cd(3);
  HxmmpiCosK_Sign_Md->Draw("COLZ");
  HxmmpiCosK_Sign_Md->GetYaxis()->SetTitleOffset(1.5);
  HxmmpiCosK_Sign_Md->GetYaxis()->SetTitle("M(#psi(2S)#pi) (GeV/c^{2})");
  HxmmpiCosK_Sign_Md->GetXaxis()->SetTitle("cos(#theta_{K})");
  TArrow *arrow3=new TArrow(hely_zmd,sqrt(m13min),hely_zmd,ZMass);
  arrow3->SetLineColor(kMagenta);
  arrow3->SetArrowSize(0.01);
  arrow3->Draw("same,<");
  TArrow *arrow02=new TArrow(-1.,ZMass,hely_zmd,ZMass);
  arrow02->SetLineColor(kRed);
  arrow02->SetArrowSize(0.01);
  arrow02->Draw("same,<");
  TLegend* leg_zmd = new TLegend(0.50,0.75,0.90,0.80);
  leg_zmd->SetTextSize(0.025) ;
  leg_zmd->SetTextAlign(13);
  leg_zmd->SetBorderSize(0.);
  leg_zmd->SetFillStyle(0);
  leg_zmd->AddEntry(arrow02 ,"mass Z(4430)","l");
  leg_zmd->AddEntry(arrow3  ,Form( "cos(#theta_{K})(Z(4430),K#pi)"),"l");
  leg_zmd->Draw("same");
  cPeak->cd(4);
/*   Hxmmpiq_Sign->SetMaximum(1.7*Hxmmpiq_Sign->GetMaximum());
  Hxmmpiq_Sign->Scale(1.3     /Hxmmpiq_Sign->Integral());
  Hxmmpiq_Sign_Ks->Scale(1  /Hxmmpiq_Sign_Ks->Integral());
  Hxmmpiq_Sign_NK->Scale(0.6 /Hxmmpiq_Sign_NK->Integral());
  Hxmmpiq_Sign_Md->Scale(0.4/Hxmmpiq_Sign_Md->Integral());
  Hxmmpiq_Sign_Tl->Scale(0.3/Hxmmpiq_Sign_Tl->Integral());
  Hxmmpiq_Sign->Draw("HIST");
  Hxmmpiq_Sign_Ks->SetLineColor(kBlue);
  Hxmmpiq_Sign_Md->SetLineColor(kRed);
  Hxmmpiq_Sign_NK->SetLineColor(kMagenta);
  Hxmmpiq_Sign_Tl->SetLineColor(kGreen);
  Hxmmpiq_Sign_Ks->Draw("HIST same");
  Hxmmpiq_Sign_Md->Draw("HIST same");
  Hxmmpiq_Sign_Tl->Draw("HIST same");
  Hxmmpiq_Sign_NK->Draw("HIST same");
 */ 
  Hxmmpi_Sign   ->Scale(1.3  /Hxmmpi_Sign->Integral());
  Hxmmpi_Sign_Ks->Scale(1.1  /Hxmmpi_Sign_Ks->Integral());
  Hxmmpi_Sign_NK->Scale(0.6  /Hxmmpi_Sign_NK->Integral());
  Hxmmpi_Sign_Md->Scale(0.4  /Hxmmpi_Sign_Md->Integral());
  Hxmmpi_Sign_Tl->Scale(0.2  /Hxmmpi_Sign_Tl->Integral()/(Hxmmpi_Sign_Tl->GetBinWidth(2)/Hxmmpi_Sign->GetBinWidth(2)));
  Hxmmpi_Sign->SetMaximum(1.5*Hxmmpi_Sign->GetMaximum());
  Hxmmpi_Sign->Draw("HIST");
  Hxmmpi_Sign->GetXaxis()->SetTitleOffset(0.8);
  Hxmmpi_Sign->GetXaxis()->SetTitle("M(#psi(2S)#pi) (GeV/c^{2})");
  Hxmmpi_Sign->GetYaxis()->SetTitle("arbitrary scale");
  Hxmmpi_Sign->GetYaxis()-> SetLabelSize(0);
  TArrow *arrowZ=new TArrow(ZMass,Hxmmpi_Sign->GetMinimum(),ZMass,0.7*Hxmmpi_Sign->GetMaximum());
  TPaveLabel *labelZ= new TPaveLabel(ZMass,0.72*Hxmmpi_Sign->GetMaximum(),ZMass+0.2,0.79*Hxmmpi_Sign->GetMaximum(),"mass Z(4430)");
  arrowZ->SetLineColor(kRed);
//  arrowZ->SetLineWidth(2);
//  arrowZ->SetAngle(40);
  arrowZ->SetArrowSize(0.01);
  arrowZ->Draw("same,<");
  labelZ->Draw("same");
  Hxmmpi_Sign_Ks->SetLineColor(kBlue);
  Hxmmpi_Sign_Md->SetLineColor(kRed);
  Hxmmpi_Sign_NK->SetLineColor(kMagenta);
  Hxmmpi_Sign_Tl->SetLineColor(kGreen);
  Hxmmpi_Sign_Ks->Draw("HIST same");
  Hxmmpi_Sign_Md->Draw("HIST same");
  Hxmmpi_Sign_Tl->Draw("HIST same");
  Hxmmpi_Sign_NK->Draw("HIST same");
//   TArrow *arrowz=new TArrow(ZMass,Hxmmpi_Sign->GetMinimum(),ZMass,Hxmmpi_Sign->GetMaximum());
//   arrowz->SetLineColor(kRed);
//   arrowz->Draw("same");
  Hxmmpi_Sign_NK->Draw("HIST same");
  TLegend* leg_mmpi = new TLegend(0.10,0.55,0.70,0.80);
  leg_mmpi->SetTextSize(0.025) ;
  leg_mmpi->SetTextAlign(13);
  leg_mmpi->SetBorderSize(0.);
  leg_mmpi->SetFillStyle(0);
  leg_mmpi->AddEntry(Hxmmpi_Sign ,"all events","l");
  leg_mmpi->AddEntry(Hxmmpi_Sign_Ks ,Form( "|mass K#pi-K^{*0}|<0.1 GeV/c^{2}"),"l");
  leg_mmpi->AddEntry(Hxmmpi_Sign_Md ,Form( "1<(mass K#pi)^{2}<1.75 GeV^{2}/c^{4}"),"l");
  leg_mmpi->AddEntry(Hxmmpi_Sign_Tl ,Form( "(mass K#pi)^{2}>1.75 GeV^{2}/c^{4}"),"l");
  leg_mmpi->AddEntry(Hxmmpi_Sign_NK ,Form( "K^{*0} veto"),"l");
  leg_mmpi->Draw("same");
  cPeak->cd(5);
  HxCosK_Sign->GetXaxis()->SetTitle("cos(#theta_{K})");
  HxCosK_Sign->GetYaxis()->SetTitle("arbitrary scale");
  HxCosK_Sign->GetYaxis()-> SetLabelSize(0);
  HxCosK_Sign   ->Scale(1.3/HxCosK_Sign->Integral());
  HxCosK_Sign_Ks->Scale(1.1/HxCosK_Sign_Ks->Integral());
  HxCosK_Sign_NK->Scale(0.6/HxCosK_Sign_NK->Integral());
  HxCosK_Sign_Md->Scale(0.4/HxCosK_Sign_Md->Integral());
  HxCosK_Sign_Tl->Scale(0.2/HxCosK_Sign_Tl->Integral()/(HxCosK_Sign_Tl->GetBinWidth(2)/HxCosK_Sign->GetBinWidth(2)));
  HxCosK_Sign->SetMaximum(1.3*HxCosK_Sign->GetMaximum());
  HxCosK_Sign->Draw("HIST");
  HxCosK_Sign_Ks->SetLineColor(kBlue);
  HxCosK_Sign_Ks->Draw("HIST same");
  HxCosK_Sign_Md->SetLineColor(kRed);
  HxCosK_Sign_Md->Draw("HIST same");
  HxCosK_Sign_Tl->SetLineColor(kGreen);
  HxCosK_Sign_Tl->Draw("HIST same");
  HxCosK_Sign_NK->SetLineColor(kMagenta);
  HxCosK_Sign_NK->Draw("HIST same");
  TArrow *arrowC=new TArrow(hely_zmd,HxCosK_Sign->GetMinimum(),hely_zmd,0.7*HxCosK_Sign->GetMaximum());
  arrowC->SetLineColor(kMagenta);
  arrowC->SetArrowSize(0.01);
  arrowC->Draw("same,<");
  TArrow *arrowK=new TArrow(hely_zed,HxCosK_Sign->GetMinimum(),hely_zed,0.76*HxCosK_Sign->GetMaximum());
  arrowK->SetLineColor(kRed);
  arrowK->SetArrowSize(0.01);
  arrowK->Draw("same,<");
  TPaveLabel *labelk= new TPaveLabel(hely_zed,0.79*HxCosK_Sign->GetMaximum(),hely_zed+0.8,0.86*HxCosK_Sign->GetMaximum(),"cos(#theta_{K})(Z(4430),K^{*0}(892))");
  TPaveLabel *labelc= new TPaveLabel(hely_zmd,0.72*HxCosK_Sign->GetMaximum(),hely_zmd+0.6,0.77*HxCosK_Sign->GetMaximum(),"cos(#theta_{K})(Z(4430),K#pi)");
  labelk->Draw();
  labelc->Draw();
  
  cPeak->Print(sss.str().c_str());
  sss.clear();
  sss.str("");
  sss<<PDFName<<"-clambda"<<TXTMASSW<<TYPEFILE;
  gSystem->Exec(Form("mv %s %s.tmp",sss.str().c_str(),sss.str().c_str()));
  cLambda->cd();
  HxLambdab_Sign->Draw();
  cLambda->Print(sss.str().c_str());
  sss.clear();
  sss.str("");
//   HxLambdab_Sign_c->Draw("same");
//   HxLambdab_Sign_c->SetLineColor(kMagenta);
//   HxLambdab_SBR->SetLineColor(kRed);
//   HxLambdab_SBL->SetLineColor(kBlue);
//   HxLambdab_SBL->Draw("same");
//   HxLambdab_SBR->Draw("same");
//   sss<<PDFName<<"-studies"<<TXTMASSW<<TYPEFILE;
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
//===============================================================================================================
double computeCosine (double Vx,double Vy,double Vz,double Wx, double Wy, double  Wz){

  double Vnorm = (Vx*Vx + Vy*Vy + Vz*Vz);
  double Wnorm = (Wx*Wx + Wy*Wy + Wz*Wz);
  double denom = sqrt(Vnorm*Wnorm);
  double VdotW =      Vx*Wx + Vy*Wy + Vz*Wz;
  double cosAlpha=0;
  if (Vnorm > 0. && Wnorm > 0.){
    cosAlpha = VdotW / denom;
//    cosAlpha = VdotW / (Vnorm * Wnorm);
  }else{
    std::cout<<Form("Error in computeCosine: cosAlpha=%f",cosAlpha)<<std::endl;
    cosAlpha = -99.;
  }  
    
  return cosAlpha;
}
//===============================================================================================================
//===============================================================================================================
//===============================================================================================================
double P_OnShell_Compute( double XM12, double XM1,double  XM2 ){
       double P1 = ( XM12*XM12 - ( XM1 + XM2)*( XM1 - XM2) );
              P1 = ( XM12*XM12 - ( XM1 - XM2)*( XM1 - XM2) )*P1;
	if(P1>0) {
	 P1 = sqrt( P1 )/( 2. * XM12);
	}else{
	  std::cout<<"Errore in P_OnShell_Compute: P1<0"<<std::endl;
	  exit(1);
	} 
	return P1;
}
//===============================================================================================================
double vettore(double XDmMass,double Am12,double Am13,double Am1,double Am2, double Am3, double Vemas, double GamVe) {
       double P1 = ( Vemas*Vemas - ( Am1 + Am2)*( Am1 + Am2) );
              P1 = ( Vemas*Vemas - ( Am1 - Am2)*( Am1 - Am2) )*P1;
	if(P1>0) {
	 P1 = sqrt( P1 )/( 2. * Vemas);
	}else{
	  std::cout<<"Errore in P_OnShell_Compute: P1<0=>"<<P1<<std::endl;
	  exit(1);
	} 
	double P_OnshellVe= P1 ;
        double Am_12;
        if(Am12>0){
	 Am_12=sqrt(Am12);
	}else{
	  std::cout<<"Errore in vettore: Am12<0"<<std::endl;
	  exit(1);
	} 
	
	double E1 = ( Am12 - Am2*Am2 + Am1*Am1 ) / ( 2. * Am_12 );
	double E3 = ( XDmMass*XDmMass - Am12 - Am3*Am3 ) / ( 2. * Am_12 );
	
	double cdota = .5 * ( Am1*Am1 + Am3*Am3 + 2.*E1*E3 - Am13 );

	double p1mod = E1*E1 - Am1*Am1;
	if(p1mod>0) {
	 p1mod=sqrt(p1mod);
	}else{
	  std::cout<<"Errore in vettore: p1mod<0"<<std::endl;
	  exit(1);
	} 
	double  prpi  =     ( ( XDmMass*XDmMass - ( Am_12 + Am3 )*( Am_12 + Am3 ) ) * \
	                      ( XDmMass*XDmMass - ( Am_12 - Am3 )*( Am_12 - Am3 ) )) ;
	if(prpi>0) {
	 prpi=sqrt(prpi)/( 2. * XDmMass );
	}else{
	  std::cout<<"Errore in vettore: prpi<0"<<std::endl;
	  exit(1);
	} 
			      
	double Rd    = 5.;
	double Rr    = 1.6;
	double FD    = 1. / sqrt( 1. + ( Rd * prpi )*( Rd * prpi ));
	//Blatt-Weisskopf
	double Frp   = 1. / sqrt( 1. + ( Rr * p1mod)*( Rr * p1mod) );
	double Frp0  = 1. / sqrt( 1. + ( Rr * P_OnshellVe)*( Rr * P_OnshellVe) );
	double Gamm  = GamVe * pow(( p1mod / P_OnshellVe ),3) * Frp*Frp / Frp0*Frp0 * Vemas / Am_12;
	double DenBW = ( Vemas*Vemas - Am12 )*( Vemas*Vemas - Am12 ) + (Gamm*Vemas)*(Gamm*Vemas);
	double Ang = -2*cdota;
	double AReBW = Ang * FD * Frp * ( Vemas*Vemas - Am12 ) / DenBW ;
	double AImBW = Ang * FD * Frp * Gamm*Vemas	       / DenBW;
	double ret = (AReBW*AReBW + AImBW*AImBW);
// 	if (ret<0){
// 	 std::cout<<Form("Errore in vettore:ret<0 AReBW=%f AReBW=%f ret=",AReBW,AImBW,ret)<<std::endl;
// 	 exit(1);
// 	}
	return ret ;

}
//---------------------------------------------------------------------------------------	
//---------------------------------------------------------------------------------------	
  double DalitzModel( double* x,  double *par){
  double BackgroundPlot( double m12, double m13, double DmMass, double m1 , double m2, double m3);
  double vettore(double XDmMass,double Am12,double Am13,double Am1,double Am2, double Am3, double Vemas, double GamVe);
  double m12 =x[0];
  double m13 =x[1];
  double DMass = par[0];
  double m1    = par[1];
  double m2    = par[2];
  double m3    = par[3];
  double Vemas = par[4];
  double GamVe = par[5];
  if (BackgroundPlot(  m12,  m13,  DMass,  m1 ,  m2, m3)<1.) return 0;
  if (fabs(Vemas-sqrt(m12))>0.150) return 0;
  double ret = vettore( DMass, m12, m13,  m1, m2,  m3,  Vemas,  GamVe);
  return ret;
}
//---------------------------------------------------------------------------------------	
//---------------------------------------------------------------------------------------	
  double DalitzModel2( double* x,  double *par){
  double Background( double m12, double m13, double DmMass, double m1 , double m2, double m3);
  double BackgroundPlot( double m12, double m13, double DmMass, double m1 , double m2, double m3);
  double vettore(double XDmMass,double Am12,double Am13,double Am1,double Am2, double Am3, double Vemas, double GamVe);
  double m12 =x[0];
  double m13 =x[1];
  double DMass = par[0];
  double m1    = par[1];
  double m2    = par[2];
  double m3    = par[3];
  double Ve1mas = par[4];
  double GamVe1 = par[5];
  double Ve2mas = par[6];
  double GamVe2 = par[7];
  if (BackgroundPlot(  m12,  m13,  DMass,  m1 ,  m2, m3)<1.) return 0;
  if (fabs(Ve1mas-sqrt(m12))>0.150) return 0;
//  if(m12<1.1||m12>1.8) continue 0;
//  double ret =   0;
  double A0q = BackgroundPlot(  m12,  m13,  DMass,  m1 ,  m2, m3)/10;
  double A1q = vettore( DMass, m12, m13,  m1, m2,  m3,  Ve1mas,  GamVe1)/100;
  double A2q = vettore( DMass, m13, m12,  m1, m3,  m2,  Ve2mas,  GamVe2);
  double A0 = sqrt(A0q);
  double A1 = sqrt(A1q);
  double A2 = sqrt(A2q);
//   if( A1<exp(-20) ) A1=0;
//   if( A2<exp(2) ) A2=0;
//   if( A1>exp(-20) ) A1=1;
//   if( A2>exp(2) ) A2=1;
  double ret = A2q;
//  double ret = A0q + A1q +  A2q + 2*A0*A1 + 2*A1*A2;
  
//  cout<<"A1q="<<A1q<<"  A2q="<<A2q<<"    A1="<<A1<<"  A2="<<A2<<" ret = "<<ret<<endl;
  return ret;
}
//---------------------------------------------------------------------------------------	
double lambda(double x,double y,double z){
 return x*x+y*y+z*z-2*x*y-2*y*z-2*z*x;
}
//---------------------------------------------------------------------------------------	
double cos12( double DmMass,double m1,double m2,double m3,double sm12,double sm13){
 double lambda(double x,double y,double z);
// double sigma1=sm23*sm23;
 double sigma2=sm13*sm13;
 double sigma3=sm12*sm12;
 double cos12 = 2*sigma3*(sigma2-m3*m3-m1*m1)-(sigma3+m1*m1-m2*m2)*(DmMass*DmMass-sigma3-m3*m3);
        cos12 = cos12/(sqrt(lambda(DmMass*DmMass,m3*m3,sigma3)*lambda(sigma3,m1*m1,m2*m2)));
 return cos12;
}
double BreitWigner( double* x,  double *par){

   return par[2]*TMath::BreitWigner(x[0],par[0],par[1]);

}
