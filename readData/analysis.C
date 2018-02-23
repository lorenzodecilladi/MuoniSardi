#if !defined(__CINT__) || defined(__MAKECINT__)
#include <stdio.h>
#include <endian.h>
#include "TString.h"
#include <Riostream.h>
#include "TStopwatch.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TF1.h"
#include "TFile.h"
#include "TTree.h"
#include "TVectorF.h"
#include "TStyle.h"
#include "TMath.h"
#include "TGraphErrors.h"
#endif


void analysis(TString inputFilePath="readFile.root"){

  //open readFile
  TFile   *readFile = new TFile(inputFilePath);
  if(readFile -> IsZombie()){
    cout << "File " << inputFilePath << " not found!!" << endl;
    return;
  }

  TTree   *readTree = (TTree*)readFile -> Get("readTree");
  int validStr[16];
  int CLOCKch; //OK
  int SCLinh;       //scaler inhibited
  int SCLuninh;     //scaler not inhibited
  int TDCch;   //OK
  int ADCchS1; //OK
  int ADCchSG; //OK
  int pattReg[16];
  int pattMul[16];
  TBranch *bvalidStr = readTree -> GetBranch("validStr");
  TBranch *bCLOCKch  = readTree -> GetBranch("CLOCKch" );
  TBranch *bSCLinh   = readTree -> GetBranch("SCLinh"  );
  TBranch *bSCLuninh = readTree -> GetBranch("SCLuninh");
  TBranch *bTDCch    = readTree -> GetBranch("TDCch"   );
  TBranch *bADCchS1  = readTree -> GetBranch("ADCchS1" );
  TBranch *bADCchSG  = readTree -> GetBranch("ADCchSG" );
  TBranch *bpattReg  = readTree -> GetBranch("pattReg" );
  TBranch *bpattMul  = readTree -> GetBranch("pattMul" );
  bvalidStr -> SetAddress(&validStr);
  bCLOCKch  -> SetAddress(&CLOCKch );
  bSCLinh   -> SetAddress(&SCLinh  );
  bSCLuninh -> SetAddress(&SCLuninh);
  bTDCch    -> SetAddress(&TDCch   );
  bADCchS1  -> SetAddress(&ADCchS1 );
  bADCchSG  -> SetAddress(&ADCchSG );
  bpattReg  -> SetAddress(&pattReg );
  bpattMul  -> SetAddress(&pattMul );

  TVectorF *vec     = (TVectorF*)readFile->Get("vec");
  Double_t nonInhib = (*vec)[0];
  Double_t inhib    = (*vec)[1];
  Double_t deadTime = (*vec)[2];

  Int_t nbinsTIME   = 401;
  Int_t nbinsADC    = 2049;
  TH1F *histCLOCK   = new TH1F("histCLOCK"  , "hist Scaler (TDC) [C257]", 111      ,  -0.5,   110.5);
  TH1F *histADC1    = new TH1F("histADC1"   , "ADCS1 [2249W]"      , nbinsADC ,  -0.5,  2048.5);
  TH1F *histADCG    = new TH1F("histADCG"   , "ADCSG [2249W]"      , nbinsADC ,  -0.5,  2048.5);
  TH1F *histADC1nP  = new TH1F("histADC1nP" , "ADCS1 [2249W] nP"   , nbinsADC ,  -0.5,  2048.5);
  TH1F *histADCGnP  = new TH1F("histADCGnP" , "ADCSG [2249W] nP"   , nbinsADC ,  -0.5,  2048.5);
  TH1F *histTIME    = new TH1F("histTIME"   , "hist TDC[2228A]"         , nbinsTIME,  -0.5,  4010.5);
  TH1F *histTIMEln  = new TH1F("histTIMEln" , "hist TDC[2228A] (ln)"    , nbinsTIME,  -0.5,  4010.5);

  UInt_t  nEvts     = readTree -> GetEntries();
  Float_t binW      = 10;
  Float_t muonRatio = 1.268; //N(mu+)/N(mu-)
  Float_t rho       = muonRatio;

  int validCheck[16] = {1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0};
  int validCounter   = 0;

  int vetoCounter = 0;
  
  for(UInt_t i=0; i<nEvts; i++){
    readTree   -> GetEvent(i);

    //check for validation string
    Bool_t validation = kTRUE;
    for(UInt_t j=0; j<16; j++){
      if(validStr[j]!=validCheck[j])validation=kFALSE;      
    }
    
    if(validation==kFALSE){ //exclude event if validation string is not OK
      validCounter++;
    }
    else{
      
      
      //if(pattReg[0]!=1){ //elimina eventi vetati da S2 (offline)
      //if(pattReg[0]==1){ //   solo eventi vetati da S2 (offline)
      //vetoCounter++;

      
      //if((ADCchSG>330 && ADCchSG<890) || ADCchSG>2025){ //"elettroni" + overflow
      //if(ADCchSG>890 && ADCchSG<2025){ //muoni
      //if(ADCchS1>330 && ADCchS1<1000){
      // if(TDCch<180){ //cattura (solo tempi bassi, 180 ch = 0.45 us = 450 ns
      //if(TDCch>180){ //solo decadimento (e overflow dell'ADC)
      //if(TDCch>180 && TDCch<4000){ //solo decadimento, tagliato anche l'overflow dell'ADC
      histCLOCK  -> Fill(CLOCKch);
      histTIME   -> Fill(TDCch);
      histADC1   -> Fill(ADCchS1);
      histADCG   -> Fill(ADCchSG);
      //}
      //}
    }
  }
  
  Int_t pedestalS1 = 159;
  Int_t pedestalSG = 264;
  
  for(Int_t i=0; i<nbinsADC; i++){
    histADC1nP->SetBinContent(i+1, histADC1->GetBinContent(pedestalS1+i+1));
    histADCGnP->SetBinContent(i+1, histADCG->GetBinContent(pedestalSG+i+1));
    if(histTIME->GetBinContent(i+1)>0){
      histTIMEln->SetBinContent(i+1, TMath::Log(histTIME->GetBinContent(i+1)));
      histTIMEln->SetBinError(i+1, TMath::Sqrt(histTIME->GetBinContent(i+1))/histTIME->GetBinContent(i+1) );
    }
  }

  cout << "[INFO    ]" << endl;
  cout << "[INFO    ] Events not validated: " << validCounter << endl;
  cout << "[INFO    ]" << endl;
  cout << "[INFO    ]Events not vetoed by S2 (offline, through pattern unit): " << vetoCounter << endl;
  cout << "[INFO    ]" << endl;
  
  
  //=======================================
  //---------------- TDC ------------------
  //=======================================
  TCanvas *canvas1  = new TCanvas("TDC", "TDC", 200, 10, 600, 400);
  canvas1  ->cd();
  gPad     -> SetLogy();
  gStyle   ->SetOptStat("emr");
  histTIME -> GetXaxis() -> SetTitle("# canali");
  histTIME -> GetYaxis() -> SetTitle("# eventi");
  //histTIME -> GetXaxis() -> SetRangeUser(0.,2000.);   
  
  /*
  //---------------------------------------
  cout << "---------------------------------------" << endl;
  cout << "DECAY - lifetime hardcoded"              << endl;
  TF1 *fitDec2_2   = new TF1("fitDec2_2","[0]+[1]*TMath::Exp(- 2.197*0.025*(x))");
   fitDec2_2 -> SetParLimits(1, 0, 10000000);
   histTIME  -> Fit(fitDec2_2, "ME", "", 180, 800); //400, 1200
   cout<<"Chi^2: "<<fitDec2_2->GetChisquare()<<", number of DoF: "<<fitDec2_2->GetNDF()<<" (Probability: "<<fitDec2_2->GetProb()<<")."<<endl;
   */

   
   //---------------------------------------
   cout << "---------------------------------------" << endl;
   cout << "DECAY - esponenziale semplice"           << endl;
   
   TF1 *fitDec      = new TF1("fitDec","[0]+[1]*TMath::Exp(-[2]*0.0025*(x))");
   fitDec   -> SetParameter(0, 10.        );
   fitDec   -> SetParLimits(0, 0, 10000000);
   fitDec   -> SetParLimits(1, 0, 10000000);
   histTIME -> Fit(fitDec, "ME", "",520, 2040); //400, 1200
   Double_t parDec1 = fitDec->GetParameter(1);
   Double_t parDec2 = fitDec->GetParameter(2);
   Double_t lambda0 = parDec2;
   cout<<"Chi^2: "<<fitDec->GetChisquare()<<", number of DoF: "<<fitDec->GetNDF()<<" (Probability: "<<fitDec->GetProb()<<")."<<endl;


   //---------------------------------------
   cout << "CAPTURE - esponenziale semplice" << endl;
   cout << "---------------------------------------" << endl;
   TF1 *fitCap0     = new TF1("fitCap0","[0]+[1]*TMath::Exp(-[2]*0.0025*(x))");
   fitCap0  -> SetParLimits(1, 0, 10000000);
   fitCap0  -> SetLineColor(kBlue);
   fitCap0  -> SetLineWidth(2);
   histTIME -> Fit(fitCap0, "ME+", "", 50, 110);
   cout<<"Chi^2: "<<fitCap0->GetChisquare()<<", number of DoF: "<<fitCap0->GetNDF()<<" (Probability: "<<fitCap0->GetProb()<<")."<< endl;
   Double_t parCap02 = fitCap0->GetParameter(2);
   
   histTIME -> Draw();

   
   cout << "[RESULT   ]" << endl;
   cout << "[RESULT   ] Muon mean lifetime            : (" << 1/parDec2 << "+-" << (1/(fitDec->GetParameter(2)*fitDec->GetParameter(2)))*fitDec->GetParError(2) << ") us" <<endl;
   cout << "[RESULT   ] Muon capture time             : (" << 1./parCap02 << "+-" << (1/(fitCap0->GetParameter(2)*fitCap0->GetParameter(2)))*fitCap0->GetParError(2) << ") us ERRORE SBAGLIATO" <<endl;   
   cout << "[RESULT   ]" << endl;



   //=======================================
   //---------------- Ln(TDC) ------------------
   //=======================================
   TCanvas *canvas1ln  = new TCanvas("TDCln", "TDCln", 200, 10, 600, 400);
   canvas1ln  -> cd();
   gStyle     -> SetOptStat("emr");
   histTIMEln -> GetXaxis() -> SetTitle("# canali");
   histTIMEln -> GetYaxis() -> SetTitle("# eventi");
   histTIMEln -> Draw("E");

   //---------------------------------------
   cout << "---------------------------------------" << endl;
   cout << "DECAY - linear fit"                      << endl;
   
   TF1 *fitDecLin   = new TF1("fitDec","[0]-[1]*0.0025*x");
   histTIMEln -> Fit(fitDecLin, "ME", "",520, 2040); //400, 1200
   cout<<"Chi^2: "<<fitDecLin->GetChisquare()<<", number of DoF: "<<fitDecLin->GetNDF()<<" (Probability: "<<fitDecLin->GetProb()<<")."<<endl;

   //---------------------------------------
   cout << "CAPTURE - linear fit" << endl;
   cout << "---------------------------------------" << endl;
   TF1 *fitCapLin0  = new TF1("fitCapLin0","[0]-[1]*0.0025*x");
   fitCapLin0  -> SetLineColor(kBlue);
   fitCapLin0  -> SetLineWidth(2);
   histTIMEln -> Fit(fitCapLin0, "ME+", "", 50, 110);
   cout<<"Chi^2: "<<fitCapLin0->GetChisquare()<<", number of DoF: "<<fitCapLin0->GetNDF()<<" (Probability: "<<fitCapLin0->GetProb()<<")."<< endl;

         cout << "[RESULT   ]" << endl;
	 cout << "[RESULT   ] Muon mean lifetime            : (" << 1./fitDecLin->GetParameter(1) << "+-" << (1/(fitDec->GetParameter(2)*fitDec->GetParameter(2)))*fitDec->GetParError(2) << ") us" <<endl;
	 cout << "[RESULT   ] Muon capture time             : (" << 1./fitCapLin0->GetParameter(1) << "+-" << (1/(fitCap0->GetParameter(2)*fitCap0->GetParameter(2)))*fitCap0->GetParError(2) << ") us ERRORE SBAGLIATO" <<endl;   
   cout << "[RESULT   ]" << endl;

   


   
   
   
   
   //=======================================
   //--------------- CLOCK -----------------
   //=======================================
   TCanvas *canvas2  = new TCanvas("CLOCK", "CLOCK", 200, 10, 600, 400);
   canvas2  ->cd();
   gPad     -> SetLogy();
   gStyle   ->SetOptStat("emr");
   histCLOCK -> GetXaxis()->SetTitle("#canali");  
   histCLOCK -> GetYaxis()->SetTitle("#eventi");
   //histCLOCK -> GetXaxis()->SetRangeUser(0.,500.);
   
   //---------------------------------------
   cout << "---------------------------------------" << endl;
   cout << "DECAY - esponenziale semplice"           << endl;
   
   TF1 *fitDecCL = new TF1("fitDecCL","[0]+[1]*TMath::Exp(-[2]*0.100*(x))");
   fitDecCL  -> SetParLimits(1, 0, 10000000  );
   histCLOCK -> Fit(fitDecCL, "ME", "", 13, 51); //4, 20
   Double_t parDec1CL = fitDecCL->GetParameter(1);
   Double_t parDec2CL = fitDecCL->GetParameter(2);
   Double_t lambda0CL = parDec2CL;
   cout<<"Chi^2: "<<fitDecCL->GetChisquare()<<", number of DoF: "<<fitDecCL->GetNDF()<<" (Probability: "<<fitDecCL->GetProb()<<")."<<endl;
   
   //---------------------------------------
   cout << "CAPTURE - esponenziale semplice" << endl;
   cout << "---------------------------------------" << endl;
   TF1 *fitCapCL0   = new TF1("fitCapCL0","[0]+[1]*TMath::Exp(-[2]*0.100*(x))");
   fitCapCL0 -> SetParLimits(0, 0, 10000000);
   fitCapCL0 -> SetParLimits(1, 0, 10000000);
   fitCapCL0 -> SetLineColor(kBlue);
   fitCapCL0 -> SetLineWidth(2);
   histCLOCK -> Fit(fitCapCL0, "ME+", "", 2, 5);
   cout<<"Chi^2: "<<fitCapCL0->GetChisquare()<<", number of DoF: "<<fitCapCL0->GetNDF()<<" (Probability: "<<fitCapCL0->GetProb()<<")."<< endl;
   Double_t parCap02CL  = fitCapCL0->GetParameter(2);
   
   histCLOCK -> Draw();

   
   cout << "[RESULT   ]" << endl;
   cout << "[RESULT   ] Muon mean lifetime            : (" << 1/(fitDecCL->GetParameter(2)) << "+-" << (1/(fitDecCL->GetParameter(2)*fitDecCL->GetParameter(2)))*fitDecCL->GetParError(2) << ") us" <<endl;
   cout << "[RESULT   ] Muon capture time             : (" << 1./parCap02CL << "+-" << (1/(fitCapCL0->GetParameter(2)*fitCapCL0->GetParameter(2)))*fitCapCL0->GetParError(2) << ") us ERRORE SBAGLIATO" <<endl;   
   cout << "[RESULT   ]" << endl;
   
   Int_t    ndataCLOCK = 111;
   Float_t  CLOCKsamp[111] = {0.};
   Float_t  sCLOCKsamp[111] = {0.};
   Float_t  CLOCKarray[111] = {0.};
   Float_t  sCLOCKarray[111] = {0.};
   
   //Float_t   TIMEarray[5000] = {0.};
   //Float_t   sTIMEarray[5000] = {0.};
      
   //Float_t TIMElnarray[5000] = {0.};
   //Float_t sTIMElnarray[5000] = {0.};

   //cout << histCLOCK ->GetSize()-2 << endl;
   for(Int_t i=0; i<histCLOCK ->GetSize()-2; i++){
     CLOCKsamp[i]   = 0.100 * i; //[us]
     sCLOCKsamp[i]  = 0.100/2;  //[us] così l'intervallo di errore corrisponde alla "risoluzione(???)" del TDC
     CLOCKarray[i]  = histCLOCK ->GetBinContent(i+1);
     cout << CLOCKarray[i] << endl;
     sCLOCKarray[i] = histCLOCK ->GetBinError(i+1);
   }

   /*
   for(UInt_t i=0; i<(histTIME  ->GetSize()-2); i++){
     TIMEarray[i] = histTIME  ->GetBinContent(i+1);
     sTIMEarray[i] = histTIME  ->GetBinError(i+1);
   }
   
   for(UInt_t i=0; i<(histTIMEln->GetSize()-2); i++){
     TIMElnarray[i] = histTIMEln->GetBinContent(i+1);
     sTIMElnarray[i] = histTIMEln->GetBinError(i+1);
   }
   */
   
   //TGraphErrors
   //TDC
   TCanvas *gcCLOCK  = new TCanvas("gcCLOCK", "gcCLOCK", 200, 10, 600, 400);
   gcCLOCK->cd();
   //gPad->SetLogy();
   TGraphErrors *gCLOCK  = new TGraphErrors(ndataCLOCK, CLOCKsamp, CLOCKarray, sCLOCKsamp, sCLOCKarray);
   gCLOCK->SetTitle("graph Scaler (TDC) [C257]");
   gCLOCK->SetMarkerColor(kBlue);
   gCLOCK->SetLineColor(kBlue);
   gCLOCK->SetFillColor(38);
   //gCLOCK->GetYaxis()->SetRangeUser(0., 100000);
   gCLOCK->Draw("AP");
   //gcCLOCK->();
   //TGraphErrors *gTIME   = new TGraphErrors("gTIME"  , "graph TDC[2228A]"         );
   //TGraphErrors *gTIMEln = new TGraphErrors("gTIMEln", "graph TDC[2228A] (ln)"    );
   //gPad->Modified();
   //gCLOCK->SetMinimum(0.);

   TF1 *gfitDecCL = new TF1("gfitDecCL","[0]+[1]*TMath::Exp(-[2]*x)");
   gfitDecCL -> SetParLimits(1, 0, 10000000  );
   gCLOCK    -> Fit(gfitDecCL, "ME", "", 1.3, 5.1); //4, 20
   //Double_t parDec1CL = fitDecCL->GetParameter(1);
   //Double_t parDec2CL = fitDecCL->GetParameter(2);
   //Double_t lambda0CL = parDec2CL;
   cout<<"Chi^2: "<<fitDecCL->GetChisquare()<<", number of DoF: "<<fitDecCL->GetNDF()<<" (Probability: "<<fitDecCL->GetProb()<<")."<<endl;






   
   
   

   
   /*
   TCanvas *canvas3 = new TCanvas("ADC1", "ADC1", 200, 10, 600, 400);
   // gPad->SetLogy();
   histADC1 -> GetXaxis()->SetTitle("# canali");
   histADC1 -> GetYaxis()->SetTitle("# eventi");    
   //histADC1->GetYaxis()->SetRangeUser(0.,3000.);
   histADC1 -> Draw("E");
     
   TCanvas *canvas3nP = new TCanvas("ADC1nP", "ADC1nP", 200, 10, 600, 400);
   // gPad->SetLogy();

   TF1 *fitGausS1 = new TF1 ("gaus", "gaus");
   histADC1nP -> Fit(fitGausS1, "ME", "", 20, 90);
  
   TF1 *fitLanMuS1 = new TF1 ("landau", "landau");
   histADC1nP -> Fit(fitLanMuS1, "ME+", "", 250, 550);
   
   histADC1nP -> GetXaxis()->SetTitle("# canali");
   histADC1nP -> GetYaxis()->SetTitle("# eventi");    
   histADC1nP->GetYaxis()->SetRangeUser(0.,800.);
   histADC1nP -> Draw("E");
     
   TCanvas *canvas4 = new TCanvas("ADCG", "ADCG", 200, 10, 600, 400);
   // gPad->SetLogy();
   histADCG -> GetXaxis()->SetTitle("# canali");
   histADCG -> GetYaxis()->SetTitle("# eventi");
   //histADCG->GetYaxis()->SetRangeUser(0.,1400.);
   //histADCG->GetXaxis()->SetRangeUser(0.,2025.);  //c'è un picco al canale 2030
   histADCG -> Draw("E");
   
   TCanvas *canvas4nP = new TCanvas("ADCGnP", "ADCGnP", 200, 10, 600, 400);
   // gPad->SetLogy();

   TF1 *fitGaus = new TF1 ("gaus", "gaus");
   histADCGnP -> Fit(fitGaus, "ME", "", 0, 60);
  
   TF1 *fitLanEl = new TF1 ("landau", "landau");
   histADCGnP -> Fit(fitLanEl, "ME+", "", 90, 450);

   //TF1 *fitLanMuSG = new TF1 ("landau", "landau");
   //histADCGnP -> Fit(fitLanMuSG, "", "", 750, 1400);

   histADCGnP -> GetXaxis()->SetTitle("# canali");
   histADCGnP -> GetYaxis()->SetTitle("# eventi");
   histADCGnP->GetYaxis()->SetRangeUser(0.,900.);
   histADCGnP->GetXaxis()->SetRangeUser(0.,2025.);  //c'è un picco al canale 2030
   histADCGnP -> Draw("E");
   
   */

   //  readFile -> Close();   
}
