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

#endif


void analysisLivia(TString inputFilePath="readFile.root"){

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



  Int_t nbinsTIME = 401;
  Int_t nbinsADC = 2048;
  TH1F *histCLOCK   = new TH1F("histCLOCK" , "Scaler (TDC) [C257]"  , 110,  -0.5,   110.5);
  TH1F *histADC1    = new TH1F("histADC1"  , "ADCS1 [2249W]"   ,nbinsADC,  -0.5,  2048.5);
  TH1F *histADCG    = new TH1F("histADCG"  , "ADCSG [2249W]"   ,nbinsADC,  -0.5,  2048.5);
  TH1F *histADC1nP  = new TH1F("histADC1nP", "ADCS1 [2249W] nP",nbinsADC,  -0.5,  2048.5);
  TH1F *histADCGnP  = new TH1F("histADCGnP", "ADCSG [2249W] nP",nbinsADC,  -0.5,  2048.5);
  TH1F *histTIME    = new TH1F("histTIME"   , "TDC[2228A]"   , nbinsTIME,  -0.5,  4010.5);
  TH1F *histTIMEmus = new TH1F("histTIMEmus", "TDC[2228A]"   , nbinsTIME,  -0.5,  4010.5);

  
  UInt_t  nEvts     = readTree -> GetEntries();
  Float_t binW      = 10;
  Float_t muonRatio = 1.268; //N(mu+)/N(mu-)
  Float_t rho       = muonRatio;

  /*
  //FIT COMPLESSIVO DALL'ARTICOLO DELLA MARCELLO
  TF1 *fitF1 = new TF1("fitF1","[0]+[1]*TMath::Exp(-[2]*10*x) * ( [2]*(1.268+TMath::Exp(-[2]*10*x)) +[3]*TMath::Exp(-[3]*10*x) )",0,400);
  fitF1->SetParameter(0, 10.);
  fitF1->SetParameter(1, 880000000);
  fitF1->SetParameter(2, 0.01136);
  fitF1->SetParameter(3, 0.000094);
  //FUNZIONE DI FIT CON PARAMETRI TUNATI PER ITERAZIONE...manuale...!
  TF1 *fitF2 = new TF1("fitF2","10 + (2000000000.*(0.1000/4)/(1+1.268)) * TMath::Exp(-(0.1000/4)*(10*x)/2.2) * ( (1/2.2) * (1.268+TMath::Exp(-(0.1000/4)*(10*x)/2.2))  + 0.00376*TMath::Exp(-0.00376*(0.1000/4)*(10*x)))",0.1,1000);
  */


  //FIT se avessimo solo mu- .....
  TF1 *fitBrutto = new TF1("fitBrutto","[0]+[1]*TMath::Exp(-[2]*x)", 0, 200);
  fitBrutto -> SetParLimits(1, 0, 10000000);
  

  
   for(UInt_t i=0; i<nEvts; i++){
     readTree   -> GetEvent(i);
     histCLOCK  -> Fill(CLOCKch);
     histTIME   -> Fill(TDCch);
     histADC1   -> Fill(ADCchS1);
     histADCG   -> Fill(ADCchSG);
   }
   
   Int_t pedestalS1 = 159;
   Int_t pedestalSG = 264;
   
   for(Int_t i=0; i<nbinsADC; i++){
     histADC1nP->Fill(i+1, histADC1->GetBinContent(pedestalS1+i+1));
     histADCGnP->Fill(i+1, histADCG->GetBinContent(pedestalSG+i+1));
   }
   
   
   //plots
   
   TCanvas *canvas1  = new TCanvas("TDC", "TDC", 200, 10, 600, 400);
   canvas1->cd();
   //gPad->SetLogy();
   gStyle->SetOptStat("emr");
   histTIME -> GetXaxis()-> SetTitle("# canali");
   histTIME -> GetYaxis()-> SetTitle("# eventi");
   histTIME -> GetXaxis()->SetRangeUser(0.,1025.);   
   histTIME -> Fit(fitBrutto, "R");
   histTIME -> Draw();


   /*
   //per il fit complessivo
   cout << "[RESULT   ]" << endl;
   cout << "[RESULT   ] Muon mean lifetime: (" << (0.1/4)/(fitF1->GetParameter(2)) << "+-" << (1/(fitF1->GetParameter(2)*fitF1->GetParameter(2)))*fitF1->GetParError(2) << ") us" <<endl;
   cout << "[RESULT   ]" << endl;
  
   
   TCanvas *canvas2 = new TCanvas("CLOCK", "CLOCK", 200, 10, 600, 400);
   //gPad->SetLogy();
   histCLOCK -> GetXaxis()->SetTitle("#canali");  
   histCLOCK -> GetYaxis()->SetTitle("#eventi");
   histCLOCK-> DrawCopy("hist");
   //fitF2    -> Draw("same");

   
   TCanvas *canvas3 = new TCanvas("ADC1", "ADC1", 200, 10, 600, 400);
   // gPad->SetLogy();
   histADC1 -> GetXaxis()->SetTitle("# canali");
   histADC1 -> GetYaxis()->SetTitle("# eventi");    
   //histADC1->GetYaxis()->SetRangeUser(0.,3000.);
   histADC1 -> DrawCopy("hist");

   TCanvas *canvas3nP = new TCanvas("ADC1nP", "ADC1nP", 200, 10, 600, 400);
   // gPad->SetLogy();
   histADC1nP -> GetXaxis()->SetTitle("# canali");
   histADC1nP -> GetYaxis()->SetTitle("# eventi");    
   //histADC1nP->GetYaxis()->SetRangeUser(0.,3000.);
   histADC1nP -> DrawCopy("hist");
   
   
   TCanvas *canvas4 = new TCanvas("ADCG", "ADCG", 200, 10, 600, 400);
   // gPad->SetLogy();
   histADCG -> GetXaxis()->SetTitle("# canali");
   histADCG -> GetYaxis()->SetTitle("# eventi");
   //histADCG->GetYaxis()->SetRangeUser(0.,1400.);
   //histADCG->GetXaxis()->SetRangeUser(0.,2025.);  //c'è un picco al canale 2030
   histADCG -> DrawCopy("hist");

   TCanvas *canvas4nP = new TCanvas("ADCGnP", "ADCGnP", 200, 10, 600, 400);
   // gPad->SetLogy();
   histADCGnP -> GetXaxis()->SetTitle("# canali");
   histADCGnP -> GetYaxis()->SetTitle("# eventi");
   //histADCGnP->GetYaxis()->SetRangeUser(0.,1400.);
   //histADCGnP->GetXaxis()->SetRangeUser(0.,2025.);  //c'è un picco al canale 2030
   histADCGnP -> DrawCopy("hist");
   */ 
   //  readFile -> Close();
   
}
