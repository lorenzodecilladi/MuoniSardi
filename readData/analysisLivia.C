#if !defined(__CINT__) || defined(__MAKECINT__)

#include <stdio.h>
#include <endian.h>
#include "TString.h"
#include <Riostream.h>
#include "TStopwatch.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TFile.h"
#include "TTree.h"
#include "TVectorF.h"
#include "TStyle.h"

#endif


void plot(TString inputFilePath="readFile.root"){

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



  // TH1F *histQ     = new TH1F("histQ   "  , "histQ   "   ,  16,  -0.5,    15.5);
   TH1F *histCLOCK = new TH1F("histCLOCK" , "Scaler (TDC) [C257]"  , 110,  -0.5,   110.5);
   TH1F *histADC1  = new TH1F("histADC1"  , "ADCS1 [2249W]"   ,2048,  -0.5,  2048.5);
   TH1F *histADCG  = new TH1F("histADCG"  , "ADCSG [2249W]"   ,2048,  -0.5,  2048.5);
  // TH1F *histADCsum= new TH1F("histADCsum", "histADCsum" ,2048,  -0.5,  2048.5);
   TH1F *histTIME  = new TH1F("histTIME"  , "TDC[2228A]"   , 401,  -0.5,  4000.5);
  // TH1F *histINHIB = new TH1F("histINHIB", "histINHIB"   ,   2,   0. ,     2. );
  // TH1F *histDEADT = new TH1F("histDEADT", "histDEADTIME", 101,-0.005,   1.005);
  // TH1F *histPATT  = new TH1F("histPATT" , "histPATT"    ,  32,   0.5,    32.5);

   UInt_t nEvts  = readTree -> GetEntries();
  
   for(UInt_t i=0; i<nEvts; i++){
     readTree   -> GetEvent(i);
  //   for(int j=0; j<16; j++){
  //     histQ->Fill(j,validStr[j]);
  //   }
     histCLOCK  -> Fill(CLOCKch);
     histTIME   -> Fill(TDCch);
     histADC1   -> Fill(ADCchS1);
     histADCG   -> Fill(ADCchSG);
  //   histADCsum -> Fill(ADCchS1+ADCchSG);
  //   for(int j=0; j<16; j++){
  //     histPATT -> Fill(j+1   , pattReg[j]);
  //     histPATT -> Fill(j+16+1, pattMul[j]);
  //   }
  //   histINHIB  -> Fill(0.5, nonInhib);
  //   histINHIB  -> Fill(1.5, inhib   );
  //   histDEADT  -> Fill(deadTime     );
   }
  
  //plots
    TCanvas *canvas1  = new TCanvas("TDC", "TDC", 200, 10, 600, 400);
    gPad->SetLogy();
    gStyle->SetOptStat("emr");
    histTIME -> GetXaxis()-> SetTitle("# canali");
    histTIME -> GetYaxis()-> SetTitle("# eventi");
    histTIME -> DrawCopy("hist");
  
    TCanvas *canvas2 = new TCanvas("CLOCK", "CLOCK", 200, 10, 600, 400);
    gPad->SetLogy();
    histCLOCK -> GetXaxis()->SetTitle("#canali");  
    histCLOCK -> GetYaxis()->SetTitle("#eventi");
    histCLOCK-> DrawCopy("hist");
  
    TCanvas *canvas3 = new TCanvas("ADC1", "ADC1", 200, 10, 600, 400);
    // gPad->SetLogy();
    histADC1 -> GetXaxis()->SetTitle("# canali");
    histADC1 -> GetYaxis()->SetTitle("# eventi");    
    histADC1->GetYaxis()->SetRangeUser(0.,3000.);
    histADC1->GetXaxis()->SetRangeUser(0.,2025.);
    histADC1 -> DrawCopy("hist");
  
  
    TCanvas *canvas4 = new TCanvas("ADCG", "ADCG", 200, 10, 600, 400);
    // gPad->SetLogy();
    histADCG -> GetXaxis()->SetTitle("# canali");
    histADCG -> GetYaxis()->SetTitle("# eventi");
    histADCG->GetYaxis()->SetRangeUser(0.,1400.);
    histADCG->GetXaxis()->SetRangeUser(0.,2025.);  //c'è un picco al canale 2030
    histADCG -> DrawCopy("hist");


   //canvas   -> Write();
   //vec.Write("vec");
   //readFile -> Write();
    readFile -> Close();
  
   //histADCsum->DrawCopy("hist");
  
}
