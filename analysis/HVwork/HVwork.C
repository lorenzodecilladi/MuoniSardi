#if !defined(__CINT__) || defined(__MAKECINT__)

#include <Riostream.h>
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TMath.h"
#include "TH1D.h"
#include "TF1.h"
#include "TStopwatch.h"
#include "TGraph.h"
#include "TVectorD.h"
#include "TCanvas.h"
#include <TAxis.h>                        // ci serve per manipolare gli assi dei grafici
#include "TVirtualFitter.h"               // ci serve come supporto per i fit
#include <string>
#include <TLegend.h>
#include <TLine.h>
using namespace std;

#endif


void plot(TString inputFileName = "HVfile.txt"){
  
  TStopwatch watch;
  watch.Start(kTRUE);
  
  const UInt_t ndata    = 12;
  Int_t     HVS1[ndata] = {0};
  Float_t  sHVS1[ndata] = {0.};
  Int_t countsS1[ndata] = {0};
  Float_t sCountsS1[ndata] = {0};
  Float_t rateS1[ndata] = {0};
  Int_t     HVS2[ndata] = {0};
  Float_t  sHVS2[ndata] = {0.};
  Int_t countsS2[ndata] = {0};
  Float_t sCountsS2[ndata] = {0};
  Float_t rateS2[ndata] = {0};
  Int_t     HVS3[ndata] = {0};
  Float_t  sHVS3[ndata] = {0.};
  Int_t countsS3[ndata] = {0};
  Float_t sCountsS3[ndata] = {0};
  Float_t rateS3[ndata] = {0};
  Int_t     HVSG[ndata] = {0};
  Float_t  sHVSG[ndata] = {0.};
  Int_t countsSG[ndata] = {0};
  Float_t sCountsSG[ndata] = {0};
  Float_t rateSG[ndata] = {0};

  Int_t   AQtime[ndata] = {400, 400, 200, 200, 200, 200, 200, 200, 200, 200, 200, 200};
  
  
  ifstream in("HVfile.txt");
  if(!in){cout<<"[!]\n[!] File NOT FOUND!\n[!]"<<endl;
    return;}
  
  std::string b = "null";
  in>>b>>b>>b>>b>>b>>b>>b>>b>>b;
  
  UInt_t lineCounter = 0;
  for(UInt_t i=0; i<ndata; i++){
    Int_t col1, col2, col3, col4, col5, col6, col7, col8, col9;    
    in>>col1>>col2>>col3>>col4>>col5>>col6>>col7>>col8>>col9;
    //cout<<col1<<" "<<col2<<" "<<col3<<" "<<col4<<" "<<col5<<" "<<col6<<" "<<col7<<" "<<col8<<" "<<col9<<endl;
    HVS1[i]      = col2;
    sHVS1[i]     = HVS1[i]*0.01;
    countsS1[i]  = col3;
    sCountsS1[i] = TMath::Sqrt(countsS1[i]);
    HVS2[i]      = col4;
    sHVS2[i]     = HVS2[i]*0.01;
    countsS2[i]  = col5;
    sCountsS2[i] = TMath::Sqrt(countsS2[i]);
    HVS3[i]      = col6;
    sHVS3[i]     = HVS3[i]*0.01;
    countsS3[i]  = col7;
    sCountsS3[i] = TMath::Sqrt(countsS3[i]);
    HVSG[i]      = col8;
    sHVSG[i]     = HVSG[i]*0.01;
    countsSG[i]  = col9;
    sCountsSG[i] = TMath::Sqrt(countsSG[i]);
    lineCounter++;
  }
  
  in.close();


  for(UInt_t i=0; i<ndata; i++){
    rateS1[i] = countsS1[i]/200.;
    rateS2[i] = countsS2[i]/200.;
    rateS3[i] = countsS3[i]/200.;
    rateSG[i] = countsSG[i]/200.;
  }
  
  //  while(in>>inFile){
  //    if(inFile=="break") break;
  //    cout << inFile << endl;


  //double chiquadro95[31] = {0, 3.84, 5.99, 7.81, 9.49, 11.1, 12.6, 14.1, 15.5, 16.9, 18.3, 19.7, 21.0, 22.4, 23.7, 25.0, 26.3, 27.6, 28.9, 301, 31.4, 32.7, 33.9, 35.2, 36.4, 37.7, 38.9, 40.1, 41.3, 42.6, 43.8};
  
  
  
  //TelaB
  TCanvas *canvas = new TCanvas("canvas","N(HV)", 200, 10, 600, 400);
  canvas->Divide(2,2);

  //Scint S1
  canvas->cd(1);
  gPad->SetLogy();
  TGraph *graph1 = new TGraph(ndata, HVS1, countsS1);
  graph1 -> SetMarkerStyle(20);
  graph1 -> SetMarkerColor(kBlue);
  graph1 -> SetMarkerSize(0.7);
  graph1 -> SetTitle("N(HV) scint S1");
  graph1 -> GetXaxis() -> SetTitle("HV (V)");
  graph1 -> GetYaxis() -> SetTitle("N");
  graph1 -> Draw("AP");



  //Grafico CH2 ROSSO
  canvas->cd(2);
  gPad->SetLogy();
  TGraph *graph2 = new TGraph(ndata, HVS2, countsS2);
  graph2 -> SetMarkerStyle(20);
  graph2 -> SetMarkerColor(kOrange);
  graph2 -> SetMarkerSize(0.7);
  graph2 -> SetTitle("N(HV) scint S2");
  graph2 -> GetXaxis() -> SetTitle("HV (V)");
  graph2 -> GetYaxis() -> SetTitle("N");
  graph2 -> Draw("AP");




  //Grafico CH3 ROSSO
  canvas->cd(3);
  gPad->SetLogy();
  TGraph *graph3 = new TGraph(ndata, HVS3, countsS3);
  graph3 -> SetMarkerStyle(20);
  graph3 -> SetMarkerColor(kRed);
  graph3 -> SetMarkerSize(0.7);
  graph3 -> SetTitle("N(HV) scint S3");
  graph3 -> GetXaxis() -> SetTitle("HV (V)");
  graph3 -> GetYaxis() -> SetTitle("N");
  graph3 -> Draw("AP");




  //Grafico CHG ROSSO
  canvas->cd(4);
  gPad->SetLogy();
  TGraph *graphG = new TGraph(ndata, HVSG, countsSG);
  graphG -> SetMarkerStyle(20);
  graphG -> SetMarkerColor(kGreen);
  graphG -> SetMarkerSize(0.7);
  graphG -> SetTitle("N(HV) glass scint");
  graphG -> GetXaxis() -> SetTitle("HV (V)");
  graphG -> GetYaxis() -> SetTitle("N");
  graphG -> Draw("AP");
  











  
    cout << endl;
    cout << "\\begin{table}[h!] \n \\centering \n \\caption{Numero conteggi (N) in funzione della tensione (HV) ai PMT} \n \\label{counts(HV)} \n \\begin{tabular}{cccccccccc} \n \\toprule" << endl;
  cout << "\\# & AQtime (s) &  $\\textup{HV}_{S1}$ (V)  &  $\\textup{N}_{S1}$  &  $\\textup{HV}_{S2}$ (V)  &  $\\textup{N}_{S2}$  &  $\\textup{HV}_{S3}$ (V)  &  $\\textup{N}_{S3}$  &  $\\textup{HV}_{SG}$ (V)  &  $\\textup{N}_{SG}$ \\\\" << endl;
  cout << " & acquisizione (s) \\\\" << endl;
  cout << "\\midrule" << endl;
  for(UInt_t i=0; i<ndata; i++){
    cout << i+1 << " & " << AQtime[i] << " & " << HVS1[i]<<"$\\pm$"<<sHVS1[i] << " & " << countsS1[i]<<"$\\pm$"<<sCountsS1[i] << " & " << HVS2[i]<<"$\\pm$"<<sHVS2[i] << " & " << countsS2[i]<<"$\\pm$"<<sCountsS2[i] << " & " << HVS3[i]<<"$\\pm$"<<sHVS3[i] << " & " << countsS3[i]<<"$\\pm$"<<sCountsS3[i]<< " & " << HVSG[i]<<"$\\pm$"<<sHVSG[i] << " & " << countsSG[i]<<"$\\pm$"<<sCountsS3[i] << " \\\\" << endl;
  }
  cout << "\\bottomrule \n \\end{tabular} \n \\end{table}" << endl;
  cout << endl;

  

cout << "\n\n----------------------------------------------------------------------" << endl;  
  

  

  
  
}
