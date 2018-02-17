#if !defined(__CINT__) || defined(__MAKECINT__)

#include <Riostream.h>
#include "TMath.h"
#include "TF1.h"
#include "TStopwatch.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TAxis.h"                        // ci serve per manipolare gli assi dei grafici
#include "TVirtualFitter.h"               // ci serve come supporto per i fit
#include "string"
#include "TLegend.h"
#include "TLine.h"
using namespace std;

#endif

//Dark Current for PMT X2020
void plot(){
  
  TStopwatch watch;
  watch.Start(kTRUE);

  const UInt_t ndata = 6;
  Float_t  HV[ndata] = {1950., 2000., 2050., 2100., 2150., 2200.}; //[V]
  Float_t sHV[ndata] = {   0.}; //+-1% [V]
  Float_t  DC[ndata] = {   7.,   13.,   20.,   26.,   32.,   38.}; //[mV]

  for(UInt_t i=0; i<ndata; i++){
    sHV[i]  = HV[i]*0.01; //+-1%
  }
  
  TCanvas *canvas = new TCanvas("canvas","Dark Current PMT X2020: average pulse height vs HV", 200, 10, 600, 400);


  //CON ERRORI????????
  /*
  TGraphErrors *graph = new TGraphErrors(ndata, HV, effSG, sHV, sEffSG);
  graph -> SetMarkerStyle(20);
  graph -> SetMarkerColor(kRed);
  graph -> SetMarkerSize(0.7);
  graph -> SetTitle("Dark Current PMT X2020: average pulse height vs HV");
  graph -> GetXaxis() -> SetTitle("HV (V)");
  graph -> GetYaxis() -> SetTitle("pulse height (mV)");
  graph -> Draw("AP");
  */

  TGraph *graph = new TGraph(ndata, HV, DC);
  graph -> SetMarkerStyle(20);
  graph -> SetMarkerColor(kRed);
  graph -> SetMarkerSize(0.7);
  graph -> SetTitle("Dark Current PMT X2020: average pulse height vs HV");
  graph -> GetXaxis() -> SetTitle("HV (V)");
  graph -> GetYaxis() -> SetTitle("pulse height (mV)");
  graph -> Draw("AP");


  cout << endl;
  cout << "\\begin{table}[h!] \n \\centering \n \\caption{Dark Current (DC) PMT X2020: average pulse height (PHE) vs HV} \n \\label{darkCurrent} \n \\begin{tabular}{ccc} \n \\toprule" << endl;
  cout << "\\# & $\\Delta V_{\\textup{SET}}$ (V) & Average DC PHE (mV) \\\\" << endl;
  cout << "\\midrule" << endl;
  for(UInt_t i=0; i<ndata; i++){
    //    cout << i+1 << " & " << HV[i] << "$\\pm$" << sHV[i] << " & " << DC[i] << " \\\\" << endl; //CON ERRORI??????
    cout << i+1 << " & " << HV[i] << " & " << DC[i] << " \\\\" << endl; // SENZA ERRORI
  }
  cout << "\\bottomrule \n \\end{tabular} \n \\end{table}" << endl;
  cout << endl;
  
  
}















  //  while(in>>inFile){
  //    if(inFile=="break") break;
  //    cout << inFile << endl;

  //double chiquadro95[31] = {0, 3.84, 5.99, 7.81, 9.49, 11.1, 12.6, 14.1, 15.5, 16.9, 18.3, 19.7, 21.0, 22.4, 23.7, 25.0, 26.3, 27.6, 28.9, 301, 31.4, 32.7, 33.9, 35.2, 36.4, 37.7, 38.9, 40.1, 41.3, 42.6, 43.8};  
