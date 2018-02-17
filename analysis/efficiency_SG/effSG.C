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


void plot(TString inputFileName = "effSG.txt"){
  
  TStopwatch watch;
  watch.Start(kTRUE);

  ifstream in("effSG.txt");
  if(!in){cout<<"[!]\n[!] File NOT FOUND!\n[!]"<<endl;
    return;}
  
  string      b = "null";
  UInt_t  ndata =  0;
  Float_t  time =  0.;
  
  in >> b >> b >> ndata;
  in >> b >> time; //[s]
  cout << "---------------------------------------------" << endl;
  cout << "[INFO] Acquisition time = (" << time*1000 << "+-1) ms" << endl;
  in>>b>>b>>b>>b;

  const UInt_t maxndata = 50;
  Float_t   HVSG[maxndata] = {0.};
  Float_t  sHVSG[maxndata] = {0.};
  Int_t     Ndop[maxndata] = {0};
  Float_t  sNdop[maxndata] = {0.};
  Int_t     Ntri[maxndata] = {0};
  Float_t  sNtri[maxndata] = {0.};
  Float_t   Rdop[maxndata] = {0.};
  Float_t   Rtri[maxndata] = {0.};
  Float_t  effSG[maxndata] = {0.};
  Float_t sEffSG[maxndata] = {0.};

  for(UInt_t i=0; i<ndata; i++){
    Int_t col1, col2, col3, col4;
    in>>col1>>col2>>col3>>col4;
    HVSG[i]   = col2;
    sHVSG[i]  = HVSG[i]*0.01; //+-1%
    Ndop[i]   = col3;
    sNdop[i]  = TMath::Sqrt(Ndop[i]);
    Ntri[i]   = col4;
    sNtri[i]  = TMath::Sqrt(Ntri[i]);
    Rdop[i]   = Ndop[i]/time;
    Rtri[i]   = Ntri[i]/time;
    effSG[i]  = Rtri[i]/Rdop[i];
    sEffSG[i] = TMath::Sqrt(effSG[i]*(1-effSG[i])/(static_cast<Float_t>(Ndop[i])));
  }
  
  in.close();
  
  TCanvas *canvas = new TCanvas("canvas","Efficiency (#varepsilon) of Glass Scintillator (SG)", 200, 10, 600, 400);

  //Scint SG
  TGraphErrors *graphE = new TGraphErrors(ndata, HVSG, effSG, sHVSG, sEffSG);
  graphE -> SetMarkerStyle(20);
  graphE -> SetMarkerColor(kBlue);
  graphE -> SetMarkerSize(0.7);
  graphE -> SetTitle("Glass Scintillator Efficiency (#varepsilon) vs HV");
  graphE -> GetXaxis() -> SetTitle("HV (V)");
  graphE -> GetYaxis() -> SetTitle("#varepsilon");
  graphE -> Draw("AP");



  cout << endl;
  cout << "\\begin{table}[h!] \n \\centering \n \\caption{Numero coincidenze (N) in funzione della tensione (HV) al Vetro Scintillante} \n \\label{counts(HV)} \n \\begin{tabular}{cccc} \n \\toprule" << endl;
  cout << "\\# &  $\\textup{HV}_{\\textup{mon,}SG}$ (V)  &  $\\textup{N}_{\\textup{doppie}}$  &  $\\textup{N}_{\\textup{triple}}$ \\\\" << endl;
  cout << "\\midrule" << endl;
  for(UInt_t i=0; i<ndata; i++){
    cout << i+1 << " & " << HVSG[i] << "$\\pm$" << sHVSG[i] << " & " <<  Ndop[i]<<"$\\pm$"<<sNdop[i] << " & " << Ntri[i]<<"$\\pm$"<<sNtri[i] << " \\\\" << endl;
  }
  cout << "\\bottomrule \n \\end{tabular} \n \\end{table}" << endl;
  cout << endl;
  
  
}















  //  while(in>>inFile){
  //    if(inFile=="break") break;
  //    cout << inFile << endl;

  //double chiquadro95[31] = {0, 3.84, 5.99, 7.81, 9.49, 11.1, 12.6, 14.1, 15.5, 16.9, 18.3, 19.7, 21.0, 22.4, 23.7, 25.0, 26.3, 27.6, 28.9, 301, 31.4, 32.7, 33.9, 35.2, 36.4, 37.7, 38.9, 40.1, 41.3, 42.6, 43.8};  
