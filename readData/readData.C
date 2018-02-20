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

#endif


Bool_t DBGMODE = kFALSE;

//-- MODULES -----
//
//1. CAEN C257   - scaler             -> 16ch x 32 (24) bit = 512 bits = 64 bytes = 4   righe  --> dal clock
//               - channel 15
//               - slot/station 9
//
//2. LeCroy 2251 - scaler             -> 12ch x 32 (24) bit = 384 bits = 48 bytes = 3   righe  --> non usato 
//
//3. V560N       - scaler             -> 16ch x 32      bit = 512 bits = 64 bytes = 4   righe  --> channel 0
//
//4. V560N       - scaler inhibited   -> 16ch x 32      bit = 512 bits = 64 bytes = 4   righe  --> channel 0
//
//5. 2228A       - tdc                -> 8 ch x 16 (11) bit = 128 bits = 16 bytes = 1   riga   --> channel 7
//               - 
//
//6. 2249A       - adc                -> 12ch x 16 (10) bit = 192 bits = 24 bytes = 1.5 righe  --> non usato
//                                    + 8 bytes padding             =  8 bytes = 0.5 righe 
//
//7. 2249W       - adc                -> 12ch x 16 (11) bit = 192 bits = 24 bytes = 1.5 righe  --> channel 10, 11
//                                    + 8 Bytes padding             =  8 bytes = 0.5 righe
//
//8. V259N       - patter unit        -> 16 bit pattern reg. + 16 bit mul = 4 bytes = 0.25 riga  --> channel 1 (primo bit)
//                                    Non scritto nella descrizione (perché già esauriti i 64 bytes a disposizione),
//                                    ma ci sono sicuramente altri 12 bytes di riempimento = 0.75 riga
//
//9. C211        - programmable delay -> 16ch x  8      bit = 128 bits = 16 bytes = 1   riga
//



//function declaration
int  readBuffer(unsigned char *hex, size_t bufSize, FILE * dataFile);
void skipBuffer(unsigned char *hex, size_t bufSize, FILE * dataFile);



//main
void read(TString inputPath="../../data/night2/pedestal1_20180216_162345.dat", UInt_t numEvtToAnalyse=1000, Bool_t VERB=kFALSE, Bool_t DEBUG=kFALSE){

  //open output file
  TFile *readFile = new TFile("readFile.root", "RECREATE" );

  //open output tree
  TTree *readTree = new TTree("readTree", "Read data tree");
  int validStr[16];
  int lCLOCKch;
  int lSCLinh;   //scaler inhibited
  int lSCLuninh; //scaler not inhibited
  int lTDCch;
  int lADCchS1;
  int lADCchSG;
  int pattReg[16];
  int pattMul[16];
  readTree->Branch("validStr", &validStr, "validStr[16]/I");
  readTree->Branch("CLOCKch" , &lCLOCKch , "lCLOCKch/I"     );
  readTree->Branch("SCLinh"  , &lSCLinh  , "lSCLinh/I"      );
  readTree->Branch("SCLuninh", &lSCLuninh, "lSCLuninh/I"    );
  readTree->Branch("TDCch"   , &lTDCch   , "lTDCch/I"       );
  readTree->Branch("ADCchS1" , &lADCchS1 , "lADCchS1/I"     );
  readTree->Branch("ADCchSG" , &lADCchSG , "lADCchSG/I"     );
  readTree->Branch("pattReg" , &pattReg , "pattReg[16]/I" );
  readTree->Branch("pattMul" , &pattMul , "pattMul[16]/I" );

  DBGMODE = DEBUG;
  printf("[START    ]                           \n");
  printf("[START    ] START data file reading...\n");
  printf("[START    ]                           \n");
  
  //open .dat DATA file
  FILE * dataFile;
  dataFile = fopen(inputPath,"r");
  if (dataFile == NULL){
    cout << "[ERROR] Data file NOT found!!!" << endl;
    return;
  }


  const  UInt_t MAXSIZE      = 150;
  unsigned char hex[MAXSIZE] =  "";

  TH1F *histQ     = new TH1F("histQ   "  , "histQ   "   ,  16,  -0.5,    15.5);
  TH1F *histCLOCK = new TH1F("histCLOCK" , "histCLOCK"  , 110,  -0.5,   110.5);
  TH1F *histADC1  = new TH1F("histADC1"  , "histADC1"   ,2048,  -0.5,  2048.5);
  TH1F *histADCG  = new TH1F("histADCG"  , "histADCG"   ,2048,  -0.5,  2048.5);
  TH1F *histADCsum= new TH1F("histADCsum", "histADCsum" ,2048,  -0.5,  2048.5);
  TH1F *histTIME  = new TH1F("histTIME"  , "histTIME"   , 401,  -0.5,  4000.5);
  
  int nonInhib = 0;
  int inhib    = 0;
  TH1F *histINHIB = new TH1F("histINHIB", "histINHIB"   ,   2,   0. ,     2. );
  TH1F *histDEADT = new TH1F("histDEADT", "histDEADTIME", 101,-0.005,   1.005);
  TH1F *histPATT  = new TH1F("histPATT" , "histPATT"    ,  32,   0.5,    32.5);
  

  UInt_t evtCounter = 0;

  //loop over events
  while(evtCounter < numEvtToAnalyse && !feof(dataFile)) {
    if(VERB){
      printf("[         ] --------------------------\n");
      printf("[READ EVT ] Event %i\n", evtCounter+1);
    }
    
    //Get event size
    int evtSize  = readBuffer(hex, 2, dataFile); //[bytes]
    if(VERB)printf("[READ DATA] Event size       : %d bytes \n", evtSize);
  
    //Get event number
    int evtNum   = readBuffer(hex, 4, dataFile); //[bytes]
    if(VERB)printf("[READ DATA] Event number     : %d       \n",  evtNum);
    if(DBGMODE){
      printf(      "[READ DATA] Event number     : 0x");
      for (int j=0; j<4; j++){printf("%02x", hex[j]);}
      printf("\n");
    }
    
    //Get event validation --> TO DO: skip events not passing event validation
    int evtValid  = readBuffer(hex, 2 , dataFile); //[bytes]
    //int validStr[16];
    for(int i=15; i>=0; i--){
      validStr[i] = evtValid%2;
      evtValid    = evtValid/2;
    }
    if(VERB){
      printf("[READ DATA] Event validation : ");
      for(int i=0; i<16; i++){printf("%i", validStr[i]);}
      printf("\n");
    }
    for(int i=0; i<16; i++){
      histQ->Fill(i,validStr[i]);
    }

    //Skip filling characters
    skipBuffer(hex, 8 , dataFile);
    
    //Skip descriptions for 9 modules --> 9 x 64bytes
    for(int i=0; i<9; i++) skipBuffer(hex, 64, dataFile);
    
    
    
    
    //1. CAEN C257   - scaler        -> 16ch x 32 (24) bit = 512 bits = 64 bytes = 4   righe   - channel 15
    skipBuffer(hex, 15*32/8, dataFile);
    int CLOCKch  = readBuffer(hex, 32/8, dataFile);
    if(VERB)printf("[READ DATA] Scaler CLOCK     : %d periodi da 100 ns\n", CLOCKch);
    histCLOCK -> Fill(CLOCKch);
    lCLOCKch = CLOCKch;
    
    //2. LeCroy 2251 - scaler        -> 12ch x 32 (24) bit = 384 bits = 48 bytes = 3   righe   - NON USATO
    skipBuffer(hex, 12*32/8, dataFile);
    
    //3. V560N       - scaler        -> 16ch x 32      bit = 512 bits = 64 bytes = 4   righe   - channel 0
    int SCLuninh = readBuffer(hex, 32/8, dataFile);
    nonInhib     = SCLuninh;
    skipBuffer(hex, 15*32/8, dataFile);
    if(VERB)printf("[READ DATA] Scaler non inibit: %d \n", SCLuninh);
    lSCLuninh = SCLuninh;
    
    //4. V560N       - scaler inhibited -> 16ch x 32   bit = 512 bits = 64 bytes = 4   righe   - channel 0
    int SCLinh  = readBuffer(hex,  32/8, dataFile);
    inhib       = SCLinh;
    skipBuffer(hex, 15*32/8, dataFile);
    if(VERB)printf("[READ DATA] Scaler inibito   : %d \n", SCLinh);
    lSCLinh = SCLinh;
    
    //5. 2228A       - tdc           -> 8 ch x 16 (11) bit = 128 bits = 16 bytes = 1 riga   - channel 7
    skipBuffer(hex,  7*16/8, dataFile);
    int TDCch   = readBuffer(hex, 16/8, dataFile);
    if(VERB)printf("[READ DATA] TDC channel      : %d \n", TDCch);
    histTIME->Fill(TDCch);
    lTDCch = TDCch;
    
    //6. 2249A       - adc           -> 12ch x 16 (10) bit = 192 bits = 24 bytes = 1.5 righe
    //                                  + 8 bytes padding             =  8 bytes = 0.5 righe - NON USATO
    skipBuffer(hex, 12*16/8 + 8, dataFile);
    
    //7. 2249W       - adc           -> 12ch x 16 (11) bit = 192 bits = 24 bytes = 1.5 righe
    //                                  + 8 Bytes padding             =  8 bytes = 0.5 righe - channels 10, 11
    skipBuffer(hex, 10*16/8, dataFile);
    int ADCchS1 = readBuffer(hex, 16/8, dataFile);
    int ADCchSG = readBuffer(hex, 16/8, dataFile);
    if(VERB)printf("[READ DATA] ADC (input S1)   : %d \n", ADCchS1);
    if(VERB)printf("[READ DATA] ADC (input SG)   : %d \n", ADCchSG);
    skipBuffer(hex, 8, dataFile);
    histADC1  ->Fill(ADCchS1);
    histADCG  ->Fill(ADCchSG);
    histADCsum->Fill(ADCchS1+ADCchSG);
    lADCchS1 = ADCchS1;
    lADCchSG = ADCchSG;
    
    
    //8. V259N       - pattern unit  -> 16 bit pattern reg. + 16 bit mul = 4 bytes = 0.25 riga
    //                               Non scritto nella descrizione (perché già esauriti i 64 bytes a disposizione),
    //                               ma ci sono sicuramente altri 12 bytes di riempimento = 0.75 riga
    int pattRegInt = readBuffer(hex, 16/8, dataFile);
    int pattMulInt = readBuffer(hex, 16/8, dataFile);
    //int pattReg[16];
    //int pattMul[16];
    for(int i=0; i<16; i++){
      pattReg[i] = pattRegInt%2;
      pattRegInt = pattRegInt/2;
      pattMul[i] = pattMulInt%2;
      pattMulInt = pattMulInt/2;
    }
    if(VERB){
      printf("[READ DATA] Pattern register : ");
      for(int i=0; i<16; i++) printf("%i", pattReg[i]);
      printf("\n");
      printf("[READ DATA] Pattern multipl. : ");
      for(int i=0; i<16; i++) printf("%i", pattMul[i]);
      printf("\n");
    }
    skipBuffer(hex, 12, dataFile);
    for(int i=0; i<16; i++){
      histPATT->Fill(i+1, pattReg[i]);
      histPATT->Fill(i+16+1, pattMul[i]);
    }
    
    //9. C211 - programmable delay -> 16ch x 8 bit = 128 bits = 16 bytes = 1 riga - NON USATO
    skipBuffer(hex,  16*8/8, dataFile);
    
    //last line of event
    skipBuffer(hex,      16, dataFile); 
    if(VERB)printf("\n");
    
    readTree->Fill();
    evtCounter++;
  }//end event loop
  
  cout << "[END      ] "                                  << endl;
  cout << "[END      ] End-of-File reached."              << endl;
  cout << "[END      ] Read " << evtCounter << " events." << endl;
  cout << "[END      ] "                                  << endl;
  
  //close data file
  fclose(dataFile);
  
  //dead time
  TVectorF vec(3);   //vector = (nonInhib, inhib, deadTime)
  vec[0] = nonInhib;
  vec[1] = inhib;
  Float_t deadTime = (nonInhib-inhib)/static_cast<Float_t>(nonInhib);
  vec[2] = deadTime;
  histINHIB -> Fill(0.5, nonInhib);
  histINHIB -> Fill(1.5, inhib   );
  histDEADT -> Fill(deadTime     );

  //plots
  TCanvas *canvas  = new TCanvas("canvas", "canvas", 200, 10, 600, 400);
  canvas   -> Divide(3,2);
  canvas   -> cd(1); gPad->SetLogy();
  histTIME -> DrawCopy("hist");

  canvas   -> cd(2); gPad->SetLogy();
  histCLOCK-> DrawCopy("hist");

  canvas   -> cd(3); gPad->Divide(2,1); gPad->cd(1);
  histINHIB-> SetMinimum(0.);
  histINHIB-> DrawCopy("hist");
  canvas   -> cd(3); gPad->cd(2);
  histDEADT-> DrawCopy("hist");

  canvas   -> cd(4);
  histADC1 -> DrawCopy("hist");

  canvas   -> cd(5);
  histADCG -> DrawCopy("hist");

  canvas   -> cd(6); gPad->Divide(2,1); gPad->cd(1);
  histPATT -> DrawCopy("hist");
  canvas   -> cd(6); gPad->cd(2);
  histQ    -> DrawCopy("hist");
  
  canvas   -> Write();
  vec.Write("vec");
  readFile -> Write();
  readFile -> Close();
  
  //histADCsum->DrawCopy("hist");
  
  printf("[END      ] END data file reading...     \n");
  printf("[END      ] OUTPUT in readFile.root \n");
  printf("[END      ]                         \n");
  
  return;
}









int readBuffer(unsigned char *hex, size_t bufSize, FILE * dataFile){
  UInt_t each = 0;
  size_t bytes = 0;
  
  bytes = fread(hex, 1, bufSize, dataFile);
  
  if(DBGMODE){
    cout << "[         ]"                                                   << endl;
    cout << "[DEBUG    ] number of bytes read up to now: " << bytes         << endl;
    if(bytes == bufSize) cout << "[DEBUG    ] correct number of bytes read" << endl;
    printf( "[DEBUG    ] read hexadecimal buffer: 0x");
    for (each = 0; each < bytes; each++) {
      printf( "%02x", hex[each]);
    }
    cout << endl;
    if(bufSize <= 2)      cout <<"[DEBUG    ] 16conversion of the " << bufSize << " bytes to decimal = " <<  be16toh(*(int16_t *)hex) << endl;
    else if(bufSize <= 4) cout <<"[DEBUG    ] 32conversion of the " << bufSize << " bytes to decimal = " <<  be32toh(*(int32_t *)hex) << endl;
    else if(bufSize <= 8) cout <<"[DEBUG    ] 64conversion of the " << bufSize << " bytes to decimal = " <<  be64toh(*(int64_t *)hex) << endl;
    else cout << "Buffer longer than 8 bytes" << endl;
  }  
  
  if     (bufSize <= 2) return be16toh(*(int *)hex);
  else if(bufSize <= 4) return be32toh(*(int *)hex);
  else if(bufSize <= 8) return be64toh(*(int *)hex);
  else{if(DBGMODE) cout << "TOO BIG BUFFER" << endl;
    return -1;}
}



void skipBuffer(unsigned char *hex, size_t bufSize, FILE * dataFile){
  UInt_t each = 0;
  size_t bytes = 0;
  
  bytes = fread(hex, 1, bufSize, dataFile);
  
  if(DBGMODE){
    cout << "[DEBUG    ] number of bytes skipped: " << bytes << endl;
    if(bytes == bufSize) cout << "[DEBUG    ] correct number of bytes skipped" << endl;
    printf( "[DEBUG    ] read hexadecimal buffer: 0x");
    for (each = 0; each < bytes; each++) {
      printf( "%02x", hex[each]);
    }
    cout << endl;
  }  
}
























  /*  
      union
      {
      int result;
      char b[2];
      } U;
      
      {
      U.b[0] = hex[0];
      U.b[1] = hex[1];
      printf ("%04x",U.result);
      uint16_t h_16bit = be16toh(U.result);
      cout << endl;
      printf("0x%X  ->  0x%X\n", U.result, h_16bit);        // 0x1234  ->  0x3412
      }
  */




  
  /*
    size_t lSize = 2;
    size_t result;
    char * buffer;
    
    //allocate memory to contain the read bytes
    buffer = (char*) malloc (sizeof(char)*lSize);
    
    //read bytes and copy them into the buffer
    result = fread (buffer,1,lSize,dataFile);
    
    cout << buffer << endl;



    //free (buffer);

  */
  
