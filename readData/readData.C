#if !defined(__CINT__) || defined(__MAKECINT__)

#include <stdio.h>
#include <endian.h>
#include "TString.h"
#include <Riostream.h>
#include "TStopwatch.h"
#include "TCanvas.h"
#include "TH1.h"

#endif


Bool_t DBGMODE = kFALSE;

//function declaration
int readBuffer(unsigned char *hex, size_t bufSize, FILE * dataFile);


void read(TString inputPath = "../../data/night2/pedestal1_20180216_162345.dat", int numEvtToAnalyse = 1000, Bool_t VERB = kFALSE, Bool_t DEBUG = kFALSE){

  DBGMODE = DEBUG;
  printf("[START    ] \n");
  printf("[START    ] DATA file reading...\n");
  printf("[START    ] \n");
  
  //open .dat DATA file
  FILE * dataFile;
  dataFile = fopen(inputPath,"r");
  if (dataFile == NULL){
    cout << "[ERROR] Data file NOT found!!!" << endl;
    return;
  }


  const UInt_t MAXSIZE = 100;
  unsigned char hex[MAXSIZE] = "";

  TH1I *histADC  = new TH1I("histADC" , "histADC" , 501,  0  , 10000.5);
  TH1F *histTIME = new TH1F("histTIME", "histTIME", 10001, -0.5, 10000.5);
  
  

  //while(fgetc(dataFile) != EOF){
  for(int evtCounter=0; evtCounter<numEvtToAnalyse; evtCounter++){
    
    if(VERB)printf("[READ EVT ] Event %i\n", evtCounter+1);
    
    //Get event size
    int evtSize  = readBuffer(hex, 2, dataFile); //[bytes]
    if(VERB)printf("[READ DATA] Event size       : %d bytes \n", evtSize);
    
    //Get event number
    int evtNum   = readBuffer(hex, 4 , dataFile); //[bytes]
    if(VERB)printf("[READ DATA] Event number     : %d \n", evtNum);
    if(DBGMODE){
    printf("[READ DATA] Event number     : 0x");
    for (int j=0; j<4; j++){
      printf( "%02x", hex[j]);
    }
    cout << endl;}
    
    //Get event validation
    int evtValid = readBuffer(hex, 2 , dataFile); //[bytes]
    int validStr[16];
    for(int i=15; i>=0; i--){
    validStr[i] = evtValid%2;
    evtValid = evtValid/2;
    }
    if(VERB){
      printf("[READ DATA] Event validation : ");
      for(int i=0; i<16; i++){
	printf("%i", validStr[i]);
      }
      printf("\n");
    }
    //Skip filling characters
    int fill1 = readBuffer(hex, 8 , dataFile); //[bytes]
    
    //Skip descriptions for 9 modules --> 9 x 64bytes
    int descr = 0;
    for(int i=0; i<9; i++){
      int descr = readBuffer(hex, 64, dataFile);
    }
    
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
    
    
    
    
    
    //1. CAEN C257   - scaler             -> 16ch x 32 (24) bit = 512 bits = 64 bytes = 4   righe
    //               - channel 15
    int skip1   = readBuffer(hex, 15*32/8, dataFile);
    int count1  = readBuffer(hex,    32/8, dataFile);
    if(VERB)printf("[READ DATA] Count 1          : %d \n", count1);
    
    //2. LeCroy 2251 - scaler             -> 12ch x 32 (24) bit = 384 bits = 48 bytes = 3   righe
    //               - NON USATO
    int skip2   = readBuffer(hex, 12*32/8, dataFile);
    
    //3. V560N       - scaler             -> 16ch x 32      bit = 512 bits = 64 bytes = 4   righe
    //               - channel 0
    int count3  = readBuffer(hex,    32/8, dataFile);
    int skip3   = readBuffer(hex, 15*32/8, dataFile);
    if(VERB)printf("[READ DATA] Count 3          : %d \n", count3);
    
    //4. V560N       - scaler inhibited   -> 16ch x 32      bit = 512 bits = 64 bytes = 4   righe
    //               - channel 0
    int count4  = readBuffer(hex,    32/8, dataFile);
    int skip4   = readBuffer(hex, 15*32/8, dataFile);
    if(VERB)printf("[READ DATA] Count 4          : %d \n", count4);
    
    //5. 2228A       - tdc                -> 8 ch x 16 (11) bit = 128 bits = 16 bytes = 1   riga
    //               - channel 7
    int skip5   = readBuffer(hex,  7*16/8, dataFile);
    int count5  = readBuffer(hex,    16/8, dataFile);
    if(VERB)printf("[READ DATA] Count 5          : %d \n", count5);
    histTIME->Fill(count5);
    
    //6. 2249A       - adc                -> 12ch x 16 (10) bit = 192 bits = 24 bytes = 1.5 righe
    //                                    + 8 bytes padding             =  8 bytes = 0.5 righe
    //               - NON USATO
    int skip6   = readBuffer(hex, 12*16/8 + 8, dataFile);
    
    //7. 2249W       - adc                -> 12ch x 16 (11) bit = 192 bits = 24 bytes = 1.5 righe
    //                                    + 8 Bytes padding             =  8 bytes = 0.5 righe
    //               - channels 10, 11
    int skip7a  = readBuffer(hex, 10*16/8, dataFile);
    int count7a = readBuffer(hex,    16/8, dataFile);
    int count7b = readBuffer(hex,    16/8, dataFile);
    if(VERB)printf("[READ DATA] Count 7a         : %d \n", count7a);
    if(VERB)printf("[READ DATA] Count 7b         : %d \n", count7b);
    int skip7b  = readBuffer(hex,       8, dataFile);
    histADC->Fill(count7a+count7b);
    
    //8. V259N       - patter unit        -> 16 bit pattern reg. + 16 bit mul = 4 bytes = 0.25 riga
    //                                    Non scritto nella descrizione (perché già esauriti i 64 bytes a disposizione),
    //                                    ma ci sono sicuramente altri 12 bytes di riempimento = 0.75 riga
    int count8a = readBuffer(hex,    16/8, dataFile);
    int count8b = readBuffer(hex,    16/8, dataFile);
    int pattRegInt = count8a;
    int pattMultInt = count8a;
    int pattReg[16];
    int pattMult[16];
    for(int i=15; i>=0; i--){
      pattReg[i] = pattRegInt%2;
      pattRegInt = pattRegInt/2;
      pattMult[i] = pattMultInt%2;
      pattMultInt = pattMultInt/2;
    }
    if(VERB){
      printf("[READ DATA] Pattern register : ");
      for(int i=0; i<16; i++){
	printf("%i", pattReg[i]);
      }
      printf("\n");
      printf("[READ DATA] Pattern multipl. : ");
      for(int i=0; i<16; i++){
	printf("%i", pattMult[i]);
      }
    }
    int skip8  = readBuffer(hex,      12, dataFile);
    
    //9. C211        - programmable delay -> 16ch x  8      bit = 128 bits = 16 bytes = 1   riga
    //               - NON USATO
    int skip9  = readBuffer(hex,  16*8/8, dataFile);
    
    //LAST LINE
    int skip10 = readBuffer(hex,      16, dataFile);
    
    if(VERB)printf("\n");
  }

  TCanvas *canvas = new TCanvas("canvas", "canvas", 200, 10, 600, 400);
  canvas->Divide(2);
  canvas->cd(1);
  histTIME->Draw("hist");
  canvas->cd(2);
  histADC ->Draw("hist");
  
  /*
    UInt_t each = 0;
    size_t bytes = 0;
    size_t bufSize = 16;
    
    bytes = fread (hex, 1, bufSize, dataFile);
    
    for (each = 0; each < bytes; each++) {
    printf ( "read this char as int %u and as hex %x\n", hex[each], hex[each]);
    //cout << static_cast<uint16_t>(hex[each]) << endl;0
    }
    
    cout << "[INFO] number of bytes read up to now: " << bytes << endl;
    
    cout << "[    ] conversion of first 2 bytes = size of event = " <<  be16toh(*(int *)hex) << endl;
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
  */
  
  
  
  //1. capire quanto è grande un evento
  //2. fare un ciclo sull'evento
  //. leggere pezzo per pezzo per estrarre le info
  
  

  








  
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
  cout << endl;
  
  
  fclose (dataFile);
  //free (buffer);
  
  return;
}











int readBuffer(unsigned char *hex, size_t bufSize, FILE * dataFile){
  UInt_t each = 0;
  size_t bytes = 0;
  //size_t bufSize = 16;
  
  bytes = fread(hex, 1, bufSize, dataFile);
  /*NO!!!
  for (each = 0; each < bytes; each++) {
    printf ( "read this char as int %u and as hex %02x\n", hex[each], hex[each]);
    //cout << static_cast<uint16_t>(hex[each]) << endl;0
  }
  */
  
  if(DBGMODE){
    cout << "[         ]" << endl;
    cout << "[DEBUG    ] number of bytes read up to now: " << bytes << endl;
    printf( "[DEBUG    ] read hexadecimal buffer: 0x");
    
    for (each = 0; each < bytes; each++) {
      printf( "%02x", hex[each]);
    }
    cout << endl;
    cout <<"[DEBUG    ] conversion of the " << bufSize << " bytes to decimal = " <<  be16toh(*(int *)hex) << endl;
  }  

  if(bufSize <= 2) return be16toh(*(int *)hex);
  if(bufSize <= 4) return be32toh(*(int *)hex);
  if(bufSize <= 8) return be64toh(*(int *)hex);
  else{
    if(DBGMODE) cout << "TOO BIG BUFFER" << endl;
    return -1;
  }
  
}
