#if !defined(__CINT__) || defined(__MAKECINT__)

#include <stdio.h>
#include <endian.h>
#include "TString.h"
#include <Riostream.h>
#include "TStopwatch.h"

#endif



//function declaration
int readBuffer(unsigned char *hex, size_t bufSize, FILE * dataFile);


void read(TString inputPath = "../../data/night2/pedestal1_20180216_162345.dat"){
  
  //open .dat DATA file
  FILE * dataFile;
  dataFile = fopen(inputPath,"r");
  if (dataFile == NULL){
    cout << "[ERROR] Data file NOT found!!!" << endl;
    return;
  }


  const UInt_t MAXSIZE = 100;
  unsigned char hex[MAXSIZE] = "";
  
  int prova = readBuffer(hex, 2, dataFile);
  //cout << prova << endl;


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
    
    cout <<"[] conversion of first 2 bytes = size of event = " <<  be16toh(*(int *)hex) << endl;
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
  /*
  for (each = 0; each < bytes; each++) {
    printf ( "read this char as int %u and as hex %02x\n", hex[each], hex[each]);
    //cout << static_cast<uint16_t>(hex[each]) << endl;0
  }
  */
  
  cout << "[INFO] number of bytes read up to now: " << bytes << endl;
  printf( "[INFO] read hexadecimal buffer: ");
  for (each = 0; each < bytes; each++) {
    printf( "%02x", hex[each]);
  }
  cout << endl;
  cout <<"[INFO] conversion of the " << bufSize << " bytes = " <<  be16toh(*(int *)hex) << endl;
  return be16toh(*(int *)hex);
  //return 0.;
}
