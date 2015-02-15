/***********************************************************************
 * 
 * 				Main generator program
 * 
 *//////////////////////////////////////////////////////////////////////

#include <iostream>

using namespace std;

//#define DEBUG_MAIN
#undef DEBUG_MAIN

#include "Generator.h"
#include "TConfigReader.h"
#include "TPolicyReader.h"

///Displays simple progress bar
///\arg [in] Now - how much is left
///\arg [in] Total - general number of events
///\arg [in] width  - how many dots will be used in display
int progressBar( double Now, double Total, int width = 40 );



////////////////////////////MAIN////////////////////////////////////////
int main()
{
  
  ////Create objects:
  //configuration reader
  TConfigReader * configReader = new TPolicyReader();
  //Generator
  Generator * generator = new Generator();
  
  
  cout << "Initializing generator..." << endl;
 
	generator->Init();
  
	configReader->ReadConfFile( string("Generator.dat") );
  
	int TotalEvents = configReader->GetIntValue("NumberOfEventsToGenerate");
  
  cout << "...done" << endl;
  
  cout << "Generate events..." << endl;
 
 
  //loop of generation
  for( int i = 1; i < TotalEvents + 1; i++ )
  {
  
    generator->Generate();
	progressBar( i, TotalEvents );
     
	#ifdef DEBUG_MAIN     
		cout << " event " <<  i << " generated" << endl;
	#endif
    
  };
  
  cout << endl << "... events done" << endl;
  
  //clean up
  delete generator;
  delete configReader;
  
  cout << "EXIT" << endl;
  
  return( 0 );
};

////////////////////////////////////////////////////////////////////////


int progressBar( double Now, double Total, int width )
{
  
    double fraction = Now / Total;
    // number of = to display
    int N = round(fraction * width);

    
    int i=0;
    printf("%3.0f%% [",fraction*100);
    // done part
    for ( ; i < N;i++) 
    {
        printf("=");
    }
    // remaining part (spaces)
    for ( ; i < width;i++) 
    {
        printf(" ");
    }
    //print end and info
    printf("] ");
    printf(" %.0f of %.0f", Now, Total );
    //go back at the beginning
    printf("\r");
    fflush(stdout);
    
    return(0);
    
};
