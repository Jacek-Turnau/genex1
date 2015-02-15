/*

SHORT DESCRIPTION:
Simple program to perform analysis on data written by TEventWriter class

LONG DESCRIPTION:
Particles are aviable in tables:
 pb[1..NIn]  - beam particles
 pf[1..NOut] - out particles
 NIn - number of In particles
 NOut - number of Out particles


Analysis can be done in analysis loop. Every turn of the loop is new event extracted from the file.

To perform analysis:
1) Copy root file under the variable 'RootEventFile' in Generator.dat file to this folder.
2) Copy Generator.dat file, MODEL_DATA directory and integral.log file which was used to generate root file into this folder.
3) Adjust analysis.cxx to your specific tasks - comments can be a guide

To start program type: make run 
To clean executable type: make clean


*/


#include <iostream>
#include <assert.h>
#include <cmath>
#include <complex>
#include <string>


#include "TROOT.h"
#include "TSystem.h"
#include "TApplication.h"
#include "TTree.h"
#include "TBranch.h"
#include "TFile.h"

#include "TLorentzVector.h"


#include "TConfigReader.h"
#include "TPolicyReader.h"


using namespace std;


int main()
{

	gSystem->Load("libTree.so");
	
	//For stand alone program
	//TROOT simple("simple","Test of histogramming and I/O");

	TPolicyReader * config = new TPolicyReader();
	//read generator config file and aplly policy
	config->ReadConfFile( string("Generator.dat" ) );
	//reading integral and its error without applying policy
	config->ReadConfFileNoPolicy( string("integral.log") );
	//cout << *config << endl;
	
	//show x-section and its error
	cout << "Integral = " << config->GetDoubleValue( string( "Integral" ) ) << endl;
	cout << "Error = " << config->GetDoubleValue( string( "Error" ) ) << endl;
	//cout << *config << endl;
	
	
	//Nmber of particles
	const int NIn = config->GetIntValue( string( "nip" ) );
	//const int NOut = 4;
	const int NOut = config->GetIntValue( string( "nop" ) );

	//Loretz vectors for particles
	TLorentzVector * pb[NIn+1];
	TLorentzVector * pf[NOut];
	
	//Filename with data
	string filename = config->GetStringValue( string( "RootEventFile" ) );

	TFile file( filename.c_str() );

	// create tree
	TTree *tree = (TTree*)file.Get("events");

	string basename("");
	string number(""); 
	
	
	//in particles
	basename = "In particle ";
  
	for( int i = 1; i < NIn + 1; i++ )
	{
		number = "";		
		number = static_cast<ostringstream*>( &(ostringstream() << i) )->str();
		
		cout << "Set pb[" << i << "]" << endl;
		
		pb[i] = 0;
		
		tree->SetBranchAddress( (basename + number).c_str(), &( pb[i] )  );
		//tree->ResetBranchAddresses();
	}		


	//out particles
	basename = "Out particle ";
  
	for( int i = 1; i < NOut + 1; i++ )
	{
		number = "";		
		number = static_cast<ostringstream*>( &(ostringstream() << i) )->str();
		
		cout << "Set pf[" << i << "]" << endl;
		
		pf[i] = 0;
		
		tree->SetBranchAddress((basename + number).c_str(), &( pf[i] )  );
		//tree->ResetBranchAddresses();
		
	}		

	int n = (int) tree->GetEntries();
	cout << "#Entries " << n << std::endl; 
  
	//PROLOG USER SECTION (definitions etc.):
	
	
	//END OF PROLOG USER SECTION
  
  
	//ANALYSIS LOOP
	for (int i = 0; i < n; ++i) 
	{ 
		
		tree->GetEntry(i);
		
		//DO ANALYSIS:
		
		cout << "Event " << i << endl;
		
		for(int j = 1; j < NIn + 1; j++)
		{ 
			cout << "Retrived pb[" << j <<"] = [" << pb[j]->Px() << ", " << pb[j]->Py() << ", " << pb[j]->Pz() << ", " << pb[j]->E() << "] " << std::endl;
		}
		
		for(int j = 1; j < NOut + 1; j++)
		{ 
			cout << "Retrived pf[" << j <<"] = [" << pf[j]->Px() << ", " << pf[j]->Py() << ", " << pf[j]->Pz() << ", " << pf[j]->E() << "] " << std::endl;
		}
		
		//END ANALYSIS SECTION
	
	}
	//END ANALYSIS LOOP
	
	
	//EPILOG USER SECTION (writing files etc.)
	
	//END EPILOG SECTION
	
	
	//Cleaning
	
	delete tree;
	delete config;
	
return( 0 );	
	
}




