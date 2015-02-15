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
1) Copy LHE txt file called 'LHETXTEventFile' in Generator.dat file to this folder.
2) Copy Generator.dat file, MODEL_DATA directory and integral.log file which was used to generate root file into this folder.
3) Adjust analysis.cxx to your specific tasks - comments can be a guide

To start program type: make run 
To clean executable type: make clean


NOTE:
It can be used to read DIME files for 2p->2p+2pi reaction if suitably adjusted too:   https://www.hepforge.org/downloads/dimemc



*/


#include <iostream>
#include <assert.h>
#include <cmath>
#include <complex>
#include <string>

#include <iostream>
#include <fstream>
#include <string>

#include "TROOT.h"
#include "TSystem.h"
#include "TApplication.h"
#include "TTree.h"
#include "TBranch.h"
#include "TFile.h"

#include "TLorentzVector.h"


#include <TApplication.h>
#include <TSystem.h>
#include <TFile.h>
#include <TObjArray.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TF1.h>
#include <TAxis.h>
#include <TLine.h>
#include <TMath.h>
#include <TLorentzVector.h>
#include <TPostScript.h>

#include "TConfigReader.h"
#include "TPolicyReader.h"


#undef DEBUG
//#define DEBUG

using namespace std;



int main()
{

	gSystem->Load("libTree.so");
	
	//For stand alone program
	TROOT simple("simple","Test of histogramming and I/O");
	
	
	TPolicyReader * config = new TPolicyReader();
	//read generator config file and aplly policy
	config->ReadConfFile( string("Generator.dat" ) );
	//reading integral and its error without applying policy
	config->ReadConfFileNoPolicy( string("integral.log") );
	//cout << *config << endl;
	
	//show x-section and its error
	cout << "Integral = " << config->GetDoubleValue( string( "Integral" ) ) << endl;
	cout << "Error = " << config->GetDoubleValue( string( "Error" ) ) << endl;
	
	//Nmber of particles
	//const int NIn = 2;
	const int NIn = config->GetIntValue( string( "nip" ) );
	//const int NOut = 4;
	const int NOut = config->GetIntValue( string( "nop" ) );

	//Loretz vectors for particles
	TLorentzVector pb[NIn+1];
	TLorentzVector pf[NOut+1];
	
	
	//Filename with data - TODO read from config file
	string filename = config->GetStringValue( string( "LHETXTEventFile" ) );
	

	//number of events in file
	long int nEvnets = config->GetIntValue( string( "NumberOfEventsToGenerate" ) );

	
	//PROLOG USER SECTION (definitions etc.):

		//example - central invariant mass of pi+ pi- system
		TH1F* cInvMass = new TH1F("cmass","central inv. mass;M [GeV];# events",100,0.0,5.0);	
	
	//END OF PROLOG USER SECTION
  

  
	//ANALYSIS LOOP
	
    ifstream myfile (filename.c_str());
	
	string line;

	stringstream ss; 
			
	double litter;
	double px,py, pz, E;
  
  
	if (myfile.is_open())
	{
		
		
		//for( int ev = 1; ev < 10001; ev ++ )
		for( int ev = 1; ev < nEvnets+1; ev ++ )
		{
			
			
			//This part of the code ectract from HEP format 4-momentum of particles 
			//assuming that first two are mother particles and the rest 4 are daughter particles
			//The code is used as an illustration and example for specific events 2p->2p+2pi and by no means it is portable or flexible - it works.
			//For more elaborated examples please consult :
			//http://arxiv.org/pdf/hep-ph/0609017.pdf
			//http://arxiv.org/pdf/hep-ph/0109068.pdf			
			
			ss << "";
			ss.clear();
			line = "";
			line.clear();
			
			
			//line with numbers
			getline (myfile,line);
			
			#ifdef DEBUG
			cout << line << endl;
			#endif
			
			//////line with pb1
			getline (myfile,line);
			#ifdef DEBUG
			cout << line << endl;
			#endif
		
			ss.str( line ); 
		
			//skip irrelewant data
			for (int i = 1; i < 8; i++ )
			ss >> litter; 
		
			
			ss >> px >> py >> pz >> E;
		
			#ifdef DEBUG
			cout << " pb1x= " << px << " pb1y= " << py << " pb1z= " <<  pz << " pb1E= " << E << endl;
			#endif
		
		
			pb[1].SetPxPyPzE( px, py, pz, E );
		
		
			ss << "";
			ss.clear();
			line = "";
			line.clear();
		
			
			///////line with pb2
			getline (myfile,line);
			#ifdef DEBUG
			cout << line << endl;
			#endif
			
			ss.str( line ); 
			
		
			//skip irrelewant data
			for (int i = 1; i < 8; i++ )
			ss >> litter; 

			ss >> px >> py >> pz >> E;
		
			#ifdef DEBUG
			cout << " pb2x= " << px << " pb2y= " << py << " pb2z= " <<  pz << " pb2E= "<< E << endl;
			#endif
		
		
			pb[2].SetPxPyPzE( px, py, pz, E );
		
			ss << "";
			ss.clear();
			line = "";
			line.clear();

	
			//////line with pp1
			getline (myfile,line);
			#ifdef DEBUG
			cout << line << endl;
			#endif
		
			ss.str( line ); 
		
		
			//skip irrelewant data
			for (int i = 1; i < 8; i++ )
			ss >> litter; 
		

			ss >> px >> py >> pz >> E;
		
			#ifdef DEBUG
			cout << " pp1x= " << px << " pp1y= " << py << " pp1= " <<  pz << " pp1E= "<< E << endl;
			#endif
			
			
			pf[1].SetPxPyPzE( px, py, pz, E );
			
			
			ss << "";
			ss.clear();
			line = "";
			line.clear();
				
			///////line with pp2
			getline (myfile,line);
			#ifdef DEBUG
			cout << line << endl;
			#endif
		
			ss.str( line );
		
			
		
			//skip irrelewant data
			for (int i = 1; i < 8; i++ )
			ss >> litter; 

			ss >> px >> py >> pz >> E;
		
			#ifdef DEBUG
			cout << " pp2x= " << px << " pp2y= " << py << " pp2= " <<  pz << " pp2E= "<< E << endl;
			#endif
			
			
			pf[2].SetPxPyPzE( px, py, pz, E );			
			
			ss << "";
			ss.clear();
			line = "";
			line.clear();
			
			///////line with ppi1
			getline (myfile,line);
			#ifdef DEBUG
			cout << line << endl;
			#endif
			
			ss.str( line ); 
		
	
			//skip irrelewant data
			for (int i = 1; i < 8; i++ )
			ss >> litter; 
		
			ss >> px >> py >> pz >> E;
		
			#ifdef DEBUG
			cout << " ppi1x= " << px << " ppi1y= " << py << " ppi1z= " <<  pz << " ppi1E= "<< E << endl;
			#endif
			
			pf[3].SetPxPyPzE( px, py, pz, E );
			
			
			ss << "";
			ss.clear();
			line = "";
			line.clear();
			
			
			/////line with ppi2
			getline (myfile,line);
			#ifdef DEBUG
			cout << line << endl;
			#endif
		
			ss.str( line );
		
			//skip irrelewant data
			for (int i = 1; i < 8; i++ )
			ss >> litter; 
		
			ss >> px >> py >> pz >> E;
		
			#ifdef DEBUG
			cout << " ppi2x= " << px << " ppi2y= " << py << " ppi2z= " <<  pz << " ppi2E= "<< E << endl;
			#endif
			
			pf[4].SetPxPyPzE( px, py, pz, E );
		
		
			ss << "";
			ss.clear();
			line = "";
			line.clear();
		
		////////////////////////////////////////////////////////////////
		//Write your analysis code here:
		
			//filling histograms - example of filling central invariant mass of pi+ pi-system
			double cmass2 = ( pf[3] + pf[4] ).M2();
			double cmass = sqrt(cmass2);
			cInvMass->Fill(cmass,1.0);
		
		////////////////////////////////////////////////////////////////
		
		}
		myfile.close();
	}

	//END ANALYSIS LOOP
	
	
	//EPILOG USER SECTION (writing files etc.)
	
		string PsFilename("histogram.ps");
		string RootFilename("histogram.root");
	
		//Post Scritp file for histograms
		TPostScript *ps = new TPostScript(PsFilename.c_str(),112);
		// Set style
		gROOT->SetStyle("Plain");
		gROOT->ForceStyle();
		gStyle->SetOptStat(1);
		gStyle->SetTitleBorderSize(0);
		gStyle->SetTitleSize(0.04);
		gStyle->SetTitleFont(42, "hxy");      // for histogram and axis titles
		gStyle->SetLabelFont(42, "xyz");      // for axis labels (values)
		gROOT->ForceStyle();

		//set canvas
		TCanvas* c11 = new TCanvas("c22","Dime");
       //c11->Divide(2,2);
		//set canvas style
		c11->SetFillColor(0);
		c11->UseCurrentStyle();
		c11->SetBorderMode(0);       // still leaves red frame bottom and right
		c11->SetFrameBorderMode(0);   // need this to turn off red hist frame!
		c11->UseCurrentStyle();
  
		ps->NewPage();
		c11->cd(1);
	
		cInvMass->Draw();
	
		c11->Update();
	
		ps->Close();
	
		delete ps;
	
	
		//ROOT file for histograms
		TFile * tf = new TFile( RootFilename.c_str(),"RECREATE", "Histograms form Dime");
     
     
		cInvMass-> Write();
     
		tf->Write();
		tf->Close();
     
		delete tf;

	//END EPILOG SECTION
	
	
	//CLEANING
	delete config;
	
	
return( 0 );	
	
}




