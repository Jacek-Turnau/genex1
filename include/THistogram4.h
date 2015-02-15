/***********************************************************************
 
				THistogram4
  
Class for collecting statistics/histograms of reaction a+b->a+b+c+d


***********************************************************************/

#ifndef THistogram4_H
#define THistogram4_H

/*!
  @class THistogram4

  @brief Class for collecting statistics/histograms of reaction a+b->a+b+c+d.


 */


#include <cmath>
#include <math.h>
#include <complex>
#include <string>

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

#include "Global.h"
#include "TConfigReader.h"
#include "TPolicyReader.h"
#include "TEvent.h"
#include "THistogramStrategy.h"

////////////////////////////////////////////////////////////////////////
//#define DEBUG
////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////
class THistogram4 : public THistogramStrategy
{
public:
	
	/// Constructor
	THistogram4();
			
	/// Destructor
	virtual ~THistogram4();

	/// Collect event
	/// @param event - event that contains an event
	/// @param weight - weight of event
	virtual void Fill( TEvent * event, double weight  );
	
	/// Read configuration from the file
	/// @param filename - config file name
	void ReadConfigFile( const string & filename );
	
	/// Set output file
	/// @param filename set the name of the output file
	void SetOutputFile( const string & filename );
	
	///Save histograms into file
	///@param xsection cross-section, e.g., for differential cross-section histograms
	virtual void WriteHistograms( double xsection );
	
private:

	///Number of particles cumulated
	long int Nevents;

	///Weight of event
	double EventWeight;
	
	///Determines type of reaction
    int weightStrategy;
	
	///Number of particles in final state
	int nop;
	
	///Total energy in Center of Mass
	double tecm;
	
	///proton mass [GeV/c^2]
	double mProton;
	
	///pion 	mass [GeV/c^2]  
    double mPion;  
	
	/// beam 1 && 2 fourvectors 
    TLorentzVector  pb1; 
	TLorentzVector  pb2;
	
	/// final particles four-vectors
    TLorentzVector * pf;

    double dazim34,dazim12,coll12,tav,ptprotot,rmass;
    
    double * azim;
    double * theta;
    double * pt;
    double * prap;
    double dz1,dz2;
    
	/// four-body event invariants
	double t1,t2,ta,tb,s13,s24,s14,s23,cmass,cmass2;
	
	///Helper function to extract particles from event to internal class structure
	/// @param event  event that contains an event
	void ExtractParticles( TEvent * event );
	
	///Calculates kinematic variables in 2->4 reaction
	void Kinematics2to4( void );
	
	/// Normalize azimutal angle
	double fdelta(double dazim);

	///Allocate memmory for histograms
	void AllocateHistograms( void );
	
	///Fill histograms with data
	void FillHistograms( void );
	
	///Clear memmory from histograms
	void DeallocateHistograms(void );

	///Helper function which makes differential cross section histograms and write them into file
	/// @param xsection the value of cross section
	void CreateXSectionHistograms( double xsection );
	
	///flag for creation xsection histograms
	bool xsectionHistogramsCreated;
	
	/// generator characteristics histogram array
	TH1F *HOGENE[30];
	TH2F *HOGENE2[30];
	TH2F* hEtaPhi;
	///diffferential cross-section histograms
	TH1F *XHOGENE[30];
	
	///Postscript filename for histograms
	string PsFilename;
	
	///Root file for histograms
	string RootFilename;
	
	///Postscript filename for differential cross section histograms
	string PsXFilename;
	
	///Root file for for differential cross section histograms
	string RootXFilename;
	
	/// Particle Database PDG
	TDatabasePDG * PDGDatabese;
	
	///Reader of configuration file
	TConfigReader * ConfigReader;
	
	/// Logging
	TLog * Logger;

	
};

#endif
