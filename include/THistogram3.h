/***********************************************************************
 
				THistogram3
  
Class for Collecting statistics/histograms of reaction a+b->a+b+c.


***********************************************************************/

#ifndef THistogram3_H
#define THistogram3_H

/*!
  @class THistogram3

  @brief Class for Collecting statistics/histograms of reaction a+b->a+b+c.


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
#include "TEvent.h"
#include "THistogramStrategy.h"

////////////////////////////////////////////////////////////////////////
//#define DEBUG
////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////
class THistogram3 : public THistogramStrategy
{
public:
	
	/// Constructor
	THistogram3();
			
	/// Destructor
	virtual ~THistogram3();

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
	virtual void WriteHistograms(double xsection);
	
		
private:

	///Weight of event
	double EventWeight;
	
	/// generator characteristics histogram array
	TH2F* hEtaPhi;
	TH2F* hEtaPhiPhoton;
	TH2F* hEtaPhiPP;
	
	///Postscript filename for histograms
	string PsFilename;
	
	///Root file for histograms
	string RootFilename;
	
	/// Particle Database PDG
	TDatabasePDG * PDGDatabese;
	
	///Reader of configuration file
	TConfigReader * ConfigReader;
	
	///Allocate histograms
	void AllocateHistograms( void );
	
	///Deallocate histograms
	void DeallocateHisograms( void );
	
	/// Logging
	TLog * Logger;
	
};

#endif
