/***********************************************************************
 
				THistogramStrategy
  
Interface for histogramming classes.



***********************************************************************/

#ifndef THISTOGRAMSTARTEGY_H
#define THISTOGRAMSTARTEGY_H

/*!
  @class THistogramStrategy

  @brief Interface for histogramming classes.


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

////////////////////////////////////////////////////////////////////////
//#define DEBUG
////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////
class THistogramStrategy
{
public:
	
	/// Constructor
	THistogramStrategy(){};
			
	/// Destructor
	virtual ~THistogramStrategy(){};

	/// Collect event
	/// @param event - event that contains an event
	/// @param weight - weight of event
	virtual void Fill( TEvent * event, double weight  ) = 0;
	
	/// Set output file
	/// @param filename set the name of the output file
	virtual void SetOutputFile( const string & filename ) = 0;
	
	///Save histograms into file
	///@param xsection cross-section, e.g., for differential cross-section histograms
	virtual void WriteHistograms(double xsection) = 0;
	
	
};

#endif
