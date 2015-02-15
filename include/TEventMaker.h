/***********************************************************************
*                          TEventMaker                                 
*                                                                           
* Selects event generation strategy.                                                 
*		       
*   									       
***********************************************************************/


#ifndef TEventMaker_H
#define TEventMaker_H

/*!
  \class TEventMaker

  \brief Selects event generation strategy.


 */


#include <stdlib.h>
#include <stdio.h>
#include <cmath>
#include <math.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <assert.h>
#include <complex>

#include <TSystem.h>
#include <TROOT.h>
#include <TRandom.h>
#include <TDatabasePDG.h>
#include <TApplication.h>
#include <TSystem.h>
#include <TFile.h>
#include <TObjArray.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TF1.h>
#include <TAxis.h>
#include <TLine.h>
#include <TMath.h>
#include <TLorentzVector.h>
#include <TPostScript.h>

#include "TParticlePDG.h"
#include "TDatabasePDG.h"

#include "Global.h"
#include "TConfigReader.h"
#include "TPolicyReader.h"
#include "TEvent.h"
#include "TEventMakerStrategy.h"
#include "TEventMakerCylindricalPhaseSpace.h"
#include "TEventMaker2to4.h"
#include "TEventMaker2to3.h"
#include "TEventMaker2toN.h"

using namespace std;

////////////////////////////////////////////////////////////////////////
//#define DEBUG
#undef DEBUG
////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////
class TEventMaker: public TEventMakerStrategy
{
private:
  
	///Strategy for generating event from random numbers
	TEventMakerStrategy * strategy;   
   
	///Configuration file reader
	TConfigReader * ConfigReader;

	///Reads configuration form the file
	/// @param filename filename to read configuration from
	void ReadConfigFile( const string & filename );
	
	/// Initialize strategy
	void Initialize(void);
	
	///Type of generation strategy
    int eventGenerationStrategy;
    
    /// Logging
	TLog * Logger;
   
public:
          
    ///Constructor
    TEventMaker(void);                    
    
    ///Destructor
    virtual ~TEventMaker(void);           
  
	///Calculates event and returns Event and returns Event weitgh (double)
	/// @param Xarg vector 
	/// @param nDim dimension of \ref Xarg
	double SetEvent( int nDim, double *Xarg, TEvent * event );
  
	///Check if generation failed
	bool isGenerationFailed( void );

};
#endif
