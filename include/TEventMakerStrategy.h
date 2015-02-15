/***********************************************************************
*                        TEventMakerStrategy                                 
*                                                                           
* Interface for different generation strategies for event from the set of random vectors.                                                 
*		                
*   									       
***********************************************************************/

#ifndef TEventMakerStrategy_H
#define TEventMakerStrategy_H


/*!
  \class TEventMakerStrategy

  \brief Interface for different generation strategies for event from the set of random vectors.


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

#include "Global.h"
#include "TEvent.h"


using namespace std;

////////////////////////////////////////////////////////////////////////
//#define DEBUG
#undef DEBUG
////////////////////////////////////////////////////////////////////////



class TEventMakerStrategy
{	
	
 public:
    
    ///Constructor
    TEventMakerStrategy(void){};                    
    
    ///Destructor
    virtual ~TEventMakerStrategy(void){};            
  
	///Calculates event and returns Event and returns Event weitgh (double)
	/// @param Xarg vector 
	/// @param nDim dimension of @ref Xarg
	virtual double SetEvent( int nDim, double *Xarg, TEvent * event ) = 0;
  
    /// If generation will fail returns true and false otherwise
	virtual bool isGenerationFailed( void ) = 0;
  
};
#endif
