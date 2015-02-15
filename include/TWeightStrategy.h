/***********************************************************************
*                       TWeightStrategyStrategy                                 
*                                                                                  
* Interface class for different Weight functions.                                        
*		                  
*   									       
***********************************************************************/

#ifndef TWeightStrategy_H
#define TWeightStrategy_H


/*!
  @class TWeightStrategy

  @brief Interface class for different weight functions.

 */

#include <stdlib.h>
#include <stdio.h>
#include <cmath>
#include <math.h>
#include <complex>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <assert.h>

#include <TSystem.h>
#include <TROOT.h>
#include <TRandom.h>
#include <TApplication.h>
#include <TSystem.h>
#include <TFile.h>
#include <TObjArray.h>
#include <TLorentzVector.h>
#include <TParticlePDG.h>
#include <TDatabasePDG.h>

#include "Global.h"
#include "TConfigReader.h"
#include "TPolicyReader.h"
#include "TEvent.h"


using namespace std;


////////////////////////////////////////////////////////////////////////
class TWeightStrategy
{
public:    
    
    ///Constructor
    TWeightStrategy(void){};                    
    
    ///Destructor
    virtual ~TWeightStrategy(void){};           
  
	/// @returns 0 or 1 in accordance to the fact that event is accepted by detector or not.
	virtual double getAcceptance( TEvent * event ) = 0;
  
	/// @returns true if event was accepted or false if not.
	virtual bool getIsAccepted(TEvent * event ) = 0;
 
	/// @returns weight of the event
	virtual double GetWeight( TEvent * event, double eventWeight ) = 0;
  
	///Calculate matrix element
	virtual complex<double> getMatrixElement( TEvent * event ) = 0;
   
};
#endif
