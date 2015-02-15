/***********************************************************************
*                       TWeight                                 
*                                                                           
* Choice of matrix element (dynamic part of the event weight wgt_ME).
*		            
*   									       
***********************************************************************/

#ifndef TWeight_H
#define TWeight_H

/*!
  @class TWeight

  @brief Choice of matrix element (dynamic part of the event weight wgt_ME).


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
#include "TWeightStrategy.h"
#include "TWeightGaussPt.h"
#include "TWeightUnit.h"
#include "TWeightLSpipi.h"
#include "TWeightLSKK.h"
#include "TWeightf0to2.h"
#include "TWeightf0.h"

using namespace std;

////////////////////////////////////////////////////////////////////////
//#define DEBUG
////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////
class TWeight:  public TWeightStrategy 
{
private:
     
     ///weight
     TWeightStrategy * strategy;
       
	///configuration file reader
	TConfigReader * ConfigReader;

	///reads configuration form the file
	void ReadConfigFile( const string & filename );

	///Model used
	int weightStrategy;
	
	/// Logging
	TLog * Logger;
   
public:
    
    ///Constructor
    TWeight(void);                    
    
    ///Destructor
    virtual ~TWeight(void);           
  
	/// @returns 0 or 1 in accordance to the fact that event is accepted by detector or not.
	virtual double getAcceptance( TEvent * event );
  
	/// returns true if event was accepted or false if not.
	virtual bool getIsAccepted( TEvent * event );
 
	///Returns weight of the event
	virtual double GetWeight( TEvent * event, double eventWeight );
  
	///Calculate matrix element
	virtual complex<double> getMatrixElement( TEvent * event );

};
#endif
