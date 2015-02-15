/******************************************************************************
*                          Class TAcceptance                                
*                                                                           
* Interface class for different acceptances.                                              
*		       	       
*******************************************************************************/

#ifndef TAcceptance_H
#define TAcceptance_H


/*!
  @class TAcceptance

  @brief Interface class for different acceptances.

 */



///////////////////////////////////////////////////////////////////////////////


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


#include "Global.h"
#include "TEvent.h"
#include "TWeightStrategy.h"
#include "TConfigReader.h"
#include "TPolicyReader.h"

using namespace std;

////////////////////////////////////////////////////////////////////////
//#define DEBUG
////////////////////////////////////////////////////////////////////////


class TAcceptance: public TWeightStrategy
{

private:

	///Weight 
	TWeightStrategy * weight;  
	
protected:
  
	///configuration file reader
	TConfigReader * ConfigReader;

	/// Particle Database PDG
	TDatabasePDG * PDGDatabese;

	///reads configuration form the file
	virtual void ReadConfigFile( const string & filename );
	 
public:

    ///Constructor
    TAcceptance( TWeightStrategy * wght );                   
    
    ///Destructor
    virtual ~TAcceptance(void);           
  
	/// @returns 0 or 1 in accordance to the fact that event is accepted by detector or not.
	double getAcceptance( TEvent * event );
  
	/// @returns true if event was accepted or false if not.
	bool getIsAccepted( TEvent * event );
	
	///Calculate matrix element
	/// @param event event 
	complex<double> getMatrixElement( TEvent * event );
 
	/// @returns weight of the event
	double GetWeight( TEvent * event, double eventWeight );


};
#endif
