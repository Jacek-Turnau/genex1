/***********************************************************************
*                        TWeightUnit                                
*                                                                           
* Returns dynamical part of weight corresponding to matrix element = 1.0.                                          
*		                  
*   									       
***********************************************************************/

#ifndef TWeightUnit_H
#define TWeightUnit_H


/*!
  @class TWeightUnit

  @brief Returns dynamical part of weight corresponding to matrix element = 1.0. Can be used for testing.


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


using namespace std;

////////////////////////////////////////////////////////////////////////
//#define DEBUG
////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////
class TWeightUnit: public TWeightStrategy
{
private:
    
    /// Pi = 3.14...
    static const double PI;    
 
	///configuration file reader
	TConfigReader * ConfigReader;

	/// Particle Database PDG
	TDatabasePDG * PDGDatabese;

	///reads configuration form the file
	void ReadConfigFile( const string & filename );
    
    ///Sets event 
	void SetEvent( TEvent * event, double evenWeight = 1.0 );
   
	/// calculate weight based on Cylindrical Phase Space and scattering amplitude
	double GetWeight( double eventWeight = 1.0 );

	/// Logging
	TLog * Logger;
   
public:
    
    ///Constructor
    TWeightUnit(void);                    
    
    ///Destructor
    virtual ~TWeightUnit(void);           
  
	/// @returns 0 or 1 in accordance to the fact that event is accepted by detector or not.
	virtual double getAcceptance( TEvent * event );
  
	/// @returns true if event was accepted or false if not.
	virtual bool getIsAccepted( TEvent * event );
 
	/// @returns weight of the event
	virtual double GetWeight( TEvent * event, double eventWeight );
  
	///Calculate matrix element
	virtual complex<double> getMatrixElement( TEvent * event );

};
#endif
