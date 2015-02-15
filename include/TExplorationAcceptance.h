/***********************************************************************
*                      TExplorationAcceptance                                
*                                                                           
* Decorator class for TAcceptance. It provides acceptance cuts in exploration mode.                                             
*		                  
*   									       
***********************************************************************/

#ifndef TExplorationAcceptance_H
#define TExplorationAcceptance_H


/*!
  @class TExplorationAcceptance

  @brief Decorator class for TAcceptance. It provides acceptance cuts in exploration mode.                                             

  @warning It uses the file "../ACCEPTANCE_CUTS/ExplorationCuts.cxx" for cuts definition.

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

#include "TEvent.h"
#include "TWeightStrategy.h"
#include "TAcceptance.h"


using namespace std;

////////////////////////////////////////////////////////////////////////
//#define DEBUG
////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////
class TExplorationAcceptance: public TAcceptance
{
   
public:

    
    ///Constructor
    TExplorationAcceptance( TWeightStrategy * wght );                   
    
    ///Destructor
    virtual ~TExplorationAcceptance(void);           
  
	/// @returns 0 or 1 in accordance to the fact that event is accepted by detector or not.
	double getAcceptance( TEvent * event );
  
	/// returns true if event was accepted or false if not.
	bool getIsAccepted( TEvent * event );
	
	///Calculate matrix element
	complex<double> getMatrixElement( TEvent * event );
 
	///Returns weight of the event
	double GetWeight( TEvent * event, double eventWeight );


private:
	
	/// Calculate invariants 2->4 body event.
	void Kinematics2to4(TEvent * event);
	
	void ReadConfigFile( const string & filename );
	
	int cutType;
	
	double cut0( void );
	double cut1( void );
	double cut2( void );
	double cut3( void );
	double cut4( void );
	double cut5( void );
	double cut6( void );
	double cut7( void );
	double cut8( void );
	double cut9( void );
	double cut10( void );
	double cut11( void );
	
	
    /// number of outgoing particles (2 protons + ( nop-2 ) pions)
    int nop;
    
    /// total energy in CM 
    double tecm ;		
   
    /// final particles four-vectors
    TLorentzVector * pf;
    
    /// beam 1 && 2 fourvectors 
    TLorentzVector pb1,pb2; 
    
    /// Logging
	TLog * Logger;


};
#endif
