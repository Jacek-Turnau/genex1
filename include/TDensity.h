/******************************************************************************
*                         TDesity                             
*                                                                           
* Integrand class for FOAM.                                                  
* Inherits  from abstract class TFoamIntegrand.                            
*									       
* 
*******************************************************************************/
#ifndef TDensity_H
#define TDensity_H


/*!
  \class TDensity

  \brief Integrand class for FOAM. Inherits  from abstract class TFoamIntegrand.
  
  For more information about TFoamIntegrand, see:
  http://root.cern.ch/root/html/TFoamIntegrand.html


 */

#include "TSystem.h"
#include "TROOT.h"
#include "TFoamIntegrand.h"

#include "Global.h"
#include "TEventMaker.h"
#include "TEvent.h"
#include "TWeightStrategy.h"
#include "TWeight.h"
#include "TExplorationAcceptance.h"
#include "TGenerationAcceptance.h"
#include "TConfigReader.h"
#include "TPolicyReader.h"

////////////////////////////////////////////////////////////////////////
class TDensity: public TFoamIntegrand 
{

private:
  
    /// Collect events
	TEvent * event;
	
	/// Event weight
	double eventWeight;
  
    /// TWeight - for calculating weight base on specific matrix element
	TWeightStrategy * weight;
	
	///Use different cuts during generation
	/// 0 - no, 1 - yes
	int UseExplorationCuts;
	
	///configuration file reader
	TConfigReader * ConfigReader;
	
	///Read configuration data form filename
	void ReadConfigFile( const string & filename );
	
	/// Logging
	TLog * Logger;
	
public:
    
    /// keep generator of phase space
    TEventMaker *generator;	      
    
    /// Constructor
    TDensity();                     
    
    /// Destructor
    virtual ~TDensity(void);        

    /// Dendisty function required by TFoam
    /// @returns absolute value of matrix element squared
    double Density(int, double *);      


};
#endif
