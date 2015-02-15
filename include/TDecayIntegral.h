/***********************************************************************
*                          Class TDecayIntegral                           
*                                                                           
* Calculates integral of the phase space for decay M-->1,2,...N 
* for M from [threshol,Mmax] in 100 bins. Integral for each bin
* is stored in table to be used in the interpolation routine.  
* 
***********************************************************************/

#ifndef TDecayIntegral_H
#define TDecayIntegral_H


/*!
  \class TDecayIntegral

  \brief Calculates integral of the phase space for decay M-->1,2,...N for M from [threshol,Mmax] in 100 bins. Integral for each bin is stored in table to be used in the interpolation routine. 
  
  Typical use:
  1) Create and set the upper bound Mmax for interpolation in constructor
  2) Initialize with Init() or InitError() 
  3) Get all vaues for decaying masses in the interpolation range [threshold; Mmax]


 */

#include <iostream>
#include <queue> 
#include <limits>
#include <assert.h>

#include "TSystem.h"
#include "TROOT.h"
#include "TRandom3.h"
#include <TLorentzVector.h>

#include "TIntegral.h"
#include "TDecay.h"
#include "TInterpolationTable.h"

////////////////////////////////////////////////////////////////////////
using namespace std;
////////////////////////////////////////////////////////////////////////


class TDecayIntegral
{

private:

	/// If decay integral was already calculated
	double _bInitialized;
	
	/// Random number generator
	TRandom * rndGen;
	
	/// Integral 
	TIntegral * integral;
	
	/// decay generator
	TDecay * decay;
	
	///queue for random numbers
    queue<double> rndQueue;
    
    ///Interpolation class
    TInterpolationTable * interpolation;
    
    ///Maximal range of decaying mass
    double M_max;
    
	
public:     
    
    /// Constructor
    /// @param Mmax - maximal value of decaying mass range
    TDecayIntegral( double Mmax = 5.0 );                     
    
    /// Destructor
    virtual ~TDecayIntegral(void);        

    /// @returns decay integral
    double getDecayIntegral( double M );
    
    /// Calculate decay integral
    /// @param mass - mass table
    /// @param PDGcode - PDG codes of particles (not used yet)
    /// @param N - number of particless
    /// @param NEvents - number of events used to integrate decay weight in each bin
    /// @warning Only one initialization is allowed
    void Init( double * mass, int * PDGcode, int N, int NEvents = 10000 );
    
    /// Calculate decay integral - alternative method
    /// @param mass - mass table
    /// @param PDGcode - PDG codes of particles (not used yet)
    /// @param N - number of particless
    /// @param relativeError - bound on relative error (=MCerror/MCvalue) in integrating decay function
    /// @warning relativeError < 10^-3 can lead to very long time of initialization. There is no protection from very long loops!
    /// @warning Only one initialization is allowed
    void InitError( double * mass, int * PDGcode, int N, double relativeError = 0.001 );

	/// @returns initialization flag
	bool isInitialized( void ) {return(_bInitialized); };

};
#endif
