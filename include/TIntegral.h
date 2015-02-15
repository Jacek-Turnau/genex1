/***********************************************************************

				TIntegral
  
  
  Calculates and prints phase space integral and its error 
  estimated in generation phase.


***********************************************************************/

#ifndef TINTEGRAL_H
#define TINTEGRAL_H

/*!
  @class TIntegral

  @brief Calculates and prints phase space integral and its error estimated in generation phase.


 */


#include <math.h>
#include <iostream>
#include <assert.h>

#include "TROOT.h"

#include "Global.h"


using namespace std;

////////////////////////////////////////////////////////////////////////
#define DEBUG
//#undef DEBUG
////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////
class TIntegral
{
public:	
	
	/// Constructor
	TIntegral();
		
	///Destructor	
	virtual ~TIntegral();
	
	///Add weight to MC integral
	void AddWeigt( double weight );
	
	///Returns MC integral
	long double GetIntegral( void );
	
	///Print Integral value
	void PrintIntegral( void );
	
	///Returns MC integral error
	long double GetError( void );
	
	///Returns number of steps remembered
	unsigned long GetNSteps( void );
	
	///Clear all rememebred values
	void Clear( void );
    
private:	
	
	/// Sum of weights
	long double Sf;
	
	/// Sum of squares of weights
	long double Sff;
	
	/// Number of steps performed
	unsigned long NSteps;
	
	/// Logging
	TLog * Logger;
	
};
#endif
