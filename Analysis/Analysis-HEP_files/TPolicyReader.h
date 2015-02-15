/***********************************************************************

					TPolicyReader
  
		Class that reads configuration file with prescribed policy.
Simple class that reads config file into map (key, value).

 File format:
 
 key1 = value1
 key2 = value2
 ...
 
 Lines with comments starts with # or ;



***********************************************************************/

#ifndef TPolicyReader_H
#define TPolicyReader_H

/*!
  @class TPolicyReader

  @brief Class that reads configuration file with prescribed policy.
	
	Simple class that reads config file into map (key, value).

	File format:
 
	key1 = value1
	key2 = value2
	...
 
	Lines with comments starts with # or ;
 
 @warning It applies some policy/constraints on readed data.

 */

#include <iostream>
#include <map>
#include <string>
#include <fstream> 
#include <sstream>
#include <assert.h>
#include <cmath>
#include <math.h>
#include <complex>
#include <string>

#include "TROOT.h"

#include "TConfigReader.h"
#include "TLorentzVector.h"
#include "TParticlePDG.h"
#include "TDatabasePDG.h"


////////////////////////////////////////////////////////////////////////
class TPolicyReader: public TConfigReader
{
public:	
		
	/// Constructor
	TPolicyReader(){};
	
	/// Destructor
	virtual ~TPolicyReader(){};
	
	///Read configuration file
	///@param filename configuration file
	virtual void ReadConfFile( const std::string& filename );	
		
	///Read configuration file without applying policy
	///@param filename configuration file
	virtual void ReadConfFileNoPolicy( const std::string& filename );
		
	///@return boost vector from CM frame to Lab frame if defined or zero vector otherwise
	virtual TVector3 getBoost( void );

private:

	///boost vetor from CM to LAB, and zero vector if LAB=CM
	TVector3  boost;
	
	///Read config file and apply constraints	
	void ApplyPolicy( void );	
	
	///Calculate frame 
	void ClalculateKinematics( void );

	  
};
#endif
