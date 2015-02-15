/***********************************************************************
				TEvent

Stores 4-momenta of generated particles.
***********************************************************************/

#ifndef TEVENT_H
#define TEVENT_H

/*!
  \class TEvent

  \brief Stores 4-momenta of generated particles.


 */


#include <iostream>
#include <vector> 
#include <assert.h>
#include <string>
#include <iomanip>
#include <fstream>

#include "TROOT.h"
#include "TParticlePDG.h"
#include "TDatabasePDG.h"
#include "TLorentzVector.h"
#include "TParticle.h"

#include "Global.h"
#include "TConfigReader.h"
#include "TPolicyReader.h"

using namespace std;

////////////////////////////////////////////////////////////////////////
//#define DEBUG
#undef DEBUG
////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////
class TEvent
{
public:	
	
	/// Constructor
	TEvent();
	
	/// Cosntructor 
	/// \param NIn number of incoming particles
	/// \param NOut number of outgoing particles
	TEvent( unsigned int NIn, unsigned int NOut );
		
	///Destructor	
	virtual ~TEvent();

    /// Get number of In particles
    int GetNIn();
    
    /// Get number of Out particles
    int GetNOut();
    
    /// Beam fourvectors 
    TLorentzVector * pb;
	
	/// Beam particles PDG info
	TParticlePDG ** pbInfo;
	
	/// Final particles four-vectors
    TLorentzVector * pf;
    
    /// Final particles PDG info
    TParticlePDG ** pfInfo;
    
    
    /// List all particles in event
    void PrintParticles();
    
    /// Saves event in file: filename
    void SaveEvent( string filename );
    
    /// Logging
	TLog * Logger;
    
private:

	/// Number of In particles
	int NIn;
	
	/// Number of Out particles
	int NOut;
	
	/// Particle Database PDG
	TDatabasePDG * PDGDatabese;
	
	///configuration file reader
	TConfigReader * ConfigReader;

	///reads configuration form the file
	void ReadConfigFile( const string & filename );
	
};

#endif
