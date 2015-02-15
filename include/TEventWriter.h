/***********************************************************************

				TEventWriter
  
Writes events to files for further investigation.


***********************************************************************/

#ifndef TEventWriter_H
#define TEventWriter_H
/*!
  \class TEventWriter

  \brief Writes events to files for further investigation.


 */
 
#include <iostream>
#include <assert.h>
#include <cmath>
#include <complex>
#include <string>
#include <fstream>
#include <vector>
#include <utility>

#include "TROOT.h"
#include "TSystem.h"
#include "TTree.h"
#include "TBranch.h"
#include "TFile.h"
#include "TLorentzVector.h"

#include "TEvent.h"
#include "TConfigReader.h"
#include "TPolicyReader.h"
#include "LHEF.h"

using namespace std;

////////////////////////////////////////////////////////////////////////
//#define DEBUG
#undef DEBUG
////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////
class TEventWriter
{
public:	
	
	/// Constructor
	TEventWriter( TEvent * event );
		
	///Destructor	
	virtual ~TEventWriter();
    
    /// Saves event
    /// Overwrite existing same tree if overvrite is set true
    void Save( bool overwrite = true );
    
    ///Fill tree with event
    void Fill( void );
    
private:

	///Tree for particles
	TTree * tree;
	
	///The event object to make a boost
	TEvent * eventObject;
	
	///The event pointer
	TEvent * eventPointer;
	
	///boost vetor from CM to LAB, and zero vector if LAB=CM
	TVector3  boost;
	
	///text file to save
	ofstream datafile;
	
	///flag for saving events in root file
	double SaveEventDataRoot;
	
	///root filename
	string rootFilename;
	
	///flag for saving evnets in LHE XML file
	double SaveEventDataLHE_XML;
	
	///LHE XML filename
	string lheFilename;
	
	///flag for saving evnets in LHE XML file
	double SaveEventDataLHE_TXT;
	
	///LHE txt file 
	string lheTxtBackupFile;
	
	///LHEF writer
	LHEF::Writer* writer;
	
	///Energy in CM
	double tecm;
	
	///Determines type of reaction
    int weightStrategy;
	
	///events counter
	long int num;
	
	///configuration file reader
	TPolicyReader * ConfigReader;
	
	///Parse configuration filer
	void ReadConfigFile( const string & filename );
	
	/// Logging
	TLog * Logger;
	
};
#endif
