/******************************************************************************
*                       TEventMaker2to4                                 
*                                                                           
* Generates phase space for a+b -> a+b + C + D  event.                                             
*		                 
*   									       
*******************************************************************************/

#ifndef TEventMaker2to4_H
#define TEventMaker2to4_H


/*!
  @class TEventMaker2to4

  @brief Generates phase space for a+b -> a+b + C + D  event. It uses modified Lebiedowicz & Szczurek prescription [1], see also [2].

  
  References:
  
  [1] P. Lebiedowicz, A. Szczurek ,,Exclusive \f$2p->2p + \pi^{+} + \pi^{-}\f$ reaction: from the threshold to LHC'', arXiv:0912.0190
  
  [2] L.A. Harland-Lang, V.A. Khoze , M.G. Ryskin, W.J. Stirling "The phenomenology of central exclusive production at hadron colliders'', arXiv:1204.4803 [hep-ph] 
  

 */

#include <stdlib.h>
#include <stdio.h>
#include <cmath>
#include <math.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <assert.h>
#include <complex>

#include <TSystem.h>
#include <TROOT.h>
#include <TRandom.h>
#include <TDatabasePDG.h>
#include <TApplication.h>
#include <TSystem.h>
#include <TFile.h>
#include <TObjArray.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TF1.h>
#include <TAxis.h>
#include <TLine.h>
#include <TMath.h>
#include <TLorentzVector.h>
#include <TPostScript.h>

#include "TParticlePDG.h"
#include "TDatabasePDG.h"

#include "Global.h"
#include "TConfigReader.h"
#include "TPolicyReader.h"
#include "TEvent.h"
#include "TEventMakerStrategy.h"


using namespace std;

////////////////////////////////////////////////////////////////////////
//#define DEBUG
#undef DEBUG

// corecting momenta at low energies
#define CORRECT_MOMENTA_DISPERSION
//#undef CORRECT_MOMENTA_DISPERSION
////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////
class TEventMaker2to4: public TEventMakerStrategy
{
private:
    
    /// Pi = 3.14...
    static const double PI;
   
    /// Chcek if generation of event failed
    bool generationFailed;	
   
	///Determines type of reaction
    int weightStrategy;
   
	///configuration file reader
	TPolicyReader * ConfigReader;

	/// Particle Database PDG
	TDatabasePDG * PDGDatabese;

	///reads configuration form the file
	void ReadConfigFile( const string & filename );
   
	/// store masses of particles corresponding to pcm properties;
    double * mass;
	
	///PDG id of outgouing particles
    int * idOut;   
   
    /// final particles four-vectors
    TLorentzVector * pf;
    
    /// beam 1 && 2 fourvectors 
    TLorentzVector pb1,pb2; 

	///Ranges of generation momentum and pseudorapidity
	double p_min;
	double p_max;
	double y_min;
	double y_max;
	double pmt_max;
   
	///Square of a
    inline double sq(double a){double temp=pow(a,2);return temp;}
   
  	///Frame type
	int frameType;
	///PDG id of particle A
	int idA;
	///PDG id of particle B
	int idB;
	///Energy of particle A in LAB frame
	double EA;
	///Energy of particle B in LAB frame
	double EB;
	///Total energy in CM
	double tecm;
	///Mass of A
	double mA;
	///Mass of B
	double mB;
   
	/// number of outgoing particles (2 protons + ( nop-2 ) pions)   	       
    int nop;
   
	/// Initialize internal variables of @ref TEventMaker2to4. @note Can be used to erase internal state of the object.
    void Initialize(void);
  
	/// CalculateKinematics solves energy-momentum conservation equations
	int CalculateKinematics(double p1t, double p2t, double phi1, double phi2, double y3, double y4, double pmt, double phipmt, int isol);
	
	/// Logging
	TLog * Logger;
   
public:
       
    ///Constructor
    TEventMaker2to4(void);                    
    
    ///Destructor
    virtual ~TEventMaker2to4(void);           
    
	/// Generate one event in Cylindrical Phase Space starting from 
	/// @param Xarg vector 
	/// @param nDim dimension of @ref Xarg
	void SetEvent(int nDim, double *Xarg);
  
  
	/// Generate one event in Cylindrical Phase Space starting from 
	/// @param Xarg vector 
	/// @param nDim dimension of @ref Xarg
	/// @param event event which was generated
	double SetEvent( int nDim, double *Xarg, TEvent * event );
  
	///Check if generation failed
	bool isGenerationFailed( void );
  

};
#endif
