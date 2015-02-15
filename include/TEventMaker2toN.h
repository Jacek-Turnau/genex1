/***********************************************************************
*                 TEventMaker2toN                                 
*                                                                           
* Generate event and weight of phase space for process 2 -> 2 + M, M->3+4+...+N                                            
*
* It bases on TGenPhaseSpace, which bases on GENBOD by F. James 
* 		                
*   									       
***********************************************************************/

#ifndef TEventMaker2toN_H
#define TEventMaker2toN_H


/*!
  @class TEventMaker2toN

  @brief Generate event and weight of phase space for process 2 -> 2 + M, M->3+4+...+N                                            
		
		It bases on TGenPhaseSpace, which bases on GENBOD by F. James 


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
#include <queue>          

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
#include "TDecay.h"
#include "TDecayIntegral.h"

using namespace std;

////////////////////////////////////////////////////////////////////////
//#define DEBUG
#undef DEBUG
////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////
class TEventMaker2toN: public TEventMakerStrategy
{

private:
    
    /// Pi = 3.14...
    static const double PI;
   
	///configuration file reader
	TPolicyReader * ConfigReader;

	/// Particle Database PDG
	TDatabasePDG * PDGDatabese;

	///reads configuration form the file
	void ReadConfigFile( const string & filename );
   
    /// CalculateKinematics solves energy-momentum conservation equations
	int CalculateKinematics(double p1t, double p2t, double dphi1, double dphi2, double y, double vmass2, int isol);
	
	///calculate 2-particle momentum in their cma frame
	double kcms(double s, double m1, double m2);
   
	///Check if generation failed
    bool generationFailed;
    
    ///Determines type of reaction
    int weightStrategy;
   
    /// number of all particles
    int nopf;  
    
    /// final particles four-vectors
    TLorentzVector * pf;
    
    ///central mass four-vector
    TLorentzVector PM;
    
    /// beam 1 && 2 fourvectors 
    TLorentzVector pb1,pb2;
    
        /// number of outgoing particles (2 protons + ( nop-2 ) pions)   	       
    int nop;
    /// store masses of particles corresponding to pcm properties;
    double * mass;
    
    ///PDG id of outgouing particles
    int * idOut;   
    
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
	///4-momentum of A
	TLorentzVector pA;
	///4-momentum of B
	TLorentzVector pB;	
    
    double p_min;
    double p_max;
    double y_min;
    double y_max;
    
    int centralMassOption;
	double mass_min,mass2_min,mass_max, mass2_max;
	double dmom;
	
	///weight for decay
	double wtdecay;
	
	///weight
	double weight;
	
    inline double sq(double a){double temp=pow(a,2);return temp;}
    
    ///queue for random numbers
    queue<double> rndQueue;
    
    ///decay central mass
    TDecay * Decay;
    
    /// Logging
	TLog * Logger;
	
	/// Decay weight
	TDecayIntegral * DecayIntegral;
	
	///initialize TDecayIntegral
	int bIDecay;
    
public:
 
    ///Constructor
    TEventMaker2toN(void);                    
    
    ///Destructor
    virtual ~TEventMaker2toN(void);           
         
	/// Initialize internal variables of @ref TEventMaker2toN. @note Can be used to erase internal state of the object.
	void Initialize(void);
  
	/// Generate one event in Cylindrical Phase Space starting from 
	/// @param Xarg vector 
	/// @param nDim dimension of @ref Xarg
	void SetEvent(int nDim, double *Xarg);
  
	///Calculates event and returns Event and returns Event weitgh (double)
	double SetEvent( int nDim, double *Xarg, TEvent * event );

	/// @returns weight of the last event
	double getWeight( void );
	
	/// @returns weight of the decay
	double getDecayWeight( void );
	
	/// @returns weight of the phase space
	double getDecayPhaseSpaceIntegral( double M );	
  
	bool isGenerationFailed( void );
  

};
#endif
