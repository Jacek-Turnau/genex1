/***********************************************************************
*                        TEventMakerCylindricalPhaseSpace                                 
*                                                                           
* Generate particles and phase space weight according to 
* Cylindical phase space generator given by S. Jadach Comp.Phys.Commun 9(175)297                                                 
*		             
*   									       
***********************************************************************/

#ifndef TEventMakerCylindricalPhaseSpace_H
#define TEventMakerCylindricalPhaseSpace_H


/*!
  \class TEventMakerCylindricalPhaseSpace

  \brief Generate particles and phase space weight according to Cylindical phase space generator given by S. Jadach Comp.Phys.Commun 9(175)297


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
class TEventMakerCylindricalPhaseSpace: public TEventMakerStrategy
{
private:
    
    /// false - if given event is not in acceptance; true otherwise (initally true)
    bool isAccepted;

    /// Pi = 3.14...
    static const double PI;
    
    /// average transverse momentum for generated events [GeV/c]
    double ptt;		
    
    /// dispersion of transverse momenta
    double R;
    
	/// Bad Newton iterations
	int nNewtonbad;
	/// Number of Inf given by matrix elements - to be chceked
	int nInfMik;			
	/// Chcek if generation of event failed
	bool generationFailed;	

	//for Newton iteration:
  
	/// Y value from Newton iteration:
	double yNewt;
	/// Derivative of function for Newton iteration
	double deriv;			
  
	//internal variables:
   
	/// Sum of transverse meomenta
	double sumMt;		
	/// Sum of squares of transverse momenta
	double sumPt2; 		
   
	///configuration file reader
	TPolicyReader * ConfigReader;

	/// Particle Database PDG
	TDatabasePDG * PDGDatabese;

	///reads configuration form the file
	void ReadConfigFile( const string & filename );
    
    /// number of outgoing particles (2 protons + ( nop-2 ) pions) 	       
    int nop;
    
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
    
    /// beam 1 && 2 fourvectors 
    TLorentzVector pb1,pb2;
    
    /// final particles four-vectors
    TLorentzVector * pf;
    
    //arrays to store particles of one event: 
    
    /// store particles energy-momentum vector and rapidity for i-th particle :
    /// pcm[1][i] - longitudal momentum
    /// pcm[2][i], pcm[3][i] - transversal momenta
    /// pcm[4][i] - energy
    /// pcm[5][i] - rapidity
    double * pcm[6];	

    /// store masses of particles corresponding to pcm properties;
    double * mass;
    
    /// store transverse momenta of particles (m^2 + pt^2)
    double * tmass;
  
    /// Logging
	TLog * Logger;
   
public:
       
    ///Constructor
    TEventMakerCylindricalPhaseSpace(void);                    
    
    ///Destructor
    virtual ~TEventMakerCylindricalPhaseSpace(void);           
    
	/// Initialize internal variables of @ref TEventMakerCylindricalPhaseSpace. @note Can be used to erase internal state of the object.
	void Initialize(void);
  
	/// Generate one event in Cylindrical Phase Space starting from 
	/// @param Xarg vector 
	/// @param nDim dimension of @ref Xarg
	///
	/// Adaptation of the algorithm of S. Jadach @ref [3]
	void SetEvent(int nDim, double *Xarg);
  
	///Calculates event and returns Event and returns Event weitgh (double)
	double SetEvent( int nDim, double *Xarg, TEvent * event );
  
	///Check if generation failed
	bool isGenerationFailed( void );


};
#endif
