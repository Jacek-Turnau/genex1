/***********************************************************************
*                          Class TDecayIntegral                           
*                                                                           
* Calculates integral of the phase space for decay M-->1,2,...N 
* for M from [threshol,Mmax] in 100 bins. Integral for each bin
* is stored in table to be used in the interpolation routine. 
*                      
* 
***********************************************************************/

#include "TDecayIntegral.h"




////////////////////////////////////////////////////////////////////////
TDecayIntegral::TDecayIntegral( double Mmax) : _bInitialized( false ), M_max( Mmax )
{ 
	
	rndGen 	 = new TRandom3();
	integral = new TIntegral();
	decay 	 = new TDecay();
	interpolation = new TInterpolationTable();	
	
};

////////////////////////////////////////////////////////////////////////
TDecayIntegral::~TDecayIntegral()
{
	delete rndGen;
	delete integral;
	delete decay;
	delete interpolation;
	
};

////////////////////////////////////////////////////////////////////////
double TDecayIntegral::getDecayIntegral( double M )
{
    double integral = interpolation->Interpolate( M ); 
    
    return integral;
};

////////////////////////////////////////////////////////////////////////
void TDecayIntegral::Init( double * mass, int * PDGcode, int N, int NEvents )
{
	assert( N > 1 );
	
	double M_min = 0.0;
	for( int i = 0; i < N; i++ )
	{
		M_min += mass[i];
	}
	//small shift above threshold
	M_min += numeric_limits<double>::epsilon(); 

	assert( M_max > M_min );
	
	int Nbins = 100;
	
	double x[ Nbins + 1 ];
	double y[ Nbins + 1 ];

	double deltaM = (M_max - M_min)/double(Nbins);

	double M = M_min;
	
	int Nevents = NEvents; //per bin
	
	// loop over bin / mass 
	for( int bins = 0; bins < Nbins + 1; bins++ )
	{
		M = M_min + bins * deltaM;
		
		TLorentzVector pM;
		pM.SetXYZM( 0.0, 0.0, 0.0, M );
	
	
		//integration loop
		for( int i = 1; i < Nevents + 1; i++ )
		{
			// set up random numbers
			for( int i = 0; i < 3*N-4; i++)
			{
				rndQueue.push( rndGen->Uniform() );
			}  
		
			//set up decay
			decay->SetDecay(pM, N, mass);
		
			double wtdecay = decay->Generate( rndQueue );
		
			//fill integral for the bin
			integral->AddWeigt( wtdecay );
		}
		
		//fill tables of nodes
		x[ bins ]  = M;
		y[ bins ]  = integral->GetIntegral();
		
		//Clear integral for another bin
		integral->Clear();
		
	}	
	
	//set up nodes
	interpolation->SetNodes( x, y, Nbins + 1 );
	
		
	_bInitialized = true;
	
};

////////////////////////////////////////////////////////////////////////
void TDecayIntegral::InitError( double * mass, int * PDGcode, int N, double relativeError )
{
	assert( N > 1 );
	
	double M_min = 0.0;
	for( int i = 0; i < N; i++ )
	{
		M_min += mass[i];
	}
	//small shift above threshold
	M_min += numeric_limits<double>::epsilon(); 

	assert( M_max > M_min );
	
	int Nbins = 100;
	
	double x[ Nbins + 1 ];
	double y[ Nbins + 1 ];

	double deltaM = (M_max - M_min)/double(Nbins);

	double M = M_min;
		
	// loop over bin / mass 
	for( int bins = 0; bins < Nbins + 1; bins++ )
	{
		M = M_min + bins * deltaM;
		
		TLorentzVector pM;
		pM.SetXYZM( 0.0, 0.0, 0.0, M );
	
		//relative error
		double error = 1.0;
		
		//loop counter
		int nnn = 0;
	
		//integration loop
		while( error > relativeError )
		{
			// set up random numbers
			for( int i = 0; i < 3*N-4; i++)
			{
				rndQueue.push( rndGen->Uniform() );
			}  
		
			//set up decay
			decay->SetDecay(pM, N, mass);
		
			double wtdecay = decay->Generate( rndQueue );
		
			//fill integral for the bin
			integral->AddWeigt( wtdecay );
			
			//gatherering sufficient statitics for prediction
			if( nnn > 10 )
				error = abs(integral->GetError()/integral->GetIntegral());
			
			nnn++;
		}
			
		//fill tables of nodes
		x[ bins ]  = M;
		y[ bins ]  = integral->GetIntegral();
		
		//Clear integral for another bin
		integral->Clear();
		
	}	
	
	//set up nodes
	interpolation->SetNodes( x, y, Nbins + 1 );
	
		
	_bInitialized = true;
	
};

