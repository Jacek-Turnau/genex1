/***********************************************************************
			
				TIntegral
  
  Calculates and prints phase space integral and its error 
  estimated in generation phase.

***********************************************************************/


#include "TIntegral.h"

////////////////////////////////////////////////////////////////////////
TIntegral::TIntegral() : Sf( 0.0 ), Sff( 0.0 ), NSteps( 0 )
{
	
	Logger = LoggerSubsystem::Instance();
	
	return;
};


////////////////////////////////////////////////////////////////////////
TIntegral::~TIntegral()
{
	return;
};

////////////////////////////////////////////////////////////////////////
void TIntegral::AddWeigt( double weight )
{
	Sf  = Sf + weight;
	Sff = Sff + weight * weight;
	NSteps++; 
	
	return;
	
};

////////////////////////////////////////////////////////////////////////
long double TIntegral::GetIntegral( void )
{
	assert( NSteps > 0 );
	
	long double tmp =  Sf / static_cast<long double>( NSteps );
	
	return( tmp );
};

////////////////////////////////////////////////////////////////////////
void TIntegral::PrintIntegral( void )
{
	
	
	if( NSteps > 0)
	{
	
		long double integral = GetIntegral();
		long double error    = GetError();
	
		//print results
		cout << "------------------------------------------------" << endl;
		cout << "Total number of events = " << NSteps << endl;
		cout << "------------------------------------------------" << endl;
		cout << "Integral   = " << integral << endl;
		cout << "Error    = " << error << endl;
		cout << "Error    = " << error/integral * 100.0 << " % " <<endl;
		cout << "------------------------------------------------" << endl;
	}			
	else
	{
		cout << "------------------------------------------------" << endl;
		cout << "Total number of events = " << NSteps << endl;
		cout << "------------------------------------------------" << endl;
		cout << "------------------------------------------------" << endl;
		cout << "INTEGRAL NOT DETERMINED" << endl;
		cout << "------------------------------------------------" << endl;
	}
	
	return;
};

////////////////////////////////////////////////////////////////////////
long double TIntegral::GetError( void )
{
	assert( NSteps > 0 );
	
	long double N = static_cast<long double>( NSteps );
	
	long double tmp =  sqrt( Sff - (Sf * Sf) / N )/N;
	
	return( tmp );
};

////////////////////////////////////////////////////////////////////////
unsigned long TIntegral::GetNSteps( void )
{
	assert( NSteps >= 0 );
	
	return( NSteps );	
		
};

////////////////////////////////////////////////////////////////////////
void TIntegral::Clear( void )
{
	Sf = 0.0;
	Sff = 0.0;
	NSteps = 0;
	
	return;	
};
