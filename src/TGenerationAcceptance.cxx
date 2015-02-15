/***********************************************************************
*                     TGenerationAcceptance                                
*                                                                           
* Decorator class for TAcceptance. It provides acceptance cuts in generation mode.                                                                                         
*		       
*   									       
***********************************************************************/

#include"TGenerationAcceptance.h"


////////////////////////////////////////////////////////////////////////

//include user defined cuts
#include "../ACCEPTANCE_CUTS/GenerationCuts.cxx"

////////////////////////////////////////////////////////////////////////
TGenerationAcceptance:: TGenerationAcceptance( TWeightStrategy * wght ) : TAcceptance( wght ) 
{
	
	Logger = LoggerSubsystem::Instance();
	
	ReadConfigFile( string("Generator.dat") );
	
	//allocate tables:

	pf = new TLorentzVector [nop+1];
  
};

////////////////////////////////////////////////////////////////////////
TGenerationAcceptance::~TGenerationAcceptance(void)
{

   delete [] pf;
  
};

////////////////////////////////////////////////////////////////////////
void TGenerationAcceptance::ReadConfigFile( const string & filename )
{

	//Parse variables: setting up model variables	
	
	nop = ConfigReader->GetIntValue( "nop" );
	tecm =  ConfigReader->GetDoubleValue( "tecm" );
	cutType = ConfigReader->GetIntValue( "cutType" );

	//Log configuration
	string ConfigLog = ConfigReader->GetStringValue( string("ConfigLogFile") );
	TLog * ConfigLogger = new TLog();
	ConfigLogger->setLogTime( false );
	ConfigLogger->openLogFile( ConfigLog );

	ConfigLogger->Write() << "# " << "-----------" << endl;
	ConfigLogger->Write() << "# " << "TGenerationAcceptance " << endl;
	ConfigLogger->Write() << "# " << "-----------" << endl;
	ConfigLogger->Write() << "nop = " << nop << endl;
	ConfigLogger->Write() << "tecm = " << tecm << endl;
	ConfigLogger->Write() << "cutType = " << cutType << endl;
	
	delete ConfigLogger;
	
	return;
	
};

////////////////////////////////////////////////////////////////////////
double TGenerationAcceptance::getAcceptance( TEvent * event )
{
	
	Kinematics2to4( event );
    
    
    double Acc = 0.0;
  
	switch( cutType )
	{
		case( 0 ): 
			Acc = cut0( );
			break;
		case( 1 ):
			Acc = cut1( );
			break;
		case( 2 ):
			Acc = cut2( );
			break;
		case( 3 ):
			Acc = cut3( );
			break;
		case( 4 ):
			Acc = cut4( );
			break;
		case( 5 ): 
			Acc = cut5( );
			break;		
		case( 6 ): 
			Acc = cut6( );
			break;		
		case( 7 ): 
			Acc = cut7( );
			break;		
		case( 8 ): 
			Acc = cut8( );				
			break;
		case( 9 ): 
			Acc = cut9( );				
			break;
		case( 10 ): 
			Acc = cut10( );				
			break;
		case( 11 ): 
			Acc = cut11( );				
			break;	
		default:
			cout << "NO SUCH CUT" << endl;
			break;
	}
  
      return Acc;
};

////////////////////////////////////////////////////////////////////////
bool TGenerationAcceptance::getIsAccepted( TEvent * event )
{
	
	
	if( getAcceptance( event ) == 0.0 )
		return false;
	else 
		return true;
	
	

};

////////////////////////////////////////////////////////////////////////
void TGenerationAcceptance::Kinematics2to4(TEvent * event)
{
	
	pb1 = event->pb[1];
	pb2 = event->pb[2];
	
	for(int i=1;i<nop+1;i++)
    {
		pf[i] = event->pf[i];
     }

	
	
};

////////////////////////////////////////////////////////////////////////
complex<double> TGenerationAcceptance::getMatrixElement( TEvent * event )
{
	
	if( getAcceptance( event ) == 0.0 )
		return( 0 );	
	else
		return( TAcceptance::getMatrixElement( event ) );
};

////////////////////////////////////////////////////////////////////////
double TGenerationAcceptance::GetWeight( TEvent * event, double eventWeight )
{
	
	if( getAcceptance( event ) == 0.0 )
		return( 0 );	
	else
		return( TAcceptance::GetWeight( event, eventWeight ) );
	
};
