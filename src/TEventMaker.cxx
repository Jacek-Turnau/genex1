/***********************************************************************
*                          TEventMaker                                 
*                                                                           
* Selects event generation strategy.
*		       
*   									       
***********************************************************************/

#include"TEventMaker.h"


////////////////////////////////////////////////////////////////////////
TEventMaker::TEventMaker() : eventGenerationStrategy(-1)
{
  
	Logger = LoggerSubsystem::Instance();
  
	ConfigReader = ConfigReaderSubsystem::Instance();
  
	ReadConfigFile( string("Generator.dat") );
   
	//initialize internal structure of the class
	Initialize();
  
}

////////////////////////////////////////////////////////////////////////
TEventMaker::~TEventMaker(void)
{
	
	delete strategy;
  
};

////////////////////////////////////////////////////////////////////////
void TEventMaker::ReadConfigFile( const string & filename )
{

	//Parse variables: setting up model variables	
	
	eventGenerationStrategy = ConfigReader->GetIntValue( "eventGenerationStrategy" );
	
	
	//Log configuration
	string ConfigLog = ConfigReader->GetStringValue( string("ConfigLogFile") );
	TLog * ConfigLogger = new TLog();
	ConfigLogger->setLogTime( false );
	ConfigLogger->openLogFile( ConfigLog );
 
	ConfigLogger->Write() << "# " << "-----------" << endl;
	ConfigLogger->Write() << "# " << "TEventMaker " << endl;
	ConfigLogger->Write() << "# " << "-----------" << endl;
	ConfigLogger->Write() << "eventGenerationStrategy = " << eventGenerationStrategy << endl;
	
	delete ConfigLogger;
	
	
	return;
};


////////////////////////////////////////////////////////////////////////
void TEventMaker::Initialize()
{ 
	
	switch( eventGenerationStrategy )
	{		
		case( 1 ):	strategy = new TEventMakerCylindricalPhaseSpace();
					break;
		
		case( 2 ):  strategy = new TEventMaker2to4();
					break;
							
		case( 3 ):  strategy = new TEventMaker2to3();
					break;		
					
		case( 4 ): strategy = new TEventMaker2toN();
					break;					
					
		default:	Logger->Write( logERROR ) << "TEventMaker :: NO STARTEGY - STOP" << endl;
					cout << "TEventMaker :: NO STARTEGY - STOP" << endl;
					throw( 0 ); 
					break;
	}		
			
};

////////////////////////////////////////////////////////////////////////
double TEventMaker::SetEvent( int nDim, double *Xarg, TEvent * event )
{	
	return( strategy->SetEvent( nDim, Xarg, event) );
	
};

////////////////////////////////////////////////////////////////////////
bool TEventMaker::isGenerationFailed( void )
{
	
	return( strategy->isGenerationFailed( ) );	
		
};
