/***********************************************************************
*                       THistogram                                 
*                                                                           
* Selects histogramming class depending on the number of particles in final state.
*        
*		       
*   									       
***********************************************************************/

#include"THistogram.h"



////////////////////////////////////////////////////////////////////////
THistogram::THistogram()
{
  
	ConfigReader = ConfigReaderSubsystem::Instance();
 
 	Logger = LoggerSubsystem::Instance();
  
	ReadConfigFile( string("Generator.dat") );
 
  
	switch( weightStrategy )
	{
		case( 2 ):  strategy = new THistogram4();
					break;
			
		case( 1 ):	strategy = new THistogram4();
					break;				
							
		case( 3 ): strategy = new THistogramN();
					break;
		case( 4 ): strategy = new THistogramN();
					break;
		case( 5 ): strategy = new THistogramf0to2();
					break;
		case( 6 ): strategy = new THistogramf0();
					break;
										
		default:  cout << "THistogram :: NO SUCH STARTEGY" << endl;
				  throw( 0 );
				  break;
	}
  
   
};

////////////////////////////////////////////////////////////////////////
THistogram::~THistogram(void)
{

	delete strategy;
  
};

////////////////////////////////////////////////////////////////////////
void THistogram::ReadConfigFile( const string & filename )
{
	//Parse variables: setting up model variables	

	eventGenerationStrategy = ConfigReader->GetIntValue( "eventGenerationStrategy" );
	weightStrategy = ConfigReader->GetIntValue( "weightStrategy" );
	
	
	//Log configuration
	string ConfigLog = ConfigReader->GetStringValue( string("ConfigLogFile") );
	TLog * ConfigLogger = new TLog();
	ConfigLogger->setLogTime( false );
	ConfigLogger->openLogFile( ConfigLog );

	ConfigLogger->Write() << "# " << "-----------" << endl;
	ConfigLogger->Write() << "# " << "THistogram " << endl;
	ConfigLogger->Write() << "# " << "-----------" << endl;
	ConfigLogger->Write() << "eventGenerationStrategy = " << eventGenerationStrategy << endl;
	ConfigLogger->Write() << "weightStrategy = " << weightStrategy << endl;
	cout << "weightStrategy = " << weightStrategy << endl;
	
	delete ConfigLogger;
	
	return;
};

////////////////////////////////////////////////////////////////////////
void THistogram::SetOutputFile( const string & filename )
{
	
	strategy->SetOutputFile( filename );
	
};

////////////////////////////////////////////////////////////////////////
void THistogram::WriteHistograms(double xsection)
{
	
	strategy->WriteHistograms(xsection);
	
};

////////////////////////////////////////////////////////////////////////
void THistogram::Fill( TEvent * event, double weight  )
{

	strategy->Fill( event, weight );

};
