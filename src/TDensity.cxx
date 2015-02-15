/***********************************************************************
*                       TDesity                             
*                                                                           
* Integrand class for FOAM                                                  
* Inherits  from abstract class TFoamIntegrand.
*									       
*             
*   									       
***********************************************************************/
#include "TDensity.h"




////////////////////////////////////////////////////////////////////////
TDensity::TDensity()
{
	generator = new TEventMaker();
    
	event = new TEvent();
  
	Logger = LoggerSubsystem::Instance();
  
  
  	ConfigReader = ConfigReaderSubsystem::Instance();
	 
	ReadConfigFile( string("Generator.dat") );
  
	//setting up cuts
	if( UseExplorationCuts == 1 )
	{
		weight = new TExplorationAcceptance( new TWeight() );
	}
	else if( UseExplorationCuts == 0 )
	{
		weight = new TGenerationAcceptance( new TWeight() );
	}
  
};

////////////////////////////////////////////////////////////////////////
void TDensity::ReadConfigFile( const string & filename )
{
	//Parse variables: setting up model variables	
	

	UseExplorationCuts  = ConfigReader->GetIntValue( string( "UseExplorationCuts" ) );


	//Log configuration
	string ConfigLog = ConfigReader->GetStringValue( string("ConfigLogFile") );
	TLog * ConfigLogger = new TLog();
	ConfigLogger->setLogTime( false );
	ConfigLogger->openLogFile( ConfigLog );
 
	ConfigLogger->Write() << "# " << "-----------" << endl;
	ConfigLogger->Write() << "# " << "TDensity " << endl;
	ConfigLogger->Write() << "# " << "-----------" << endl;
	ConfigLogger->Write() << "UseExplorationCuts = " << UseExplorationCuts << endl;
	
	delete ConfigLogger;

	
	return;
};

////////////////////////////////////////////////////////////////////////
TDensity::~TDensity()
{
  delete generator;
  delete event;
  delete weight;  
  

  
};

////////////////////////////////////////////////////////////////////////
double TDensity::Density(int nDim, double *Xarg)
{
      
      eventWeight = generator->SetEvent( nDim, Xarg, event );
      
      bool generationFailed = generator->isGenerationFailed();
      
      double fullWeight = 0.0;
      
      if( generationFailed == true )
      {
		
		fullWeight = 0.0;
		
	  }
	  else	
      {
        
		fullWeight = weight->GetWeight( event, eventWeight );
       
      } 
            
      return fullWeight;
};
