/***********************************************************************
* 					Generator.cxx
*                              
*   Generator main class.
*
*                           
************************************************************************/

#include "Generator.h"


////////////////////////////////////////////////////////////////////////
Generator::Generator() : isEventGenerated(false), MCwt( 0.0 ), Rprim( 0.0 ), ErrRprim( 0.0 ), TotDim( 0 ), fMaxWtRej( 1.10 ), numberOfEvents( 0 )
{
  
	Logger = LoggerSubsystem::Instance();
  
	ConfigReader = ConfigReaderSubsystem::Instance();
	ReadConfigFile( string("Generator.dat") );
 
	Logger->openLogFile( LogFile );
	Logger->Write(  ) << "Generator :: Read config file." << endl;  
  
	//create ROOT random numbers generator
  
	PseRan = new TRandom3();      //WARNING: The period of the generator = 2**19937-1, it is big, however finite number, be warned !!!
	PseRan->SetSeed(seed);
  
	if( IntegratorSetup == 2 )
	{
		IntegratorRan = new TRandom3();
		IntegratorRan->SetSeed(seed);
	}
  
  
	if( IntegratorSetup == 1 ) //Use Foam
	{
  
		gSystem->Load("libFoam");
		Density = new TDensity();  // Create integrand function
  
		FoamX    = new TFoam("FoamX");  // Create Simulator
  
		Logger->Write( logINFO ) << "*****  GenEx Generator; Foam version "<< FoamX->GetVersion()<<"    *****"<<endl;
  
		//Setup Foam  
  
		FoamX->SetkDim(        kDim);
		FoamX->SetnCells(      nCells);
		FoamX->SetnSampl(      nSampl);
		FoamX->SetnBin(        nBin);
		FoamX->SetOptRej(      OptRej);
		FoamX->SetOptDrive(    OptDrive);
		FoamX->SetEvPerBin(    EvPerBin);
		FoamX->SetMaxWtRej(    MaxWtRej);
		FoamX->SetChat(        Chat);


		FoamX->SetRho( Density );		// Set density function
		FoamX->SetPseRan( PseRan );   // Set random number generator
  
	};
  
  
	MCvect = new double[TotDim]; // vector generated in the MC run  
  
	//create phase space generator for simulations
	eventMaker = new TEventMaker();
  
	event = new TEvent();
  
	Histogram = new THistogram();  
	FullIntegral = new TIntegral();
 
  
	if( UseExplorationCuts == 1 ) //different cuts for exploration and generation
	{
  
		AcceptanceHistogram = new THistogram();
		AcceptanceHistogram->SetOutputFile( string( "genhistos" ) );
 
		AcceptanceIntegral = new TIntegral();

	} 
  
	Weight = new TGenerationAcceptance( new TWeight() );

	EventWriter = new TEventWriter( event );
 
  
};

////////////////////////////////////////////////////////////////////////
Generator::~Generator()
{
	//save all data
	Finish();	  
  
	//dealocate resources
  
	if( IntegratorSetup == 1 )
	{
		delete FoamX;
	};
  
	delete [] MCvect;
  
	delete Density;
	delete eventMaker;
	delete event;
	delete Histogram;
	delete EventWriter;
	
	delete FullIntegral;
  
	if( UseExplorationCuts == 1 )
	{
		delete AcceptanceHistogram;
		delete AcceptanceIntegral;
	}
  
	delete Weight;

	if( IntegratorSetup == 2 )
	{
		delete IntegratorRan;  
	};
  
  
	Logger->Write( logINFO ) << "Generator :: Deleting generator. END." << endl; 
  
	// delete global subsystems
	LoggerSubsystem::Destroy();
	ConfigReaderSubsystem::Destroy();
 
};

////////////////////////////////////////////////////////////////////////
void Generator::ReadConfigFile( const string & filename )
{
	
	//Reading conf file
	ConfigReader->ReadConfFile(filename);

	
	//Parse variables: setting up model variables	
	
	//option of integration
	IntegratorSetup = ConfigReader->GetIntValue( string( "IntegratorSetup" ) );
	
	UseExplorationCuts  = ConfigReader->GetIntValue( string( "UseExplorationCuts" ) );
	
	seed = ConfigReader->GetIntValue( string( "seed" ) );
	
	
	// Foam Parameters
	//kDim      = ConfigReader->GetIntValue( string( "kDim" ) );			// HYP-CUBICAL subspace												// total dimension
	nCells    = ConfigReader->GetIntValue( string( "nCells" ) );		// Number of Cells
	nSampl    = ConfigReader->GetIntValue( string( "nSampl" ) );    	// Number of MC events per cell in build-up
	nBin      = ConfigReader->GetIntValue( string( "nBin" ) );   		// Number of bins in build-up
	OptRej    = ConfigReader->GetIntValue( string( "OptRej" ) );    	// =0, weighted events;  =1, wt=1 events
	OptDrive  = ConfigReader->GetIntValue( string( "OptDrive" ) );  	// (D=2) Option, type of Drive =0,1,2 for TrueVol,Sigma,WtMax
    EvPerBin  = ConfigReader->GetIntValue( string( "EvPerBin" ) ); 		// Maximum events (equiv.) per bin in buid-up
    MaxWtRej  = ConfigReader->GetDoubleValue( string( "MaxWtRej" ) );   // Maximum wt for rejection, for OptRej=1
	Chat      = ConfigReader->GetIntValue( string( "Chat" ) );   		// Chat level

	LogFile                = ConfigReader->GetStringValue( string( "LogFile" ) );
	backupInterval         = ConfigReader->GetIntValue( string( "backupInterval" ) );
	
	weightStrategy = ConfigReader->GetIntValue( "weightStrategy" );

	kDim = ConfigReader->GetIntValue( string( "kDim" ) );

	TotDim    = kDim;		

	if( OptRej == 0 )
	{
			Logger->Write( logWARNING ) << "Generator :: OptRej = 0. Events are weighted." << endl;  
	}
		
	//Log configuration
	string ConfigLog = ConfigReader->GetStringValue( string("ConfigLogFile") );
	TLog * ConfigLogger = new TLog();
	ConfigLogger->setLogTime( false );
	ConfigLogger->openLogFile( ConfigLog );

	ConfigLogger->Write() << "# " << "-----------" << endl;
	ConfigLogger->Write() << "# " << "Generator " << endl;
	ConfigLogger->Write() << "# " << "-----------" << endl;
	ConfigLogger->Write() << "seed = " << seed << endl;
	ConfigLogger->Write() << "kDim = " << kDim << endl;
	ConfigLogger->Write() << "nCells = " << nCells << endl;
	ConfigLogger->Write() << "nSampl = " << nSampl << endl;
	ConfigLogger->Write() << "nBin = " << nBin << endl;
	ConfigLogger->Write() << "OptRej = " << OptRej << endl;
	ConfigLogger->Write() << "OptDrive = " << OptDrive << endl;
	ConfigLogger->Write() << "EvPerBin = " << EvPerBin << endl;
	ConfigLogger->Write() << "MaxWtRej = " << MaxWtRej << endl;
	ConfigLogger->Write() << "Chat = " << Chat << endl;		
	ConfigLogger->Write() << "LogFile = " << LogFile << endl;
	ConfigLogger->Write() << "backupInterval = " << backupInterval << endl;
	ConfigLogger->Write() << "UseExplorationCuts = " << UseExplorationCuts << endl;
	ConfigLogger->Write() << "IntegratorSetup = " << IntegratorSetup << endl;
	ConfigLogger->Write() << "ConfigLogFile = " << ConfigLog << endl;

	//set up config log as possible Generator.dat
	ConfigLogger->Write() << "readConfigLog = 1" << endl;
	//set up additional data that helps to use log as Generator.dat
	ConfigLogger->Write() << "NumberOfEventsToGenerate = " << ConfigReader->GetIntValue( "NumberOfEventsToGenerate" ) << endl;
	

	delete ConfigLogger;
	
	
	return;
};

////////////////////////////////////////////////////////////////////////
int Generator::Init(void )
{ 
  
  

	if( IntegratorSetup == 1 ) //use Foam  
	{
	  
		Logger->Write( logINFO ) << "Start Foam initialization. " << endl;
	  
		// Make TFoam exploration
		FoamX->Initialize(); 
	
		Logger->Write( logINFO ) << "Foam initialized. " << endl;
		long nCalls=FoamX->GetnCalls( );
		Logger->Write( logINFO ) << "--Initialization done with " << nCalls <<" calls in exploration phase" << endl; 
		Logger->Write( logINFO ) << "--Ready for generation." << endl; 
  
		// Getting norm of the integral
		FoamX->GetIntNorm( Rprim, ErrRprim );

  
	} else if ( IntegratorSetup == 2 )
	{
		//Do nothing
	
	};
	
	FullIntegral->Clear();
  
	if( UseExplorationCuts == 1 )
	{
		AcceptanceIntegral->Clear(); 
	};    
  
	return( 0 );
 
};

////////////////////////////////////////////////////////////////////////

int Generator::Generate(void )
{
  
	double acceptance = 0.0;
  
	bool generationFailed = true;
   
  
	if( IntegratorSetup == 1 ) // generation using FOAM
	{  
 
		// generate MC event
		FoamX->MakeEvent();           
		FoamX->GetMCvect( MCvect );   
		FoamX->GetMCwt( MCwt );
    
		eventMaker->SetEvent( TotDim, MCvect, event );
      
		generationFailed = eventMaker->isGenerationFailed();
    
		if( generationFailed == true )
		{
			acceptance = 0.0;
		}	
		else
		{
			acceptance = Weight->getAcceptance( event );
		}
    
		FullIntegral->AddWeigt( MCwt );
      
		if( UseExplorationCuts == 1 )
		{
			AcceptanceIntegral->AddWeigt( MCwt * acceptance );
		};
          
      
		#ifdef DEBUG        
			cout<<"MC weight = " << Rprim * MCwt <<",  ";
			cout<<"MCvector = ";
			for ( int k=0 ; k < TotDim; k++ ) cout << MCvect[k] << " "; 
			cout<< endl;
		#endif   
      
      
		Histogram->Fill( event, MCwt );
		
		
		EventWriter->Fill();
    
      
		if( UseExplorationCuts == 1 )
		{
			AcceptanceHistogram->Fill( event, MCwt * acceptance );
		};
      
	}else if( IntegratorSetup == 2)  //generation using uniform random numbers
	{
	
		//Fill uniformly random vector
		for( int j = 0; j < TotDim; j++ )
		{
			MCvect[ j ] = IntegratorRan->Uniform();
		}	
		
		//Generate Event 
		double wght = eventMaker->SetEvent( TotDim, MCvect, event );
	
		//Calculate Weight
		generationFailed = eventMaker->isGenerationFailed();
      
		double MCwt = 0.0;
      
		if( generationFailed == true )
		{
			MCwt = 0.0;
		}
		else	
		{        
			MCwt = Weight->GetWeight( event, wght );
			acceptance = Weight->getAcceptance( event );
		} 
	

		FullIntegral->AddWeigt( MCwt );     		

		Histogram->Fill( event, MCwt );

		EventWriter->Fill();
		
     
		if( UseExplorationCuts == 1 )
		{
			AcceptanceIntegral->AddWeigt( acceptance * MCwt );
			AcceptanceHistogram->Fill( event, MCwt * acceptance );
		};	
	}


	//save data on disk every #bakupInterval events
	numberOfEvents++;
	if( numberOfEvents % backupInterval == 0 )
	{
		backup();
			
	}
	  
	//here the weight of mc event is == 1
	MCwt = 1.0; 

	return( 0 );
	
};

////////////////////////////////////////////////////////////////////////
void Generator::backup(void)
{

	EventWriter->Save( );	
	
	Logger->Write( logINFO ) << "Number of steps = " << numberOfEvents << endl;
	Logger->Write( logINFO ) <<  "Events saved" << endl;
	
};


////////////////////////////////////////////////////////////////////////
void Generator::getIntegral( double & integral, double & error )
{
	
	integral = 0.0;
	error    = 0.0;

	if( IntegratorSetup == 1)
	{   
		FoamX->GetIntegMC( integral, error );  
	}
	else if( IntegratorSetup == 2 ) // USe random number generator
	{
		integral = FullIntegral->GetIntegral();
		error    = FullIntegral->GetError();
	}
	
};

////////////////////////////////////////////////////////////////////////
void Generator::getAcceptanceIntegral( double & integral, double & error )
{

	if( UseExplorationCuts != 1 )
	{
		Logger->Write( logERROR ) << "Generator::getAcceptanceIntegral::UseExplorationCuts != 1" << endl;
		throw string( "Generator::getAcceptanceIntegral::UseExplorationCuts != 1" );
	}
	
	
	integral = 0.0;
	error    = 0.0;

	if( IntegratorSetup == 1)
	{   
		double tmpIntegral = 0.0;
		double tmpError    = 0.0;
	
		FoamX->GetIntegMC( tmpIntegral, tmpError ); 	
	
		//Weight = 1 therefore accceptanceIntegral = fullIntegral * frequency of hits inside acceptance region
		integral = tmpIntegral * AcceptanceIntegral->GetIntegral();
		//Binomial/Bernoulli error
		error = tmpError * AcceptanceIntegral->GetIntegral();
		
		//error    = sqrt( AcceptanceIntegral->GetIntegral() *( 1.0 - AcceptanceIntegral->GetIntegral() ) / double( AcceptanceIntegral->GetNSteps() ) );	
	}
	else if( IntegratorSetup == 2 ) // USe random number generator
	{
		Logger->Write( logERROR ) << "Generator::getAcceptanceIntegral:: IntegratorSetup == 2 - Option not supported!" << endl;
		throw string( "Generator::getAcceptanceIntegral:: IntegratorSetup == 2 - Option not supported!" );
	}

};

////////////////////////////////////////////////////////////////////////
void Generator::Finish( void )
{  
     
	TLog integralLog;
	integralLog.openLogFile( string( "integral.log" ) );
	integralLog.setLogTime( false );
    
    // Print Foam statistics on integration
    if( IntegratorSetup == 1)
    {
		Double_t norm = 0.0;
		Double_t rErr = 0.0;
		
		FoamX->Finalize( norm, rErr );
	}


	double integral = 0.0;
	double error    = 0.0;
		
	getIntegral( integral, error );
	
	cout << "INTEGRAL" << endl;	
	cout << "------------------------------------------------" << endl;
	cout << "Total number of events = " <<  numberOfEvents << endl;
	cout << "------------------------------------------------" << endl;
	cout << "Integral   = " << integral << endl;
	cout << "Error    = " << error << endl;
	cout << "Error    = " << error/integral * 100.0 << " % " <<endl;
	cout << "------------------------------------------------" << endl;
		
	integralLog.Write() << "#INTEGRAL" << endl;	
	integralLog.Write() << "------------------------------------------------" << endl;
	integralLog.Write() << "#Total number of events = " <<  numberOfEvents << endl;
	integralLog.Write() << "#------------------------------------------------" << endl;
	integralLog.Write() << "Integral   = " << integral << endl;
	integralLog.Write() << "Error    = " << error << endl;
	integralLog.Write() << "#Error    = " << error/integral * 100.0 << " % " <<endl;
	integralLog.Write() << "#------------------------------------------------" << endl;
			
	//make plots
	Histogram->WriteHistograms( integral );
		
	
	if ( UseExplorationCuts == 1 )
	{
		
		getAcceptanceIntegral( integral, error );
		
		cout << "EXPLORATION PHASE INTEGRAL" << endl;	
		cout << "------------------------------------------------" << endl;
		cout << "Total number of events = " << AcceptanceIntegral->GetNSteps() << endl;
		cout << "------------------------------------------------" << endl;
		cout << "Integral   = " << integral << endl;
		cout << "Error    = " << error << endl;
		cout << "Error    = " << error/integral * 100.0 << " % " <<endl;
		cout << "Effectivenes of generation = " << AcceptanceIntegral->GetIntegral() << endl;
		cout << "------------------------------------------------" << endl;
		
		integralLog.Write() << "#EXPLORATION PHASE INTEGRAL" << endl;	
		integralLog.Write() << "#------------------------------------------------" << endl;
		integralLog.Write() << "#Total number of events = " << AcceptanceIntegral->GetNSteps() << endl;
		integralLog.Write() << "------------------------------------------------" << endl;
		integralLog.Write() << "Integral   = " << integral << endl;
		integralLog.Write() << "Error    = " << error << endl;
		integralLog.Write() << "#Error    = " << error/integral * 100.0 << " % " <<endl;
		integralLog.Write() << "#Effectivenes of generation = " << AcceptanceIntegral->GetIntegral() << endl;
		integralLog.Write() << "------------------------------------------------" << endl;
	
		AcceptanceHistogram->WriteHistograms( integral );
		
	}
  
	EventWriter->Save( );
	
};
