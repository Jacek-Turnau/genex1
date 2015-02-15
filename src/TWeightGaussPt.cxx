/***********************************************************************
*                       TWeightGaussPt                                
*                                                                           
* Returns weight corresponding to gaussian distribution of transverse momenta
* of final state particles.                                                
*		       
*   									       
************************************************************************/


#include"TWeightGaussPt.h"


////////////////////////////////////////////////////////////////////////

const double TWeightGaussPt::PI =  M_PI;


////////////////////////////////////////////////////////////////////////
TWeightGaussPt::TWeightGaussPt()
{
  
	ConfigReader = ConfigReaderSubsystem::Instance();
	
	Logger = LoggerSubsystem::Instance();
	
	PDGDatabese = TDatabasePDG::Instance();
  
	ReadConfigFile( string("Generator.dat") );
  
	//allocate tables:
  
	pf = new TLorentzVector [nop+1];
  
  
  
};

////////////////////////////////////////////////////////////////////////
TWeightGaussPt::~TWeightGaussPt(void)
{
	     
   delete [] pf;

  
};

////////////////////////////////////////////////////////////////////////
void TWeightGaussPt::ReadConfigFile( const string & filename )
{
	
	//Parse variables: setting up model variables	

	nop = ConfigReader->GetIntValue( "nop" );
	alpha =  ConfigReader->GetDoubleValue( "GaussPt::Alpha" );
	norm =  ConfigReader->GetDoubleValue( "GaussPt::Norm" );
	
	//Log configuration
	string ConfigLog = ConfigReader->GetStringValue( string("ConfigLogFile") );
	TLog * ConfigLogger = new TLog();
	ConfigLogger->setLogTime( false );
	ConfigLogger->openLogFile( ConfigLog );

	ConfigLogger->Write() << "# " << "-----------" << endl;
	ConfigLogger->Write() << "# " << "TWeightGaussPt " << endl;
	ConfigLogger->Write() << "# " << "-----------" << endl;
	
	ConfigLogger->Write() << "nop = " << nop << endl;
	ConfigLogger->Write() << "GaussPt::Alpha = " << alpha << endl;
	ConfigLogger->Write() << "GaussPt::Norm = " << norm << endl;
	
	delete ConfigLogger;

	
	return;
};

////////////////////////////////////////////////////////////////////////
double TWeightGaussPt::getAcceptance( TEvent * event )
{
       return( 1.0 );
};

////////////////////////////////////////////////////////////////////////
bool TWeightGaussPt::getIsAccepted( TEvent * event )
{
	bool isAccepted = true;
	
    return isAccepted;
    
};

////////////////////////////////////////////////////////////////////////
complex<double> TWeightGaussPt::getMatrixElement( TEvent * event )
{
	return( 0.0 );	
};

////////////////////////////////////////////////////////////////////////
void TWeightGaussPt::SetEvent( TEvent * event, double eventWeight )
{	
};

////////////////////////////////////////////////////////////////////////
double TWeightGaussPt::GetWeight( TEvent * event, double eventWeight )
{
	
	for( int i = 1; i < nop+1; i++)
		pf[i] = event->pf[i];
	
	double weight = GetWeight( eventWeight );
	
	return( weight );
		
};

////////////////////////////////////////////////////////////////////////
double TWeightGaussPt::GetWeight( double eventWeight )
{
	double wt = eventWeight;
	
	double ptsq = 0.0;
	double pt = 0.0;
	
	for( int i = 1; i < nop + 1; i++ )
	{
		pt = pf[i].Pt();
		ptsq = ptsq + pt * pt;	
	}
	
	wt = wt * exp( -alpha * ptsq );
	
	//Gauss normalization
	double gaussNorm = sqrt( alpha ) / sqrt( M_PI ); 
	wt *= pow( gaussNorm, nop );
	
	wt *= norm;
	
	return wt;
	
};
