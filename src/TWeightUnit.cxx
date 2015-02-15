/***********************************************************************
*                        TWeightUnit                       
*                                                                           
* Returns dynamical part of weight corresponding to matrix element = 1.0.     
*		       
*   									       
***********************************************************************/

#include"TWeightUnit.h"


////////////////////////////////////////////////////////////////////////

const double TWeightUnit::PI =  M_PI;


////////////////////////////////////////////////////////////////////////
TWeightUnit::TWeightUnit()
{
  
	ConfigReader = ConfigReaderSubsystem::Instance();
	
	Logger = LoggerSubsystem::Instance();
	
	PDGDatabese = TDatabasePDG::Instance();
  
	ReadConfigFile( string("Generator.dat") );
 
  
}

TWeightUnit::~TWeightUnit(void)
{
	
};

////////////////////////////////////////////////////////////////////////
void TWeightUnit::ReadConfigFile( const string & filename )
{

	
	//Parse variables: setting up model variables	
	//nothing to do
	
	return;
};

////////////////////////////////////////////////////////////////////////
double TWeightUnit::getAcceptance( TEvent * event )
{
    
       return( 1.0 );

};

////////////////////////////////////////////////////////////////////////
bool TWeightUnit::getIsAccepted( TEvent * event )
{
	bool isAccepted = true;
	
    return isAccepted;

};

////////////////////////////////////////////////////////////////////////
complex<double> TWeightUnit::getMatrixElement( TEvent * event )
{
	
	return( 1.0 );
	
};


////////////////////////////////////////////////////////////////////////
void TWeightUnit::SetEvent( TEvent * event, double eventWeight )
{
	
};

////////////////////////////////////////////////////////////////////////
double TWeightUnit::GetWeight( TEvent * event, double eventWeight )
{
	
	
	double weight = GetWeight( eventWeight );
	
	return( weight );
	
	
};

////////////////////////////////////////////////////////////////////////
double TWeightUnit::GetWeight( double eventWeight )
{
	double wt = eventWeight;

	double wtf = 1.0*wt;
	
	
	return wtf;
	
};
