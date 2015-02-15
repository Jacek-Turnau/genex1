/***********************************************************************
*                          Class TAcceptance                                
*                                                                           
*  Interface class for different acceptances.                                               
*   									       
************************************************************************/
#include"TAcceptance.h"


////////////////////////////////////////////////////////////////////////
TAcceptance:: TAcceptance( TWeightStrategy * wght ) 
{
  
	weight = wght;  
	
	ConfigReader = ConfigReaderSubsystem::Instance();
 
	PDGDatabese = TDatabasePDG::Instance();
	
};

////////////////////////////////////////////////////////////////////////
TAcceptance::~TAcceptance(void)
{

	delete weight;
  
};

////////////////////////////////////////////////////////////////////////
void TAcceptance::ReadConfigFile( const string & filename )
{
	
};

////////////////////////////////////////////////////////////////////////
double TAcceptance::getAcceptance( TEvent * event )
{
      
      return 1.0;

};

////////////////////////////////////////////////////////////////////////
bool TAcceptance::getIsAccepted( TEvent * event )
{
      return true;
};

////////////////////////////////////////////////////////////////////////
complex<double> TAcceptance::getMatrixElement( TEvent * event )
{
	
	return( weight->getMatrixElement( event ) * getAcceptance( event ) );
	
};

////////////////////////////////////////////////////////////////////////
double TAcceptance::GetWeight( TEvent * event, double eventWeight )
{
	
	return( weight->GetWeight( event, eventWeight ) * getAcceptance( event ) );
	
};
