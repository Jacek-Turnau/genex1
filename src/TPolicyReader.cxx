/***********************************************************************

					TPolicyReader
  

		Class that reads configuration file with prescribed policy.
	Simple class that reads config file into map (key, value).

 File format:
 
 key1 = value1
 key2 = value2
 ...
 
 Lines with comments starts with # or ;

***********************************************************************/

#include "TPolicyReader.h"

		

////////////////////////////////////////////////////////////////////////
void TPolicyReader::ReadConfFile( const std::string& filename )
{

	TConfigReader::ReadConfFile( filename );
	
	boost.SetXYZ( 0.0, 0.0, 0.0 );
	
	
	ApplyPolicy();
	
	ClalculateKinematics();
	
		//read subfile if required
	if( ! GetBoolValue(std::string( "readConfigLog" ) ) )
	{
		//recover name of the model file
		std::string base;
		std::string result;
		std::stringstream buffer;
	
		base = std::string("Model");
		buffer << base;
		buffer << GetIntValue(std::string( "weightStrategy" ) );
		buffer << std::string( "ConfFile" );
		buffer >> result;
	
		//read model file configuration
		TConfigReader::ReadConfFile( GetStringValue( result ) );
	}
	
  	
	
};


////////////////////////////////////////////////////////////////////////
void TPolicyReader::ReadConfFileNoPolicy( const std::string& filename )
{
	TConfigReader::ReadConfFile( filename );
};

////////////////////////////////////////////////////////////////////////
void TPolicyReader::ApplyPolicy( void )
{
	
	//helper veriables for conversion number -> string
	std::string result;
	std::stringstream buffer;
	
	///////////////////kDim - dimension of monte carlo vector
	
	////Set dimension of random number space
	int eventGenerationStrategy  = TConfigReader::GetIntValue( std::string( "eventGenerationStrategy" ) );

	//Apply constraints on some variables

	int kDim;

	switch( eventGenerationStrategy )
	{
		case 4: kDim = 3*TConfigReader::GetIntValue( std::string( "TEventMaker2toN::nop" ) ) - 4; 
				break; 
				
		case 3: kDim = 5;
				break;
		
		case 2: kDim = 8;
				break;
		
		case 1: kDim = 3*TConfigReader::GetIntValue( std::string( "TEventMakerCylindricalPhaseSpace::nop" ) );
				break;
		
		default: std::cout << "Invalid value 'eventGenerationStrategy' - STOP";
				assert(1==0);
				break;
		
	};
	
	result.clear();
	buffer.clear();
	
	buffer << kDim;
	buffer >> result;

	TConfigReader::insert( std::pair<std::string,std::string>( std::string("kDim"), result ) );

	
	////////////////////////nop nip - incoming/outgoing particles
	
	int weightStrategy  = TConfigReader::GetIntValue( std::string( "weightStrategy" ) );
	
	int nip = 0;
	int nop = 0;
	
	switch( weightStrategy )
	{

		case 1: 
		case 2:
		case 3:
		case 4:
				nip = 2;
				nop = 4;
				break;
		case 5:
				nip = 2;
				nop = 4;
				break;
		case 6:
				nip = 2;
				nop = 3;
				break;
				
		default: std::cout << "Invalid value 'weightStrategy' - STOP";
				assert(1==0);
				break;
		
	};
	
	result.clear();
	buffer.clear();
	
	buffer << nip;
	buffer >> result;

	TConfigReader::insert( std::pair<std::string,std::string>( std::string("nip"), result ) );
	
	result.clear();
	buffer.clear();
	
	buffer << nop;
	buffer >> result;

	TConfigReader::insert( std::pair<std::string,std::string>( std::string("nop"), result ) );
	
	
	
	
};

////////////////////////////////////////////////////////////////////////
TVector3 TPolicyReader::getBoost( void )
{
	return( boost );
};

////////////////////////////////////////////////////////////////////////
void  TPolicyReader::ClalculateKinematics( void )
{
	
	double EA, EB;
	
	//Reading frame setup
	int frameType = TConfigReader::GetIntValue( std::string("frameType") );
	
	//for LAB frame calculate tecm and boost vector
	if( frameType == 2 )
	{
		EA = TConfigReader::GetDoubleValue( std::string("EA") );
		EB = TConfigReader::GetDoubleValue( std::string("EB") );
	
	
		int idA = TConfigReader::GetIntValue( std::string("idA") );
		int idB = TConfigReader::GetIntValue( std::string("idB") );
	
		TDatabasePDG * PDGDatabese = TDatabasePDG::Instance();
	
		///setup frame
		double mA = PDGDatabese->GetParticle(idA)->Mass();
		double mB = PDGDatabese->GetParticle(idB)->Mass();
	
		double pzA = sqrt( EA*EA - mA*mA );
		double pzB = -sqrt( EB*EB - mB*mB );
		
		TLorentzVector pA, pB;
		
		pA.SetPxPyPzE( 0.0, 0.0, pzA, EA );
		pB.SetPxPyPzE( 0.0, 0.0, pzB, EB );
		
		boost = (pA+pB).BoostVector();
		double tecm = (pA+pB).M();
	
		//helper veriables for conversion number -> string
		std::string result;
		std::stringstream buffer;
		
		buffer << tecm;
		buffer >> result;
	
		TConfigReader::insert( std::pair<std::string,std::string>( std::string("tecm"), result ) );
	
	}
	
	
};
