/***********************************************************************

				TEvent
  
Stores 4-momenta of generated particles.

***********************************************************************/

#include "TEvent.h"


////////////////////////////////////////////////////////////////////////
TEvent::TEvent(): NIn(2), NOut(4)
{
	
	PDGDatabese = TDatabasePDG::Instance();
	
	Logger = LoggerSubsystem::Instance();
	
	ConfigReader = ConfigReaderSubsystem::Instance();
  
	ReadConfigFile( string("Generator.dat") );	
	
	
	pb = new TLorentzVector [NIn+1];
	pf = new TLorentzVector [NOut+1];
	
	pbInfo = new TParticlePDG * [NIn+1];
	pfInfo = new TParticlePDG * [NOut+1];
	
	
	return;
};

////////////////////////////////////////////////////////////////////////
TEvent::TEvent(unsigned int NIn, unsigned int NOut): NIn(NIn), NOut(NOut)
{
	
	PDGDatabese = TDatabasePDG::Instance();
	
	ConfigReader = ConfigReaderSubsystem::Instance();	
	
	pb = new TLorentzVector [NIn+1];
	pf = new TLorentzVector [NOut+1];
	
	pbInfo = new TParticlePDG * [NIn+1];
	pfInfo = new TParticlePDG * [NOut+1];
	
	
	return;
};

////////////////////////////////////////////////////////////////////////
void TEvent::ReadConfigFile( const string & filename )
{
	//Parse variables: setting up model variables	

	NOut = ConfigReader->GetIntValue( "nop" );
	NIn  = ConfigReader->GetIntValue( "nip" );
	
	//Log configuration
	string ConfigLog = ConfigReader->GetStringValue( string("ConfigLogFile") );
	TLog * ConfigLogger = new TLog();
	ConfigLogger->setLogTime( false );
	ConfigLogger->openLogFile( ConfigLog );
 
	ConfigLogger->Write() << "# " << "-----------" << endl;
	ConfigLogger->Write() << "# " << "TEvent " << endl;
	ConfigLogger->Write() << "# " << "-----------" << endl;
	ConfigLogger->Write() << "NOut = " << NOut << endl;
	ConfigLogger->Write() << "NIn  = " << NIn << endl;
	
	delete ConfigLogger;
		
	return;
};


////////////////////////////////////////////////////////////////////////
TEvent::~TEvent()
{
	delete [] pb;
	delete [] pf;
	
	delete [] pbInfo;
	delete [] pfInfo;
	
	return;
};

////////////////////////////////////////////////////////////////////////
int TEvent::GetNIn()
{
	assert( NIn >= 0 );
	return( NIn );
};
    
////////////////////////////////////////////////////////////////////////
int TEvent::GetNOut()
{
	assert( NOut >= 0 );
	return( NOut );
};

////////////////////////////////////////////////////////////////////////
void TEvent::PrintParticles()
{
	
	cout << "----------------" << endl;
	cout << "IN PARTICLES:" << endl;
	cout << "----------------" << endl;
	if( NIn  > 0 )
	{
		for (int i = 1;  i < NIn+1; i ++ )
		{
			cout << "Particle " << i << " = " << " 4-momentum = [" << setw(10) << pb[i].E() << ", " << setw(10) << pb[i].Px() << ", " << setw(10) << pb[i].Py() << ", " << setw(10) << pb[i].Pz() << "];"<< endl;
		};
	};
	
	cout << "----------------" << endl;
	cout << "OUT PARTICLES:" << endl;
	cout << "----------------" << endl;
	if( NOut  > 0 )
	{
	
		for (int i = 1;  i < NOut+1; i ++ )
		{
			cout << "Particle " << i << " = " << " 4-momentum = [" << setw(10) << pf[i].E() << ", " << setw(10) << pf[i].Px() << ", " << setw(10) << pf[i].Py() << ", " << setw(10) << pf[i].Pz() << "];"<< endl;
		};
	};
	
	cout << "----------------" << endl;
	
};

////////////////////////////////////////////////////////////////////////
void TEvent::SaveEvent( string filename )
{
  
  
  //open file with append
  ofstream myfile;
  myfile.open( filename.c_str(), fstream::app );
  

  myfile << "----------------" << endl;
  myfile << "IN PARTICLES:" << endl;
  myfile << "----------------" << endl;
	if( NIn  > 0 )
	{
		for (int i = 1;  i < NIn+1; i ++ )
		{
			myfile << "Particle " << i << " = " << " 4-momentum = [" << setw(10) << pb[i].E() << ", " << setw(10) << pb[i].Px() << ", " << setw(10) << pb[i].Py() << ", " << setw(10) << pb[i].Pz() << "];"<< endl;
		};
	};
	
  myfile << "----------------" << endl;
  myfile << "OUT PARTICLES:" << endl;
  myfile << "----------------" << endl;
	if( NOut  > 0 )
	{
	
		for (int i = 1;  i < NOut+1; i ++ )
		{
			myfile << "Particle " << i << " = " << " 4-momentum = [" << setw(10) << pf[i].E() << ", " << setw(10) << pf[i].Px() << ", " << setw(10) << pf[i].Py() << ", " << setw(10) << pf[i].Pz() << "];"<< endl;
		};
	};
	
  myfile << "----------------" << endl;
  myfile << "END OF PARTICLE LIST" << endl;
  myfile << "----------------" << endl;
  
  
  myfile.close();
  
  return;
};
