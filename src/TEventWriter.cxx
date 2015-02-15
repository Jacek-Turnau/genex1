/***********************************************************************

			TEventWriter
  
Writes events to files for further investigation.


***********************************************************************/

#include "TEventWriter.h"


////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////
TEventWriter::TEventWriter( TEvent * event )
{
	
	ConfigReader = new TPolicyReader();
	
	Logger = LoggerSubsystem::Instance();
  
	ReadConfigFile( string("Generator.dat") );
	
	eventPointer = event;
	
	eventObject = new TEvent();
	
	if( SaveEventDataLHE_XML )
	{
	
		writer = new LHEF::Writer( lheFilename );
		const int protonPDGid    = 2212;
		writer->heprup.IDBMUP = make_pair( protonPDGid, protonPDGid );
		writer->heprup.EBMUP = make_pair( tecm/2.0 , tecm /2.0);
		writer->heprup.IDWTUP = 3; //ewents with unit weight
		writer->heprup.generators.push_back( make_pair( string("GenEx"), string("1.0") ) );

		writer->initComments () << "GenEx particle generator" << endl;

		writer->init();
	}	
		
	if( SaveEventDataRoot )
	{
	
		tree = new TTree("events","Event tree");
	
		string basename("");
		string number(""); 
	
		//in particles
		basename = "In particle ";
	
		for( int i = 1; i < event->GetNIn() + 1; i++ )
		{
			number = "";		
			number = static_cast<ostringstream*>( &(ostringstream() << i) )->str();
		
			tree->Branch( ( basename + number ).c_str() ,"TLorentzVector", &(eventPointer->pb[i]) );
			//tree->Branch( ( basename + number ).c_str() ,"TParticlePDG ", &(eventPointer->pbInfo[i]) );
			//tree->Branch( ( basename + number ).c_str() ,"TLorentzVector", &(eventObject->pb[i]) );
			//tree->Branch( ( basename + number ).c_str() ,"TParticlePDG ", &(eventObject->pbInfo[i]) );
		
		}
	
		//out particles
		basename = "Out particle ";
	
		for( int i = 1; i < event->GetNOut() + 1; i++ )
		{
			number = "";		
			number = static_cast<ostringstream*>( &(ostringstream() << i) )->str();
		
			tree->Branch( ( basename + number ).c_str() ,"TLorentzVector", &(eventPointer->pf[i]) );
			//tree->Branch( ( basename + number ).c_str() ,"TParticlePDG ", &(eventPointer->pfInfo[i]) );
			//tree->Branch( ( basename + number ).c_str() ,"TLorentzVector", &(eventObject->pf[i]) );
			//tree->Branch( ( basename + number ).c_str() ,"TParticlePDG ", &(eventObject->pfInfo[i]) );
		}
	
	}
	
	if( SaveEventDataLHE_TXT )
	{
		//open file to write date
		datafile.open( lheTxtBackupFile.c_str() );
		eventObject = event;
	}
	
	num = 0;

	
	return;
};
////////////////////////////////////////////////////////////////////////
void TEventWriter::ReadConfigFile( const string & filename )
{
	
	//Reading conf file
	ConfigReader->ReadConfFile(filename);
	
	#ifdef DEBUG
		cout << "-----------" << endl;
		cout << "FILE: " << filename << endl;
		cout << "-----------" << endl;
		cout << *ConfigReader << endl;
		cout << "-----------" << endl;
	#endif	
	
	//Parse variables: setting up model variables	
	
	//option of integration
	
	SaveEventDataRoot = ConfigReader->GetBoolValue( string( "SaveEventDataRoot" ) );
	rootFilename = ConfigReader->GetStringValue( string( "RootEventFile" ) );
	SaveEventDataLHE_XML = ConfigReader->GetBoolValue( string( "SaveEventDataLHE_XML" ) );
	lheFilename  = ConfigReader->GetStringValue( string( "LHEXMLEventFile" ) );
	SaveEventDataLHE_TXT = ConfigReader->GetBoolValue( string( "SaveEventDataLHE_TXT" ) );
	lheTxtBackupFile = ConfigReader->GetStringValue( string( "LHETXTEventFile" ) );
	tecm = ConfigReader->GetDoubleValue( string( "tecm" ) );
	weightStrategy = ConfigReader->GetIntValue( "weightStrategy" );
	boost = ConfigReader->getBoost();


	//Log configuration
	string ConfigLog = ConfigReader->GetStringValue( string("ConfigLogFile") );
	TLog * ConfigLogger = new TLog();
	ConfigLogger->setLogTime( false );
	ConfigLogger->openLogFile( ConfigLog );
 
	ConfigLogger->Write() << "# " << "-----------" << endl;
	ConfigLogger->Write() << "# " << "TEventWriter " << endl;
	ConfigLogger->Write() << "# " << "-----------" << endl;
	ConfigLogger->Write() << "SaveEventDataRoot = " << SaveEventDataRoot << endl;
	ConfigLogger->Write() << "RootEventFile = " << rootFilename << endl;
	ConfigLogger->Write() << "SaveEventDataLHE_XML = " << SaveEventDataLHE_XML << endl;
	ConfigLogger->Write() << "LHEXMLEventFile = " << lheFilename << endl;
	ConfigLogger->Write() << "SaveEventDataLHE_TXT = " << SaveEventDataLHE_TXT << endl;
	ConfigLogger->Write() << "LHETXTEventFile = " << lheTxtBackupFile << endl;
	ConfigLogger->Write() << "tecm = " << tecm << endl;
	ConfigLogger->Write() << "weightStrategy = " << weightStrategy << endl;
	
	delete ConfigLogger;



	#ifdef DEBUG
		cout << "-----------" << endl;
		cout << "INITIALIZED VARIABLES" << endl;
		cout << "-----------" << endl;
		
		cout << "EventsBackupFile = " << rootFilename << endl;
		cout << "LHEBackupFile = " << lheFilename << endl;
		cout << "tecm = " << tecm << endl;
		cout << "weightStrategy  = " << weightStrategy << endl;
		cout << "-----------" << endl;
	#endif	
	
	
	return;
};


////////////////////////////////////////////////////////////////////////
TEventWriter::~TEventWriter()
{
	if( SaveEventDataRoot )
	{
		delete tree;
	}
	
	if( SaveEventDataLHE_TXT )
	{
		datafile.close();
	}
	
	if( SaveEventDataLHE_XML )
	{
		delete writer;
	}
	
	
};

////////////////////////////////////////////////////////////////////////
void TEventWriter::Fill( void )
{
	
	//update particle info
	for( int i = 1; i < eventObject->GetNIn() + 1; i++ )
	{
		eventObject->pb[i] = eventPointer->pb[i];
		eventObject->pb[i].Boost( boost );
		eventObject->pbInfo[i] = eventPointer->pbInfo[i];	
		
		eventPointer->pb[i].Boost(boost);
		
	}
	
	for( int i = 1; i < eventObject->GetNOut() + 1; i++ )
	{
		eventObject->pf[i] = eventPointer->pf[i];
		eventObject->pf[i].Boost( boost );
		eventObject->pfInfo[i] = eventPointer->pfInfo[i];	
		
		eventPointer->pf[i].Boost(boost);
		
	}
	
	//fill root tree
	if( SaveEventDataRoot )
	{
		tree->Fill();
	}
	
	//undo boost on event
	for( int i = 1; i < eventObject->GetNIn() + 1; i++ )
	{
		eventPointer->pb[i].Boost(-boost);	
	}
	
	for( int i = 1; i < eventObject->GetNOut() + 1; i++ )
	{
		eventPointer->pf[i].Boost(-boost);	
	}
	
	///////
	
	int stat = 0;	 
	int m1   = 0;  // mother beginning entries
	int m2   = 0;  // mother ending entries
	int d1   = 0;  // daughter beggining entries
	int d2   = 0;  // daughter ending entries
	
	//space-time vertex position
	double vx = 0.0;
	double vy = 0.0;
	double vz = 0.0;
	double vt = 0.0;
	
	int nup;
	
	int i, nevhep;
	
	///////fill text file
	if( SaveEventDataLHE_TXT )
	{
	
		//space-time vertex position - it is fixed at the center
		vx = 0.0;
		vy = 0.0;
		vz = 0.0;
		vt = 0.0;

		// HEPMC status code; "0 - null entry, 1 - final state (not decayed or fragmented), 2 -	decayed or fragmented; Event generators do not all use status code 2 in the same way. The status code 2 should be used only for physical particles that have been decayed. It is expected that status code 1 will be changed to status 2 if a particle that is not decayed by an event generator is subsequently decayed, e.g., during detector simulation." From: http://lcgapp.cern.ch/project/simu/HepMC/205/status.html
		stat = 0;	 
		m1   = 0;  // mother beginning entries
		m2   = 0;  // mother ending entries
		d1   = 3;  // daughter beggining entries
		d2   = 6;  // daughter ending entries
	
		num += 1;
		nup = eventObject->GetNIn() + eventObject->GetNOut();

		//start record
	
		datafile << num << "\t" << nup << endl;
	
		//In particles
	
		stat = -1;
	
		i = 1;
		nevhep = 1;
	
		datafile << nevhep << "\t" << eventObject->pbInfo[i]->PdgCode() << "\t" << stat << "\t" << m1 << "\t" << m2 << "\t" << d1 << "\t" << d2 << "\t" << eventObject->pb[i].Px() << "\t" << eventObject->pb[i].Py() << "\t" << eventObject->pb[i].Pz() << "\t" << eventObject->pb[i].E() << "\t" << eventObject->pb[i].M() << "\t" <<  vx << "\t" <<  vy << "\t" <<  vz << "\t" <<  vt << endl;
	
		i = 2;
		nevhep = 2;
	
		datafile << nevhep << "\t" << eventObject->pbInfo[i]->PdgCode() << "\t" << stat << "\t" << m1 << "\t" << m2 << "\t" << d1 << "\t" << d2 << "\t" << eventObject->pb[i].Px() << "\t" << eventObject->pb[i].Py() << "\t" << eventObject->pb[i].Pz() << "\t" << eventObject->pb[i].E() << "\t" << eventObject->pb[i].M() << "\t" <<  vx << "\t" <<  vy << "\t" <<  vz << "\t" <<  vt << endl;
	
		//Out particles
	
		stat = 1;
	
		//adding final particles:
		m1 = 1;
		m2 = 1;
		d1 = 0;
		d2 = 0;
	
		i=1;
		nevhep = 3;
	
		datafile << nevhep << "\t" << eventObject->pfInfo[i]->PdgCode() << "\t" << stat << "\t" << m1 << "\t" << m2 << "\t" << d1 << "\t" << d2 << "\t" << eventObject->pf[i].Px() << "\t" << eventObject->pf[i].Py() << "\t" << eventObject->pf[i].Pz() << "\t" << eventObject->pf[i].E() << "\t" << eventObject->pf[i].M() << "\t" <<  vx << "\t" <<  vy << "\t" <<  vz << "\t" <<  vt << endl;
	
		m1 = 2;
		m2 = 2;	
	
		i = 2;
		nevhep = 4;
	
		datafile << nevhep << "\t" << eventObject->pfInfo[i]->PdgCode() << "\t" << stat << "\t" << m1 << "\t" << m2 << "\t" << d1 << "\t" << d2 << "\t" << eventObject->pf[i].Px() << "\t" << eventObject->pf[i].Py() << "\t" << eventObject->pf[i].Pz() << "\t" << eventObject->pf[i].E() << "\t" << eventObject->pf[i].M() << "\t" <<  vx << "\t" <<  vy << "\t" <<  vz << "\t" <<  vt << endl;
	
		m1 = 1;
		m2 = 2;
	
		//i = 3;
		nevhep = 5;
	
		for( int i = 3; i < eventObject->GetNOut() + 1; i++ )
		{
	
			datafile << nevhep << "\t" << eventObject->pfInfo[i]->PdgCode() << "\t" << stat << "\t" << m1 << "\t" << m2 << "\t" << d1 << "\t" << d2 << "\t" << eventObject->pf[i].Px() << "\t" << eventObject->pf[i].Py() << "\t" << eventObject->pf[i].Pz() << "\t" << eventObject->pf[i].E() << "\t" << eventObject->pf[i].M() << "\t" <<  vx << "\t" <<  vy << "\t" <<  vz << "\t" <<  vt << endl;
		
			nevhep++;
	
		}
	}
	
	
	//// Fill LHA structure
	if( SaveEventDataLHE_XML )
	{
	
		writer->hepeup.NUP = eventObject->GetNOut() + eventObject->GetNIn();
    
		// HEPMC status code; "0 - null entry, 1 - final state (not decayed or fragmented), 2 -	decayed or fragmented; Event generators do not all use status code 2 in the same way. The status code 2 should be used only for physical particles that have been decayed. It is expected that status code 1 will be changed to status 2 if a particle that is not decayed by an event generator is subsequently decayed, e.g., during detector simulation." From: http://lcgapp.cern.ch/project/simu/HepMC/205/status.html
		stat = 0;	 
		m1   = 0;  // mother beginning entries
		m2   = 0;  // mother ending entries
		d1   = 3;  // daughter beggining entries
		d2   = 6;  // daughter ending entries
	
		vector<double> p;
	
		//In particles	 
		stat = -1;
	
		//beam proton1
		i=1;
		writer->hepeup.IDUP.push_back( eventObject->pbInfo[i]->PdgCode() );
		writer->hepeup.ISTUP.push_back( stat );
		writer->hepeup.MOTHUP.push_back( make_pair( m1, m2) );
		p.clear();
		p.push_back( eventObject->pb[i].Px() );
		p.push_back( eventObject->pb[i].Py() );
		p.push_back( eventObject->pb[i].Pz() );
		p.push_back( eventObject->pb[i].E()  );
		p.push_back( eventObject->pb[i].M()  );
		writer->hepeup.PUP.push_back( p );
		writer->hepeup.VTIMUP.push_back( 9 );
		writer->hepeup.SPINUP.push_back( 9 );
	
		//beam proton2
		i=2;
		writer->hepeup.IDUP.push_back( eventObject->pbInfo[i]->PdgCode() );
		writer->hepeup.ISTUP.push_back( stat );
		writer->hepeup.MOTHUP.push_back( make_pair( m1, m2) );
		p.clear();
		p.push_back( eventObject->pb[i].Px() );
		p.push_back( eventObject->pb[i].Py() );
		p.push_back( eventObject->pb[i].Pz() );
		p.push_back( eventObject->pb[i].E()  );
		p.push_back( eventObject->pb[i].M()  );
		writer->hepeup.PUP.push_back( p );
		writer->hepeup.VTIMUP.push_back( 9 );
		writer->hepeup.SPINUP.push_back( 9 );
	
		//Out particles
	
		stat = 1;
	
		//adding final particles:
		m1 = 1;
		m2 = 2;
		d1 = 0;
		d2 = 0;
	
		for( int i = 1; i < eventObject->GetNOut() +1; i++ )
		{

			writer->hepeup.IDUP.push_back( eventObject->pfInfo[i]->PdgCode() );
			writer->hepeup.ISTUP.push_back( stat );
			writer->hepeup.MOTHUP.push_back( make_pair( m1, m2) );
			writer->hepeup.ISTUP.push_back( stat );
			p.clear();
			p.push_back( eventObject->pf[i].Px() );
			p.push_back( eventObject->pf[i].Py() );
			p.push_back( eventObject->pf[i].Pz() );
			p.push_back( eventObject->pf[i].E()  );
			p.push_back( eventObject->pf[i].M()  );
			writer->hepeup.PUP.push_back( p );
			writer->hepeup.VTIMUP.push_back( 9 );
			writer->hepeup.SPINUP.push_back( 9 );

		}	
	
		writer->writeEvent();	 
	 
		writer->hepeup.clear(); 
		writer->hepeup.IDUP.clear(); 
		writer->hepeup.ISTUP.clear(); 
		writer->hepeup.MOTHUP.clear(); 
		writer->hepeup.ISTUP.clear(); 
		writer->hepeup.PUP.clear(); 
		writer->hepeup.VTIMUP.clear(); 
		writer->hepeup.SPINUP.clear(); 
		
		
	}
	
};

////////////////////////////////////////////////////////////////////////
void TEventWriter::Save( bool overwrite )
{
	if( SaveEventDataRoot )
	{
    
		TFile file( rootFilename.c_str(), "Update" ); 
  
		file.SetCompressionLevel(1); 
  
		file.cd();
	
		if( overwrite )
		{
			tree->Write("",TObject::kOverwrite);
		}
		else
		{	
			tree->Write();
		}

		file.Close();
 
	}
};
