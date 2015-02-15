/***********************************************************************
 
				THistogram3
  
Class for Collecting statistics/histograms of reaction a+b->a+b+c.


***********************************************************************/

#include "THistogram3.h"


////////////////////////////////////////////////////////////////////////
THistogram3::THistogram3()
{
	
	PDGDatabese = TDatabasePDG::Instance();
	
	ConfigReader = ConfigReaderSubsystem::Instance();
	
	Logger = LoggerSubsystem::Instance();
	
	ReadConfigFile( string("Generator.dat") );
	
	AllocateHistograms();
	
	
};

////////////////////////////////////////////////////////////////////////
THistogram3::~THistogram3()
{
	
	DeallocateHisograms();
	
};

////////////////////////////////////////////////////////////////////////
void THistogram3::AllocateHistograms( void )
{
	 
  hEtaPhi = new TH2F("EtaPhi3", "Eta vs Phi 1,2,CM", 50,0,0, 50,0,0);
  hEtaPhi->SetBit(TH1::kCanRebin);
  hEtaPhi->SetXTitle("#eta");
  hEtaPhi->SetYTitle("#phi"); 
  
  hEtaPhiPhoton = new TH2F("EtaPhi1", "Eta vs Phi CM", 50,0,0, 50,0,0);
  hEtaPhiPhoton->SetBit(TH1::kCanRebin);
  hEtaPhiPhoton->SetXTitle("#eta");
  hEtaPhiPhoton->SetYTitle("#phi"); 
  
  
  hEtaPhiPP = new TH2F("EtaPhi2", "Eta vs Phi 1 + 2", 50,0,0, 50,0,0);
  hEtaPhiPP->SetBit(TH1::kCanRebin);
  hEtaPhiPP->SetXTitle("#eta");
  hEtaPhiPP->SetYTitle("#phi"); 
	
};

////////////////////////////////////////////////////////////////////////
void THistogram3::DeallocateHisograms( void )
{
	   
  delete hEtaPhi;
  delete hEtaPhiPhoton;
  delete hEtaPhiPP;
	
};

////////////////////////////////////////////////////////////////////////
void THistogram3::ReadConfigFile( const string & filename )
{
	
	//Parse extract model constants: setting up model variables	
	PsFilename = ConfigReader->GetStringValue( string( "PsFilename" ) );
	RootFilename = ConfigReader->GetStringValue( string( "RootFilename" ) );
	
	//Log configuration
	string ConfigLog = ConfigReader->GetStringValue( string("ConfigLogFile") );
	TLog * ConfigLogger = new TLog();
	ConfigLogger->setLogTime( false );
	ConfigLogger->openLogFile( ConfigLog );

	ConfigLogger->Write() << "# " << "-----------" << endl;
	ConfigLogger->Write() << "# " << "THistogram3 " << endl;
	ConfigLogger->Write() << "# " << "-----------" << endl;
	ConfigLogger->Write() << "PsFilename = " << PsFilename << endl;
	ConfigLogger->Write() << "RootFilename = " << RootFilename << endl;
	
	delete ConfigLogger;

	
	return;
};

////////////////////////////////////////////////////////////////////////
void THistogram3::SetOutputFile( const string & filename )
{
	
	PsFilename = string( filename ) + string(".ps");
	RootFilename = string( filename ) + string(".root");
	
};


////////////////////////////////////////////////////////////////////////
void THistogram3::WriteHistograms( double xsection )
{
 
  //Post Scritp file for histograms
  TPostScript *ps = new TPostScript(PsFilename.c_str(),112);
  // Set style
  gROOT->SetStyle("Plain");
  gROOT->ForceStyle();
  gStyle->SetOptStat(1);
  gStyle->SetTitleBorderSize(0);
  gStyle->SetTitleSize(0.04);
  gStyle->SetTitleFont(42, "hxy");      // for histogram and axis titles
  gStyle->SetLabelFont(42, "xyz");      // for axis labels (values)
  gROOT->ForceStyle();

  //set canvas
      TCanvas* c22 = new TCanvas("c22","generator tests");
      TCanvas* c31 = new TCanvas("c31","generator tests");
      TCanvas* c21 = new TCanvas("c21","generator tests");
       c22->Divide(2,2);
       c31->Divide(3,1);
       c21->Divide(2,1);
  //set canvas style
  c22->SetFillColor(0);
  c22->UseCurrentStyle();
  c22->SetBorderMode(0);       // still leaves red frame bottom and right
  c22->SetFrameBorderMode(0);   // need this to turn off red hist frame!
  c22->UseCurrentStyle();
  c21->UseCurrentStyle();
  c31->UseCurrentStyle();
     
     ps->NewPage();
     c22->cd(0);
     hEtaPhi->Draw();
     c22->Update();
     
     ps->NewPage();
     c22->cd(0);
     hEtaPhiPhoton->Draw();
     c22->Update();
     
     ps->NewPage();
     c22->cd(0);
     hEtaPhiPP->Draw();
     c22->Update();
     
     ps->Close();
    
    delete ps;
    delete c22;
    delete c21;
    delete c31;
     
     //ROOT file for histograms
     TFile * tf = new TFile( RootFilename.c_str(),"RECREATE", "Histograms form generator");
     
     
     hEtaPhiPhoton->Write();
     hEtaPhi->Write();
	 hEtaPhiPP->Write();     
     
     
     tf->Write();
     tf->Close();
     
     delete tf;
     
     
};

////////////////////////////////////////////////////////////////////////
void THistogram3::Fill( TEvent * event, double weight  )
{
	//simple check - more elaborate checks could decrease efficeincy
	assert( event->GetNIn() == 2 );
	assert( event->GetNOut() == 3 );
	
	TLorentzVector pb1 = event->pb[1];
	TLorentzVector pb2 = event->pb[2];
	
	
	TLorentzVector Ps = event->pf[1];
	TLorentzVector P  = event->pf[2];
	TLorentzVector Pk = event->pf[3];
	
	
    
    hEtaPhi->Fill( Pk.Eta(), Pk.Phi(), weight);
    hEtaPhi->Fill( Ps.Eta(), Ps.Phi(), weight);
    hEtaPhi->Fill( P.Eta(), P.Phi(), weight);
	
	hEtaPhiPhoton->Fill( Pk.Eta(), Pk.Phi(), weight);
	
	hEtaPhiPP->Fill( Ps.Eta(), Ps.Phi(), weight);
	hEtaPhiPP->Fill( P.Eta(), P.Phi(), weight);

	return;
	
};
