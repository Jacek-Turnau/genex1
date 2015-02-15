/***********************************************************************
 
				Class THistogramN
  
Class for collecting statistics/histograms of reaction a+b->a+b+d1+d2+...+ (N-2)d.

***********************************************************************/

#include "THistogramN.h"



////////////////////////////////////////////////////////////////////////
THistogramN::THistogramN()
{
	
	PDGDatabese = TDatabasePDG::Instance();
	
	ConfigReader = ConfigReaderSubsystem::Instance();
	
	Logger = LoggerSubsystem::Instance();
	
	ReadConfigFile( string("Generator.dat") );
	
	AllocateHistograms();
	
};

////////////////////////////////////////////////////////////////////////
THistogramN::~THistogramN()
{
	
	DeallocateHistograms();	

};

////////////////////////////////////////////////////////////////////////
void THistogramN::ReadConfigFile( const string & filename )
{

	//Parse extract model constants: setting up model variables	
	PsFilename = ConfigReader->GetStringValue( string( "PsFilename" ) );
	RootFilename = ConfigReader->GetStringValue( string( "RootFilename" ) );
	
	PsXFilename = ConfigReader->GetStringValue( string( "PsXSectionFilename" ) );
	RootXFilename = ConfigReader->GetStringValue( string( "RootXSectionFilename" ) );
	
	nop = ConfigReader->GetIntValue( "nop" );
	
	//Log configuration
	string ConfigLog = ConfigReader->GetStringValue( string("ConfigLogFile") );
	TLog * ConfigLogger = new TLog();
	ConfigLogger->setLogTime( false );
	ConfigLogger->openLogFile( ConfigLog );

	ConfigLogger->Write() << "# " << "-----------" << endl;
	ConfigLogger->Write() << "# " << "THistogramN " << endl;
	ConfigLogger->Write() << "# " << "-----------" << endl;
	ConfigLogger->Write() << "PsFilename = " << PsFilename << endl;
	ConfigLogger->Write() << "RootFilename = " << RootFilename << endl;
	ConfigLogger->Write() << "tecm = " << tecm << endl;
	ConfigLogger->Write() << "nop = " << nop << endl;
	ConfigLogger->Write() << "PsXSectionFilename = " << PsXFilename << endl;
	ConfigLogger->Write() << "RootXSectionFilename = " << RootXFilename << endl;
	
	delete ConfigLogger;

	
	return;
};

////////////////////////////////////////////////////////////////////////
void THistogramN::SetOutputFile( const string & filename )
{
	
	PsFilename = string( filename ) + string(".ps");
	RootFilename = string( filename ) + string(".root");
	
};

////////////////////////////////////////////////////////////////////////
void THistogramN::AllocateHistograms( void )
{

  hEtaPhi = new TH2F("EtaPhi3", "Eta vs Phi 1,2,CM", 50,0,0, 50,0,0);
  hEtaPhi->SetBit(TH1::kCanRebin);
  hEtaPhi->SetXTitle("#eta");
  hEtaPhi->SetYTitle("#phi"); 
  
  hEtaPhiCentral = new TH2F("EtaPhi1", "Eta vs Phi CM", 50,0,0, 50,0,0);
  hEtaPhiCentral->SetBit(TH1::kCanRebin);
  hEtaPhiCentral->SetXTitle("#eta");
  hEtaPhiCentral->SetYTitle("#phi"); 
  
  
  hEtaPhiPP = new TH2F("EtaPhi2", "Eta vs Phi 1 + 2", 50,0,0, 50,0,0);
  hEtaPhiPP->SetBit(TH1::kCanRebin);
  hEtaPhiPP->SetXTitle("#eta");
  hEtaPhiPP->SetYTitle("#phi"); 
  
};

////////////////////////////////////////////////////////////////////////
void THistogramN::FillHistograms( void )
{
  
};

////////////////////////////////////////////////////////////////////////
void THistogramN::WriteHistograms( double xsection )
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
     hEtaPhiCentral->Draw();
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
     
     
     hEtaPhiCentral->Write();
     hEtaPhi->Write();
	 hEtaPhiPP->Write();     
     
     
     tf->Write();
     tf->Close();
     
     delete tf;
     
};

////////////////////////////////////////////////////////////////////////
void THistogramN::DeallocateHistograms(void )
{
  
  delete hEtaPhi;
  delete hEtaPhiCentral;
  delete hEtaPhiPP;
    
  return;
  
};

////////////////////////////////////////////////////////////////////////
void THistogramN::Fill( TEvent * event, double weight  )
{
	//simple check - more elaborate checks could decrease efficeincy
	assert( event->GetNIn() == 2 );
	assert( event->GetNOut() > 2 );
	
	TLorentzVector pb1 = event->pb[1];
	TLorentzVector pb2 = event->pb[2];
	
	
	TLorentzVector pf[ nop+1];
	
	for( int i = 1; i < nop+1; i++ )
		pf[i] = event->pf[i];
	
    for( int i = 1; i< nop+1; i++ )
    {
		hEtaPhi->Fill( pf[i].Eta(), pf[i].Phi(), weight);
	}
    
    for( int i = 2; i < nop+1; i++ )
    {
		hEtaPhiCentral->Fill( pf[i].Eta(), pf[i].Phi(), weight);
	}
	
	hEtaPhiPP->Fill( pf[1].Eta(), pf[1].Phi(), weight);
	hEtaPhiPP->Fill( pf[2].Eta(), pf[2].Phi(), weight);
	

	return;
	
};
