/***********************************************************************
 
				THistogramf0
  
Class for collecting statistics/histograms of reaction a+b->a+b+c+d


***********************************************************************/

#include "THistogramf0.h"


////////////////////////////////////////////////////////////////////////
THistogramf0::THistogramf0() : xsectionHistogramsCreated(false)
{
	
	PDGDatabese = TDatabasePDG::Instance();
	
	ConfigReader = ConfigReaderSubsystem::Instance();
	
	Logger = LoggerSubsystem::Instance();
	
	ReadConfigFile( string("Generator.dat") );
	
	//particles fourmomentum
	pf = new TLorentzVector [nop+1];
	  

	azim  = new double [nop+1];
	theta = new double [nop+1];
	pt    = new double [nop+1];
	prap  = new double [nop+1];
	
	
	Nevents = 0;
	
	AllocateHistograms();
	
	
};

////////////////////////////////////////////////////////////////////////
THistogramf0::~THistogramf0()
{

	delete [] azim;
	delete [] theta;
	delete [] pt;
	delete [] prap;	
	
	delete [] pf;
	
	DeallocateHistograms();
	
	//if differential cross-section histograms were created
	if( xsectionHistogramsCreated == true )
	{
	  for( int i = 0; i < NXhisto; i++){
			delete XHOGENE[i];
			cout << " delete XHOGENE " << i << endl;}
	}
	
	
};

////////////////////////////////////////////////////////////////////////
void THistogramf0::ReadConfigFile( const string & filename )
{
	//Parse extract model constants: setting up model variables	
	PsFilename = ConfigReader->GetStringValue( string( "PsFilename" ) );
	RootFilename = ConfigReader->GetStringValue( string( "RootFilename" ) );
	nop = ConfigReader->GetIntValue( string( "nop" ) );
	tecm =  ConfigReader->GetDoubleValue( "tecm" );
        mres = ConfigReader->GetDoubleValue( "mres" );
        gres = ConfigReader->GetDoubleValue( "gres" );
	
	PsXFilename = ConfigReader->GetStringValue( string( "PsXSectionFilename" ) );
	RootXFilename = ConfigReader->GetStringValue( string( "RootXSectionFilename" ) );
	
	//Log configuration
	string ConfigLog = ConfigReader->GetStringValue( string("ConfigLogFile") );
	TLog * ConfigLogger = new TLog();
	ConfigLogger->setLogTime( false );
	ConfigLogger->openLogFile( ConfigLog );

	ConfigLogger->Write() << "# " << "-----------" << endl;
	ConfigLogger->Write() << "# " << "THistogramf0 " << endl;
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
void THistogramf0::SetOutputFile( const string & filename )
{
	
	PsFilename = string( filename ) + string(".ps");
	RootFilename = string( filename ) + string(".root");
	
};

////////////////////////////////////////////////////////////////////////
void THistogramf0::AllocateHistograms( void )
{ 
 
  HOGENE[0]=new TH1F("pt-3","transverse momentum of 3 ",100,0.0,2.0); 
  HOGENE[1]=new TH1F("pt-1&2","transverse momentum of 1 and 2 ",100,0.0,2.0);
  HOGENE[2]=new TH1F("prap-3"," rapidity 3",100,-6.0,6.0); 
  HOGENE[3]=new TH1F("prap-1&2","rapidity 1 and 2 ",100,-10.0,10.0); 
  HOGENE[4]=new TH1F("azimcor","azimuthal corr. 1-2",100,0.0,180.0);
  HOGENE[5]=new TH1F("coll12","1+2 collinearity;[mrad];# events",20,0.0,6.0);
  HOGENE[6]=new TH1F("t1","momentum transfer t_{1};[GeV^{2}];# events",100,-1.0,0.0);
  HOGENE[7]=new TH1F("t2","momentum transfer t_{2};[GeV^{2}];# events",100,-1.0,0.0);
  HOGENE[8]=new TH1F("s13","log10 s13",100,0.0,5.0);
  HOGENE[9]=new TH1F("s23","log10 s24",100,0.0,5.0);
  HOGENE[10]=new TH1F("xf"," Feynman x_{F};x_{F};# events",100,0.9,1.0);
  Nhisto=11;
  NXhisto=2;
  XM[0]=0;
  XM[1]=2;
  for( int i=0;i<Nhisto;i++ ) HOGENE[i]->Sumw2();
  
   
};

////////////////////////////////////////////////////////////////////////
void THistogramf0::FillHistograms( void )
{

  double wgt = EventWeight;
  
  HOGENE[0]->Fill(pf[3].Pt(),wgt);
  HOGENE[1]->Fill(pf[1].Pt(),wgt);
  HOGENE[1]->Fill(pf[2].Pt(),wgt);
  HOGENE[2]->Fill(prap[3],wgt);
  HOGENE[3]->Fill(prap[1],wgt);
  HOGENE[3]->Fill(prap[2],wgt);
  HOGENE[4]->Fill(dazim12,wgt);
  HOGENE[5]->Fill(coll12,wgt);
       HOGENE[6]->Fill(t1,wgt); 
       HOGENE[7]->Fill(t2,wgt); 
       HOGENE[8]->Fill(log10(s13),wgt);        
       HOGENE[9]->Fill(log10(s23),wgt);
       HOGENE[10]->Fill(dz1,wgt);
       HOGENE[10]->Fill(dz2,wgt);
 };

////////////////////////////////////////////////////////////////////////
void THistogramf0::WriteHistograms( double xsection )
{
		
	//create differential cross-section histograms
	//transform mb to nb
	//CreateXSectionHistograms( 1e6 * xsection );
  //cross section in mikro-barns
	CreateXSectionHistograms( 389.0*xsection );
	
  
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
       c22->cd(1);
       HOGENE[0]->Draw();
       c22->cd(2);
       HOGENE[1]->Draw();
       c22->cd(3);
       HOGENE[2]->Draw();
       c22->cd(4);
       HOGENE[3]->Draw();
     c22->Update();
     ps->NewPage();
       c22->cd(1);
       HOGENE[4]->Draw();
       c22->cd(2);
       HOGENE[5]->Draw();
       c22->cd(3);
       HOGENE[6]->Draw();
       c22->cd(4);
       HOGENE[7]->Draw();
     c22->Update();
     ps->NewPage();
       c22->cd(1);
       HOGENE[8]->Draw();
       c22->cd(2);
       HOGENE[9]->Draw();
       c22->cd(3);
       HOGENE[10]->Draw();
      c22->Update();
      cout << " HOGENE written " << endl;       
     //write x-section histograms
       ps->NewPage();
       c22->cd(1);
       XHOGENE[0]->Draw();
       c22->cd(2);
       XHOGENE[1]->Draw();
       c22->Update();
      cout << " XHOGENE written " << endl;       
    
     
     ps->Close();
    
    delete ps;
    delete c22;
    delete c21;
    delete c31;
     
     //ROOT file for histograms
     TFile * tf = new TFile( RootFilename.c_str(),"RECREATE", "Histograms from generator");
     
     for( int i = 0; i < Nhisto; i++ )  {
		HOGENE[ i ] -> Write();
		cout << " write histogram " << i << endl;}     
     for( int i = 0; i < NXhisto; i++ )  {
		XHOGENE[ i ] -> Write();
		cout << " write X-histogram " << i << endl;}     
      
     tf->Write();
     tf->Close();
     
     delete tf;
     
};

////////////////////////////////////////////////////////////////////////
void THistogramf0::DeallocateHistograms(void )
{
	
    for( int i = 0; i < Nhisto; i++ )     delete HOGENE[i];
    //for( int i = 0; i < NXhisto; i++ )     delete XHOGENE[i];
    return;
};

////////////////////////////////////////////////////////////////////////
void THistogramf0::Fill( TEvent * event, double weight  )
{
	//simple check - more elaborate checks could decrease efficeincy
	assert( event->GetNIn() == 2 );
	assert( event->GetNOut() >= 3 );
	
	EventWeight = weight;

	ExtractParticles( event );
	
	Kinematics2to3();
	
	FillHistograms();
	
	Nevents++;
	

	return;
	
};

////////////////////////////////////////////////////////////////////////
void THistogramf0::ExtractParticles( TEvent * event )
{
	
	assert( nop >= 3 );
	
	pb1 = event->pb[1];
	pb2 = event->pb[2];
	
	
	if( pb1.Z()*event->pf[1].Z() > 0 )
	{
			pf[1]=event->pf[1];
			pf[2]=event->pf[2];
	}
	else
	{
			pf[1]=event->pf[2];
			pf[2]=event->pf[1];
	}
	
	pf[3]=event->pf[3];

};

////////////////////////////////////////////////////////////////////////
void THistogramf0::Kinematics2to3(void )
{

	
	//radian to angle   
    double r2a = 180.0 / M_PI;

    for(int i=1;i<nop+1;i++)
    { 
		prap[i] = pf[i].Rapidity();
		azim[i] = pf[i].Phi();
		theta[i] = pf[i].Theta();
		pt[i] = pf[i].Pt(); 
    }
    
    //tecm = ( pb1 + pb2 ).M();
    
    double inp1 = sqrt( tecm * tecm / 4.0 - mProton * mProton );
    
    //four body kinematics kinematics
    TLorentzVector ptot = pf[1] + pf[2] + pf[3];
    double stot = ptot.M2();
    t1 = ( pb1 - pf[1] ).M2();
    t2 = ( pb2 - pf[2] ).M2();
    s13 = ( pf[1] + pf[3]).M2();
    s23 = ( pf[2] + pf[3]).M2();
    dazim12 = ( azim[1] - azim[2] );
    dazim12 = r2a*fdelta( dazim12 );
    coll12 = 1000.0 * ( M_PI - pf[1].Angle( pf[2].Vect()));//collinearity between protons
    dz1 = pf[1].Pz()/inp1;
    dz2 = pf[2].Pz()/inp1;
};

////////////////////////////////////////////////////////////////////////
double THistogramf0::fdelta(double dazim)
{
	double x = fmod( dazim + M_PI, 2.0*M_PI );
    if (x < 0)
        x += 2.0*M_PI;
    return x - M_PI;
    
};

////////////////////////////////////////////////////////////////////////
void THistogramf0::CreateXSectionHistograms( double xsection)
{
	
	xsectionHistogramsCreated = true;
			
	double binWidth = 0.0;
	double x = 0.0;
	double val = 0.0;
	double error = 0.0;
	int ihisto;
	double tmp = 0.0;  
        for( int i=0; i < NXhisto;i++){
	  ihisto=XM[i];
		TH1F * tmpHist = (TH1F*) HOGENE[ihisto]->Clone();
		TH1F * tmpHist2 = (TH1F*) HOGENE[ihisto]->Clone();
		tmpHist->Reset("ice");
		tmpHist->GetSumw2()->Set(0);
		
		int nbinsx = tmpHist2->GetNbinsX();
		//multiply by luminosity
		tmpHist2->Scale( 1.0/double( Nevents ) );
		tmpHist2->Scale( xsection );
		
		for( int j = 1; j < nbinsx; j++ ){
		  binWidth= tmpHist2->GetBinWidth(j);
		  val = tmpHist2->GetBinContent( j );
		  error = tmpHist2->GetBinError( j );
		  tmpHist->SetBinContent(j,val/binWidth);
		  tmpHist->SetBinError(j,error/binWidth);}
		XHOGENE[i] = tmpHist;
		cout << " XHOGENE " << ihisto << endl;
	}

	return;
	
};
