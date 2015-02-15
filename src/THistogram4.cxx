/***********************************************************************
 
				THistogram4
  
Class for collecting statistics/histograms of reaction a+b->a+b+c+d


***********************************************************************/

#include "THistogram4.h"


////////////////////////////////////////////////////////////////////////
THistogram4::THistogram4() : xsectionHistogramsCreated(false)
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
THistogram4::~THistogram4()
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
		for( int i = 0; i < 17; i++)
			delete XHOGENE[i];
	}
	
	
};

////////////////////////////////////////////////////////////////////////
void THistogram4::ReadConfigFile( const string & filename )
{
	//Parse extract model constants: setting up model variables	
	PsFilename = ConfigReader->GetStringValue( string( "PsFilename" ) );
	RootFilename = ConfigReader->GetStringValue( string( "RootFilename" ) );
	nop = ConfigReader->GetIntValue( string( "nop" ) );
	tecm =  ConfigReader->GetDoubleValue( "tecm" );
	
	PsXFilename = ConfigReader->GetStringValue( string( "PsXSectionFilename" ) );
	RootXFilename = ConfigReader->GetStringValue( string( "RootXSectionFilename" ) );
	
	//Log configuration
	string ConfigLog = ConfigReader->GetStringValue( string("ConfigLogFile") );
	TLog * ConfigLogger = new TLog();
	ConfigLogger->setLogTime( false );
	ConfigLogger->openLogFile( ConfigLog );

	ConfigLogger->Write() << "# " << "-----------" << endl;
	ConfigLogger->Write() << "# " << "THistogram4 " << endl;
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
void THistogram4::SetOutputFile( const string & filename )
{
	
	PsFilename = string( filename ) + string(".ps");
	RootFilename = string( filename ) + string(".root");
	
};

////////////////////////////////////////////////////////////////////////
void THistogram4::AllocateHistograms( void )
{ 
 
  HOGENE[0]=new TH1F("pt-3+4","transverse momentum 3+4 ",100,0.0,2.0); 
  HOGENE[1]=new TH1F("pt-1+2","transverse momentum of 1+2 ",100,0.0,2.0);
  HOGENE[2]=new TH1F("prap-3+4"," rapidity 3+4 ",100,-10.0,10.0); 
  HOGENE[3]=new TH1F("prap-1+2","rapidity 1+2 ",100,-10.0,10.0); 
  HOGENE[4]=new TH1F("azimcor","azimuthal corr. 1+2",100,0.0,180.0);
  HOGENE[5]=new TH1F("coll12","1+2 collinearity;[mrad];# events",20,0.0,6.0);
  HOGENE[6]=new TH1F("t1","momentum transfer t_{1};[GeV^{2}];# events",100,-1.0,0.0);
  HOGENE[7]=new TH1F("t2","momentum transfer t_{2};[GeV^{2}];# events",100,-1.0,0.0);
  HOGENE[8]=new TH1F("ta","momentum transfer t_{a};[GeV^{2}];# events",100,-2.0,0.0);
  HOGENE[9]=new TH1F("tb","momentum transfer t_{b};[GeV^{2}];# events",100,-2.0,0.0);
  HOGENE[10]=new TH1F("s13","log10 s13",100,0.0,3.0);
  HOGENE[11]=new TH1F("s24","log10 s24",100,0.0,3.0);
  HOGENE[12]=new TH1F("s14","log10 s14",100,0.0,3.0);
  HOGENE[13]=new TH1F("s23","log10 s23",100,0.0,3.0);
  HOGENE[14]=new TH1F("xf"," Feynman x_{F};x_{F};# events",100,0.9,1.0);
  HOGENE[15]=new TH1F("cmass","central inv. mass;[GeV];# events",100,0.0,5.0);
  HOGENE[16]=new TH1F("pmt","central tranversal momentum difference mag.;|p_{3t}-p_{4t}| [GeV/c];# events",100,0.0,5.0);
  HOGENE2[0]=new TH2F("ptprotot","total transverse momentum 1+2 vs. av. t ",30,0.0,0.1,60,0.0,0.8);
  HOGENE2[1]=new TH2F("coll12vst","1+2 collinearity vs cmass ",30,0.0,2.0,60,0.0,6.0);
  HOGENE2[2]=new TH2F("azimmass","azimuthal correlation vs. cmass",20,0.0,2.0,60,0.0,180.0);
  HOGENE2[3]=new TH2F("rmassvt","central mass aproximation  vs. av. t",30,0.0,0.1,60,0.95,3.0);

 
 
 
 
  for( int i=0;i<17;i++ ) HOGENE[i]->Sumw2();
  
  
  hEtaPhi = new TH2F("EtaPhi", "Eta vs Phi", 50,0,0, 50,0,0);
  hEtaPhi->SetBit(TH1::kCanRebin);
  hEtaPhi->SetXTitle("#eta");
  hEtaPhi->SetYTitle("#phi"); 
  
};

////////////////////////////////////////////////////////////////////////
void THistogram4::FillHistograms( void )
{

  double wgt = EventWeight;
  
  HOGENE[0]->Fill(pf[3].Pt(),wgt);
  HOGENE[0]->Fill(pf[4].Pt(),wgt);
  HOGENE[1]->Fill(pf[1].Pt(),wgt);
  HOGENE[1]->Fill(pf[2].Pt(),wgt);
  HOGENE[2]->Fill(prap[3],wgt);
  HOGENE[2]->Fill(prap[4],wgt);
  HOGENE[3]->Fill(prap[1],wgt);
  HOGENE[3]->Fill(prap[2],wgt);
  HOGENE[4]->Fill(dazim12,wgt);
  HOGENE[5]->Fill(coll12,wgt);
       HOGENE[6]->Fill(t1,wgt); 
       HOGENE[7]->Fill(t2,wgt); 
       HOGENE[8]->Fill(ta,wgt);
       HOGENE[9]->Fill(tb,wgt);        
       HOGENE[10]->Fill(log10(s13),wgt);        
       HOGENE[11]->Fill(log10(s24),wgt);
       HOGENE[12]->Fill(log10(s14),wgt);
       HOGENE[13]->Fill(log10(s23),wgt);
       HOGENE[14]->Fill(dz1,wgt);
       HOGENE[14]->Fill(dz2,wgt);
       HOGENE[15]->Fill(cmass,wgt);
       HOGENE[16]->Fill( (pf[3]-pf[4]).Pt(),wgt);
       HOGENE2[0]->Fill(tav,ptprotot,wgt);
       HOGENE2[1]->Fill(cmass,coll12,wgt);
       HOGENE2[2]->Fill(cmass,dazim12,wgt);
       HOGENE2[3]->Fill(tav,rmass,wgt);
       
  
  for( int i = 1; i < 5; i++ )
		hEtaPhi->Fill( pf[i].Eta(), pf[i].Phi(), wgt );
       
       

};

////////////////////////////////////////////////////////////////////////
void THistogram4::WriteHistograms( double xsection )
{
		
	//create differential cross-section histograms
	//transform mb to nb
	CreateXSectionHistograms( 1e6 * xsection );
	
  
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
       HOGENE[6]->Draw();
       c22->cd(2);
       HOGENE[7]->Draw();
       c22->cd(3);
       HOGENE[8]->Draw();
       c22->cd(4);
       HOGENE[9]->Draw();
     c22->Update();
     ps->NewPage();
       c22->cd(1);
       HOGENE[4]->Draw();
       c22->cd(2);
       HOGENE[5]->Draw();
       c22->cd(3);
       HOGENE[14]->Draw();
       c22->cd(4);
       HOGENE[15]->Draw();
     c22->Update();
     ps->NewPage();
       c22->cd(1);
       HOGENE2[0]->Draw();
       c22->cd(2);
       HOGENE2[1]->Draw();
       c22->cd(3);
       HOGENE2[2]->Draw();
       c22->cd(4);
       HOGENE2[3]->Draw();
     c22->Update();
     
     ps->NewPage();
     c22->cd(1);
     HOGENE[16]->Draw();
     c22->Update();
     
     ps->NewPage();
     c22->cd(1);
     hEtaPhi->Draw();
     c22->Update();
     
     //write x-section histograms
       ps->NewPage();
       c22->cd(1);
       XHOGENE[0]->Draw();
       c22->cd(2);
       XHOGENE[1]->Draw();
       c22->cd(3);
       XHOGENE[2]->Draw();
       c22->cd(4);
       XHOGENE[3]->Draw();
	   c22->Update();
	   ps->NewPage();
       c22->cd(1);
       XHOGENE[6]->Draw();
       c22->cd(2);
       XHOGENE[7]->Draw();
       c22->cd(3);
       XHOGENE[8]->Draw();
       c22->cd(4);
       XHOGENE[9]->Draw();
       c22->Update();
       ps->NewPage();
       c22->cd(1);
       XHOGENE[4]->Draw();
       c22->cd(2);
       XHOGENE[5]->Draw();
       c22->cd(3);
       XHOGENE[14]->Draw();
       c22->cd(4);
       XHOGENE[15]->Draw();
       c22->Update();
     
       ps->NewPage();
       c22->cd(1);
       XHOGENE[16]->Draw();
       c22->Update();
     
     
     ps->Close();
    
    delete ps;
    delete c22;
    delete c21;
    delete c31;
     
     //ROOT file for histograms
     TFile * tf = new TFile( RootFilename.c_str(),"RECREATE", "Histograms from generator");
     
     for( int i = 0; i < 17; i++ )  
		HOGENE[ i ] -> Write();
     
     for( int i = 0; i < 4; i++ )  
		HOGENE2[ i ] -> Write();
     
     hEtaPhi->Write();
     
    for( int i = 0; i < 17; i++ )  
		XHOGENE[ i ] -> Write();
     
     tf->Write();
     tf->Close();
     
     delete tf;
     
};

////////////////////////////////////////////////////////////////////////
void THistogram4::DeallocateHistograms(void )
{
	
    for( int i = 0; i < 17; i++ )
      delete HOGENE[i];
    
    for( int i = 0; i < 4; i++ )
      delete HOGENE2[i];
    
    delete hEtaPhi ;
    
    return;
};

////////////////////////////////////////////////////////////////////////
void THistogram4::Fill( TEvent * event, double weight  )
{
	//simple check - more elaborate checks could decrease efficeincy
	assert( event->GetNIn() == 2 );
	assert( event->GetNOut() >= 4 );
	
	EventWeight = weight;

	ExtractParticles( event );
	
	Kinematics2to4();
	
	FillHistograms();
	
	Nevents++;
	

	return;
	
};

////////////////////////////////////////////////////////////////////////
void THistogram4::ExtractParticles( TEvent * event )
{
	
	assert( nop >= 4 );
	
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
	pf[4]=event->pf[4];

};

////////////////////////////////////////////////////////////////////////
void THistogram4::Kinematics2to4(void )
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
    TLorentzVector ptot = pf[1] + pf[2] + pf[3] + pf[4];
    double stot = ptot.M2();
    t1 = ( pb1 - pf[1] ).M2();
    t2 = ( pb2 - pf[2] ).M2();
    ta = ( pb1 - pf[1] - pf[3] ).M2();
    tb = ( pb1 - pf[1] - pf[4] ).M2();
    s13 = ( pf[1] + pf[3] ).M2();
    s24 = ( pf[2] + pf[4] ).M2();
    s14 = ( pf[1] + pf[4] ).M2();
    s23 = ( pf[2] + pf[3] ).M2();
    cmass2 = ( pf[3] + pf[4] ).M2();
    cmass = sqrt(cmass2);
    dazim34 = ( azim[3] - azim[4] );
    dazim12 = ( azim[1] - azim[2] );
    dazim34 = r2a*fdelta( dazim34 );
    dazim12 = r2a*fdelta( dazim12 );
    coll12 = 1000.0 * ( M_PI - pf[1].Angle( pf[2].Vect()));//collinearity between protons
    dz1 = pf[1].Pz()/inp1;
    dz2 = pf[2].Pz()/inp1;
    ptprotot=(pf[1]+pf[2]).Pt();
    double mapprox2=stot*(1.0-abs(dz1))*(1.0-abs(dz2));//approximate central mass
      rmass=sqrt(mapprox2)/cmass;
      tav=abs((t1+t2)/2.0);
	
};

////////////////////////////////////////////////////////////////////////
double THistogram4::fdelta(double dazim)
{
	double x = fmod( dazim + M_PI, 2.0*M_PI );
    if (x < 0)
        x += 2.0*M_PI;
    return x - M_PI;
    
};

////////////////////////////////////////////////////////////////////////
void THistogram4::CreateXSectionHistograms( double xsection )
{
	
	//if differential cross-section histograms were created
	if( xsectionHistogramsCreated == true )
	{
		for( int i = 0; i < 17; i++ )
			delete XHOGENE[i];
	}
	
	
	xsectionHistogramsCreated = true;
			
	double binWidth = 0.0;
	double x = 0.0;
	double val = 0.0;
	
	double tmp = 0.0;
	
	//beginning part not modified
	for( int i = 0; i < 5; i++ )
	{
		
		TH1F * tmpHist = (TH1F*) HOGENE[i]->Clone();
		TH1F * tmpHist2 = (TH1F*) HOGENE[i]->Clone();
		tmpHist->Reset("ice");
		tmpHist->GetSumw2()->Set(0);
		//tmpHist->Sumw2();
		
		int nbinsx = tmpHist2->GetNbinsX();
		
		for( int j = 1; j < nbinsx; j++ )
		{
		
			x = tmpHist2->GetBinCenter( j );
			binWidth = tmpHist->GetXaxis()->GetBinWidth( tmpHist->GetXaxis()->FindBin( x ) );
			
			val = tmpHist2->GetBinContent( j );
			
			tmp = val/( binWidth );
			
			tmpHist->Fill( x, tmp );
	
		}
		
		XHOGENE[i] = tmpHist;
	}
	
	
	
	//first set of plots
	for( int i = 5; i < 10; i++)
	{
		
		TH1F * tmpHist = (TH1F*) HOGENE[i]->Clone();
		TH1F * tmpHist2 = (TH1F*) HOGENE[i]->Clone();
		tmpHist->Reset("ice");
		tmpHist->GetSumw2()->Set(0);
		//tmpHist->Sumw2();
		
		//multiply by luminosity
		tmpHist2->Scale( 1.0/double( Nevents ) );
		tmpHist2->Scale( xsection );
		
		int nbinsx = tmpHist2->GetNbinsX();
		
		for( int j = 1; j < nbinsx; j++ )
		{
		
			x = tmpHist2->GetBinCenter( j );
			binWidth = tmpHist->GetXaxis()->GetBinWidth( tmpHist->GetXaxis()->FindBin( x ) );
			
			val = tmpHist2->GetBinContent( j );
			
			tmp = val/( binWidth );
			
			tmpHist->Fill( x, tmp );
	
		}
		
		XHOGENE[i] = tmpHist;
	}
	
	
	//middle part not modified
	for( int i = 10; i < 14; i++)
	{
		
		TH1F * tmpHist = (TH1F*) HOGENE[i]->Clone();
		TH1F * tmpHist2 = (TH1F*) HOGENE[i]->Clone();
		tmpHist->Reset("ice");
		tmpHist->GetSumw2()->Set(0);
		//tmpHist->Sumw2();
		
		int nbinsx = tmpHist2->GetNbinsX();
		
		for( int j = 1; j < nbinsx; j++ )
		{
		
			x = tmpHist2->GetBinCenter( j );
			binWidth = tmpHist->GetXaxis()->GetBinWidth( tmpHist->GetXaxis()->FindBin( x ) );
			
			val = tmpHist2->GetBinContent( j );
			
			tmp = val/( binWidth );
			
			tmpHist->Fill( x, tmp );
	
		}
		
		XHOGENE[i] = tmpHist;
	}
	
	//second set of plots
	for( int i = 14; i < 17; i++)
	{
		
		TH1F * tmpHist = (TH1F*) HOGENE[i]->Clone();
		TH1F * tmpHist2 = (TH1F*) HOGENE[i]->Clone();
		tmpHist->Reset("ice");
		tmpHist->GetSumw2()->Set(0);
		
		//multiply by luminosity
		tmpHist2->Scale( 1.0/double( Nevents ) );
		tmpHist2->Scale( xsection );
		
		int nbinsx = tmpHist2->GetNbinsX();
		
		for( int j = 1; j < nbinsx; j++ )
		{
			
			x = tmpHist2->GetBinCenter( j );
			binWidth = tmpHist->GetXaxis()->GetBinWidth( tmpHist->GetXaxis()->FindBin( x ) );
			
			val = tmpHist2->GetBinContent( j );
			
			tmp = val/( binWidth );
			
			tmpHist->Fill( x, tmp );
	
		}
		
		XHOGENE[i] = tmpHist;
	}
	
	
		XHOGENE[5]->GetXaxis()->SetTitle("#phi_{12} [mrad]");
		XHOGENE[5]->GetYaxis()->SetTitle("#frac{d#sigma}{d#phi_{12}}[ nbarn / mrad]");

		XHOGENE[6]->GetXaxis()->SetTitle("t_{1} [GeV^{2}]");
		XHOGENE[6]->GetYaxis()->SetTitle("#frac{d#sigma}{dt_{1}}[ nbarn / GeV^{2}]");

 		XHOGENE[7]->GetXaxis()->SetTitle("t_{2} [GeV^{2}]");
		XHOGENE[7]->GetYaxis()->SetTitle("#frac{d#sigma}{dt_{2}}[ nbarn / GeV^{2}]");

 		XHOGENE[8]->GetXaxis()->SetTitle("t_{a} [GeV^{2}]");
		XHOGENE[8]->GetYaxis()->SetTitle("#frac{d#sigma}{dt_{a}}[ nbarn / GeV^{2}]");
		
 		XHOGENE[9]->GetXaxis()->SetTitle("t_{b} [GeV^{2}]");
		XHOGENE[9]->GetYaxis()->SetTitle("#frac{d#sigma}{dt_{b}}[ nbarn / GeV^{2}]");
		
 		XHOGENE[14]->GetXaxis()->SetTitle("x_{F}");
		XHOGENE[14]->GetYaxis()->SetTitle("#frac{d#sigma}{dx_{F}}");		

 		XHOGENE[15]->GetXaxis()->SetTitle("M}[GeV]");
		XHOGENE[15]->GetYaxis()->SetTitle("#frac{d#sigma}{dM_{CM}}[ nbarn / GeV]");		
		
		XHOGENE[16]->GetXaxis()->SetTitle("|p_{3t}-p_{4t}| [GeV/c]");
		XHOGENE[16]->GetYaxis()->SetTitle("#frac{d#sigma}{d|p_{3t}-p_{4t}|}[ nbarn c/ GeV]");		


	return;
	
};
