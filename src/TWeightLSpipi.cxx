/***********************************************************************
*                         TWeightLSpipi                                 
*                                                                           
* Returns weight corresponding to matrix element for process pp-->pppi+pi- 
* Lebiedowicz & Szczurek Phys.Rev.D81(2010)036003        
*		       
*  									       
***********************************************************************/
#include"TWeightLSpipi.h"

////////////////////////////////////////////////////////////////////////

const double TWeightLSpipi::PI =  M_PI;

////////////////////////////////////////////////////////////////////////
TWeightLSpipi::TWeightLSpipi()
{
  
	ConfigReader = ConfigReaderSubsystem::Instance();
	
	PDGDatabese = TDatabasePDG::Instance();
	
	Logger = LoggerSubsystem::Instance();
  
	ReadConfigFile( string("Generator.dat") );
  
	mProton = PDGDatabese->GetParticle( string("proton").c_str() )->Mass();
	mPion = PDGDatabese->GetParticle( string("pi+").c_str() )->Mass();
  
  
	//allocate tables:
	azim  = new double [nop+1];
	theta = new double [nop+1];
	pt    = new double [nop+1];
	prap  = new double [nop+1];
  
	pf = new TLorentzVector [nop+1];
  
  
}

TWeightLSpipi::~TWeightLSpipi(void)
{
	 
	delete [] azim;
	delete [] theta;
	delete [] pt;
	delete [] prap;
  
   
	delete [] pf;
 	
  
};

////////////////////////////////////////////////////////////////////////
void TWeightLSpipi::ReadConfigFile( const string & filename )
{
	//Parse model file
	
	//Parse variables: setting up model variables	
	ConfigReader->GetDoubleArray( rintc, "rintc", 3 );
	ConfigReader->GetDoubleArray( rslope, "rslope", 3 );
	ConfigReader->GetDoubleArray( areal, "areal", 3 );
	ConfigReader->GetDoubleArray( Cregge, "Cregge", 3 );
	ConfigReader->GetDoubleArray( B0, "B0", 3 );
	S0 = ConfigReader->GetDoubleValue( "S0" );
	aa = ConfigReader->GetDoubleValue( "aa" );
	W0 = ConfigReader->GetDoubleValue( "W0" );
	sigma_tot = ConfigReader->GetDoubleValue( "sigma_tot" );
	B_elastic = ConfigReader->GetDoubleValue( "B_elastic" );
	Beff = ConfigReader->GetDoubleValue( "Beff" );
	Lambda2 = ConfigReader->GetDoubleValue( "Lambda2" );
	
	//Parse remaining data from generator
	nop = ConfigReader->GetIntValue( "nop" );
	tecm =  ConfigReader->GetDoubleValue( "tecm" );
	
	
	
	//Log configuration
	string ConfigLog = ConfigReader->GetStringValue( string("ConfigLogFile") );
	TLog * ConfigLogger = new TLog();
	ConfigLogger->setLogTime( false );
	ConfigLogger->openLogFile( ConfigLog );

	ConfigLogger->Write() << "# " << "-----------" << endl;
	ConfigLogger->Write() << "# " << "TWeightLSpipi " << endl;
	ConfigLogger->Write() << "# " << "-----------" << endl;
	
	ConfigLogger->Write() << "tecm = " << tecm << endl;
	ConfigLogger->Write() << "nop = " << nop << endl;
	
	
	for( int i = 0; i < 3; i++ )
	{
		ConfigLogger->Write() << "rintc" << i << " = " << rintc[i] << endl; 
	}
	for( int i = 0; i < 3; i++ )
	{
		ConfigLogger->Write() << "rslope" << i << " = " << rslope[i] << endl; 
	}
	for( int i = 0; i < 3; i++ )
	{
		ConfigLogger->Write() << "areal" << i << " = " << areal[i] << endl; 
	}
	for( int i = 0; i < 3; i++ )
	{
		ConfigLogger->Write() << "Cregge" << i << " = " << Cregge[i] << endl; 
	}
	for( int i = 0; i < 5; i++ )
	{
		ConfigLogger->Write() << "B0" << i << " = " << B0[i] << endl; 
	}
	ConfigLogger->Write() << "S0 = " << S0 << endl; 
	ConfigLogger->Write() << "aa = " << aa << endl;
	ConfigLogger->Write() << "W0 = " << W0 << endl;
	ConfigLogger->Write() << "sigma_tot = " << sigma_tot << endl;
	ConfigLogger->Write() << "B_elastic = " << B_elastic << endl;
	ConfigLogger->Write() << "Beff = " << Beff << endl;
	ConfigLogger->Write() << "Lambda2 = " << Lambda2 << endl;

	delete ConfigLogger;

	
	return;
};

////////////////////////////////////////////////////////////////////////
void TWeightLSpipi::Kinematics2to4(TEvent * event)
{
	
	assert( nop >= 4 );
	
	//Transform momentum to the CM where Matrix Element should be calculated
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
	

	for(int i=1;i<nop+1;i++)
    {
		//prap[i]=pf[i].PseudoRapidity();
		prap[i]=pf[i].Rapidity();
		azim[i]=pf[i].Phi();
		theta[i]=pf[i].Theta();
		pt[i]=pf[i].Pt(); 
     }
  
     
      //set fourvectors calculate pseudorapidities, angles... 
      double inp1=sqrt(tecm*tecm/4.0-mProton*mProton);
      double r2a = 180.0 / M_PI;
      
      //fourbody kinematics kinematics
      TLorentzVector ptot=pf[1]+pf[2]+pf[3]+pf[4];
      //double stot=ptot.M2();
      t1=(pb1-pf[1]).M2();
      t2=(pb2-pf[2]).M2();
      ta=(pb1-pf[1]-pf[3]).M2();
      tb=(pb1-pf[1]-pf[4]).M2();
      s13=(pf[1]+pf[3]).M2();
      s24=(pf[2]+pf[4]).M2();
      s14=(pf[1]+pf[4]).M2();
      s23=(pf[2]+pf[3]).M2();
      cmass2=(pf[3]+pf[4]).M2();
      cmass=sqrt(cmass2);
	  dazim34 = ( azim[3] - azim[4] );
	  dazim12 = ( azim[1] - azim[2] );
	  dazim34 = r2a*fdelta( dazim34 );
	  dazim12 = r2a*fdelta( dazim12 );
      coll12=1000.0*(PI-pf[1].Angle(pf[2].Vect()));//collinearity between protons
      dz1=pf[1].Pz()/inp1;
      dz2=pf[2].Pz()/inp1;

};
////////////////////////////////////////////////////////////////////////
double TWeightLSpipi::getAcceptance( TEvent * event )
{      
       return( 0 );
};

////////////////////////////////////////////////////////////////////////
bool TWeightLSpipi::getIsAccepted( TEvent * event )
{
	bool isAccepted = true;
	
      return isAccepted;
};

////////////////////////////////////////////////////////////////////////
double TWeightLSpipi::fdelta(double dazim)
{
	double x = fmod( dazim + M_PI, 2.0*M_PI );
    if (x < 0)
        x += 2.0*M_PI;
    return x - M_PI;
    
};

////////////////////////////////////////////////////////////////////////
complex<double> TWeightLSpipi::MatrixElement( void )
{
	const double mpi2 = mPion * mPion;
    complex<double> ME;
    //matrix element
    ME=FFJT(ta) * FFJT(ta) * reggeamp(-1,t1,s13) * reggeamp(1,t2,s24) / (ta-mpi2) +
	  FFJT(tb) * FFJT(tb) * reggeamp(1,t1,s14) * reggeamp(-1,t2,s23) / (tb-mpi2);
	
	return( ME );
	
};

////////////////////////////////////////////////////////////////////////
complex<double> TWeightLSpipi::getMatrixElement( TEvent * event )
{
	SetEvent( event );

	
	return( MatrixElement() );
	
};

////////////////////////////////////////////////////////////////////////
double TWeightLSpipi::FJT(void)
{
       complex<double> ME = MatrixElement();
       
       double MESQ=abs(ME) * abs(ME);

       return MESQ;
  
};

////////////////////////////////////////////////////////////////////////
complex<double> TWeightLSpipi::reggeamp (int isign,double tt,double ss)
{
  //using namespace std;
   complex<double> Rgg;
   double alpha[3],Bslope[3];
   complex <double> phase[3];
   complex <double> iunit( 0.0, 1.0 );
   double sgn = 1.0;
   complex <double> temp( 0.0, 0.0 );
   for( int ii=0; ii<3; ii++ )
   {
      sgn = 1.0;
      if( ii==2 )sgn = isign * 1.0;
      alpha[ii] = rintc[ii] + tt * rslope[ii];
      Bslope[ii] = B0[ii];
      phase[ii] = areal[ii] + iunit;
      Rgg = sgn * ss * phase[ii] * Cregge[ii] * exp( log( ss / S0 ) * ( alpha[ii] - 1.0 ) + Bslope[ii] * tt / 2.0);

	  double W = sqrt(ss);
      double cutoff = 0.0;
      if(W>5.0)
      {
        cutoff = 1.0;
      }
      else
      {
        double expcut=exp((W-W0)/aa);
        cutoff=expcut/(1.0+expcut);
       
      }
     
      temp=temp+Rgg*cutoff;
    }

       return temp;
  
};

////////////////////////////////////////////////////////////////////////
double TWeightLSpipi::FFJT (double tt)
{
       double temp=exp( ( tt - mPion * mPion ) / Lambda2 );
       return temp;
};

////////////////////////////////////////////////////////////////////////
void TWeightLSpipi::SetEvent( TEvent * event, double eventWeight )
{
	
	Kinematics2to4( event );
	
};

////////////////////////////////////////////////////////////////////////
double TWeightLSpipi::GetWeight( TEvent * event, double eventWeight )
{
	
	Kinematics2to4( event );
	double weight = GetWeight( eventWeight );
	
	return( weight );

};

////////////////////////////////////////////////////////////////////////
double TWeightLSpipi::GetWeight( double eventWeight )
{
	double wt = eventWeight;
	
	double wtf = wt*FJT();
	
	
	return wtf;
	
};
