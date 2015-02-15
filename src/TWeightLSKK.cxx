/***********************************************************************
*                          Class TWeightLSKK                                 
*                                                                           
* Returns weight corresponding to matrix element for process pp-->ppK+K- 
* Lebiedowicz & Szczurek Phys.Rev.D85(2012)014026        
*		       
*  									       
***********************************************************************/
#include"TWeightLSKK.h"

////////////////////////////////////////////////////////////////////////

const double TWeightLSKK::PI =  M_PI;


////////////////////////////////////////////////////////////////////////
TWeightLSKK::TWeightLSKK() : Imag( 0.0, 1.0 )
{
  
	ConfigReader = ConfigReaderSubsystem::Instance();
	
	PDGDatabese = TDatabasePDG::Instance();
  
	Logger = LoggerSubsystem::Instance();
  
	ReadConfigFile( string("Generator.dat") );
  
	mProton = PDGDatabese->GetParticle( string("proton").c_str() )->Mass();
	mPion = PDGDatabese->GetParticle( string("pi+").c_str() )->Mass();
	mKPlus = PDGDatabese->GetParticle( string("K+").c_str() )->Mass();
  
	//allocate tables:
	azim  = new double [nop+1];
	theta = new double [nop+1];
	pt    = new double [nop+1];
	prap  = new double [nop+1];
  
	pf = new TLorentzVector [nop+1];
  
  
}

TWeightLSKK::~TWeightLSKK(void)
{

 
   delete [] azim;
   delete [] theta;
   delete [] pt;
   delete [] prap;
 
   delete [] pf;
 	
  
};

////////////////////////////////////////////////////////////////////////
void TWeightLSKK::ReadConfigFile( const string & filename )
{
	//Parse model file
	
	//Parse variables: setting up model variables	
	reggeize =  ConfigReader->GetBoolValue( "reggeize" );
	
	ConfigReader->GetDoubleArray( rintc, "rintc", 5 );
	ConfigReader->GetDoubleArray( rslope, "rslope", 5 );
	ConfigReader->GetDoubleArray( CNN, "CNN", 5 );
	ConfigReader->GetDoubleArray( CKN, "CKN", 5 );
	ConfigReader->GetDoubleArray( CKK, "CKK", 5 );
	ConfigReader->GetDoubleArray( B0, "B0", 5 );
	ConfigReader->GetComplexDoubleArray( eta, "eta", 5 );
	S0 = ConfigReader->GetDoubleValue( "S0" );
	aa = ConfigReader->GetDoubleValue( "aa" );
	W0 = ConfigReader->GetDoubleValue( "W0" );
	sigma_tot = ConfigReader->GetDoubleValue( "sigma_tot" );
	B_elastic = ConfigReader->GetDoubleValue( "B_elastic" );
	Beff = ConfigReader->GetDoubleValue( "Beff" );
	Lambda2 = ConfigReader->GetDoubleValue( "Lambda2" );
	
	LambdaInt = ConfigReader->GetDoubleValue( "LambdaInt" );
	alphaKPrim = ConfigReader->GetDoubleValue( "alphaKPrim" );
	
	//Parse remaining data from generator
	nop = ConfigReader->GetIntValue( "nop" );
	tecm =  ConfigReader->GetDoubleValue( "tecm" );
	
	
	//Log configuration
	string ConfigLog = ConfigReader->GetStringValue( string("ConfigLogFile") );
	TLog * ConfigLogger = new TLog();
	ConfigLogger->setLogTime( false );
	ConfigLogger->openLogFile( ConfigLog );

	ConfigLogger->Write() << "# " << "-----------" << endl;
	ConfigLogger->Write() << "# " << "TWeightLSKK " << endl;
	ConfigLogger->Write() << "# " << "-----------" << endl;
	
	ConfigLogger->Write() << "tecm = " << tecm << endl;
	ConfigLogger->Write() << "nop = " << nop << endl;
	
	
	ConfigLogger->Write() << "reggeize = " << reggeize << endl;
	
	for( int i = 0; i < 5; i++ )
	{
		ConfigLogger->Write() << "rintc" << i << " = " << rintc[i] << endl; 
	}
	for( int i = 0; i < 5; i++ )
	{
		ConfigLogger->Write() << "rslope" << i << " = " << rslope[i] << endl; 
	}
	for( int i = 0; i < 5; i++ )
	{
		ConfigLogger->Write() << "CNN" << i << " = " << CNN[i] << endl; 
	}
	for( int i = 0; i < 5; i++ )
	{
		ConfigLogger->Write() << "CKN" << i << " = " << CKN[i] << endl; 
	}
	for( int i = 0; i < 5; i++ )
	{
		ConfigLogger->Write() << "CKK" << i << " = " << CKK[i] << endl; 
	}
	for( int i = 0; i < 5; i++ )
	{
		ConfigLogger->Write() << "B0" << i << " = " << B0[i] << endl; 
	}
	for( int i = 0; i < 5; i++ )
	{
		ConfigLogger->Write() << "eta" << i << " = " << eta[i] << endl; 
	}
	ConfigLogger->Write() << "S0 = " << S0 << endl; 
	ConfigLogger->Write() << "aa = " << aa << endl;
	ConfigLogger->Write() << "W0 = " << W0 << endl;
	ConfigLogger->Write() << "sigma_tot = " << sigma_tot << endl;
	ConfigLogger->Write() << "B_elastic = " << B_elastic << endl;
	ConfigLogger->Write() << "Beff = " << Beff << endl;
	ConfigLogger->Write() << "Lambda2 = " << Lambda2 << endl;
	ConfigLogger->Write() << "LambdaInt = " << LambdaInt << endl;
	ConfigLogger->Write() << "alphaKPrim = " << alphaKPrim << endl;
	
	delete ConfigLogger;

	
	return;
};

////////////////////////////////////////////////////////////////////////
void TWeightLSKK::Kinematics2to4(TEvent * event)
{
	assert( nop >= 4 );
	

	pb1 = event->pb[1];
	pb2 = event->pb[2];
	
	pf[1]=event->pf[1];
	pf[2]=event->pf[2];
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
double TWeightLSKK::getAcceptance( TEvent * event )
{      
       return( 0 );
};

////////////////////////////////////////////////////////////////////////
bool TWeightLSKK::getIsAccepted( TEvent * event )
{
	bool isAccepted = true;
	
      return isAccepted;
};

////////////////////////////////////////////////////////////////////////
double TWeightLSKK::fdelta(double dazim)
{
	double x = fmod( dazim + M_PI, 2.0*M_PI );
    if (x < 0)
        x += 2.0*M_PI;
    return x - M_PI;
    
};
////////////////////////////////////////////////////////////////////////
complex<double> TWeightLSKK::MatrixElement( void )
{
	const double mpi2 = mPion * mPion;
    complex<double> ME;
    //matrix element
    ME=FFJT(ta) * FFJT(ta) * reggeamp(-1,t1,s13) * reggeamp(1,t2,s24) / (ta-mpi2) +
	  FFJT(tb) * FFJT(tb) * reggeamp(1,t1,s14) * reggeamp(-1,t2,s23) / (tb-mpi2);
	
	return( ME );
	
};

////////////////////////////////////////////////////////////////////////
complex<double> TWeightLSKK::getMatrixElement( TEvent * event )
{
	SetEvent( event );

	
	return( MatrixElement() );
	
};

////////////////////////////////////////////////////////////////////////
double TWeightLSKK::FJT(void)
{
    
    complex<double> ME;
      
    //Matrix element with reggeized propagator  
	TLorentzVector pi1(pf[3]);
	TLorentzVector pi2(pf[4]);
	
	double sPi = ( pi1 + pi2 ).M2();
	
	//matrix element
	if( reggeize == true )
    {   
		const double mKPlus2 = mKPlus * mKPlus;
		 
        ME=FFJT(ta) * FFJT(ta) * reggeamp(-1,t1,s13) * reggeamp(1,t2,s24) / (ta-mKPlus2) +
		FFJT(tb) * FFJT(tb) * reggeamp(1,t1,s14) * reggeamp(-1,t2,s23) / (tb-mKPlus2);
	 }
	 else
	 { 
		ME= reggeamp(-1,t1,s13) * reggeamp(1,t2,s24) * ReggeizedKaon( ta, sPi, ta ) +
		reggeamp(1,t1,s14) * reggeamp(-1,t2,s23) * ReggeizedKaon( tb, sPi, tb );
     }
     
    	  
    double MESQ=abs(ME) * abs(ME);

    return MESQ;
  
};


////////////////////////////////////////////////////////////////////////
complex <double> TWeightLSKK::A( int particle, double s, double t )
{
  
  double alpha = rintc[particle] + t * rslope[particle];
  
  complex<double> tmp = s * eta[particle] * CKN[particle] * exp( log( s / S0 ) * ( alpha - 1.0 ) + B0[particle] * t / 2.0);
  
  return( tmp );
  
}

////////////////////////////////////////////////////////////////////////
complex<double> TWeightLSKK::reggeamp (int isign,double tt,double ss)
{

   complex <double> temp( 0.0, 0.0 );
   
   // K+- p -> K+- p
   temp = A( POMERON, ss, tt ) + A( F2, ss, tt ) + A( A2, ss, tt ) + double( isign ) * ( A( OMEGA, ss, tt ) + A( RHO, ss, tt ) );
   
   return temp;
  
};

////////////////////////////////////////////////////////////////////////
double TWeightLSKK::FFJT (double tt)
{
     double temp=exp( ( tt - mKPlus * mKPlus ) / Lambda2 );
     return temp;
};

////////////////////////////////////////////////////////////////////////
double TWeightLSKK::alphaK( double t )
{
	double tmp = alphaKPrim * ( t - mKPlus * mKPlus );
	
	return( tmp );
};

////////////////////////////////////////////////////////////////////////
double TWeightLSKK::betaM( double s )
{
	double tmp = exp( -( s - 4.0 * mKPlus * mKPlus ) / ( LambdaInt * LambdaInt ) );

	return( tmp );
};

////////////////////////////////////////////////////////////////////////
double TWeightLSKK::betaR( double s )
{
	double tmp = 1.0 - betaM( s );

	return( tmp );
};

////////////////////////////////////////////////////////////////////////
complex<double> TWeightLSKK::PKReggeized( double t, double s , double tt )
{  
	complex<double> factor1 = alphaKPrim * ( 1.0 + exp( - Imag * M_PI * alphaK( t ) ) ) /  2.0 ;
	
	complex<double> factor2 = -exp( lgamma( - alphaK( t ) ) + log( s / S0 ) *  alphaK( t ) + 2.0 * ( tt - mKPlus * mKPlus ) / Lambda2 );
	
	complex<double> factor3 = 1.0;
	
	complex<double> tmp = factor1 * factor2 * factor3;

	return( tmp );

};

////////////////////////////////////////////////////////////////////////
complex<double> TWeightLSKK::ReggeizedKaon( double t, double s, double tt )
{
 
	complex<double> tmp = betaM( s ) * ( 1.0 / ( t - mKPlus * mKPlus ) ) + betaR( s ) * PKReggeized( t, s, tt );

	
	return( tmp );
};

////////////////////////////////////////////////////////////////////////
void TWeightLSKK::SetEvent( TEvent * event, double eventWeight )
{
	
	Kinematics2to4( event );
	
};

////////////////////////////////////////////////////////////////////////
double TWeightLSKK::GetWeight( TEvent * event, double eventWeight )
{
	
	Kinematics2to4( event );
	double weight = GetWeight( eventWeight );
	
	return( weight );

};

////////////////////////////////////////////////////////////////////////
double TWeightLSKK::GetWeight( double eventWeight )
{
	double wt = eventWeight;
	
	double wtf = wt*FJT();
	
	
	return wtf;
	
};
