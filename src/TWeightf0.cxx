/***********************************************************************
*                         TWeightf0                                
*                                                                           
* Returns weight corresponding to matrix element for process pp-->ppf0(1500) 
* Machado        
*		       
*  									       
***********************************************************************/
#include"TWeightf0.h"

////////////////////////////////////////////////////////////////////////

const double TWeightf0::PI =  M_PI;
///est genex1
////////////////////////////////////////////////////////////////////////
TWeightf0::TWeightf0()
{
  
	ConfigReader = ConfigReaderSubsystem::Instance();
	
	PDGDatabese = TDatabasePDG::Instance();
	
	Logger = LoggerSubsystem::Instance();
  
	ReadConfigFile( string("Generator.dat") );
  
	mProton = PDGDatabese->GetParticle( string("proton").c_str() )->Mass();
  
	//allocate tables:
	azim  = new double [nop+1];
	theta = new double [nop+1];
	pt    = new double [nop+1];
	prap  = new double [nop+1];
  
	pf = new TLorentzVector [nop+1];
  
  
}

TWeightf0::~TWeightf0(void)
{
	 
	delete [] azim;
	delete [] theta;
	delete [] pt;
	delete [] prap;
  
   
	delete [] pf;
 	
  
};

////////////////////////////////////////////////////////////////////////
void TWeightf0::ReadConfigFile( const string & filename )
{
	//Parse model file
	
	//Parse variables: setting up model variables	
	ConfigReader->GetDoubleArray( rintc, "rintc", 2 );
	ConfigReader->GetDoubleArray( rslope, "rslope", 2 );
	ConfigReader->GetDoubleArray( areal, "areal", 2 );
	ConfigReader->GetDoubleArray( Cregge, "Cregge", 2 );
	ConfigReader->GetDoubleArray( B0, "B0", 2 );
	Bglu = ConfigReader->GetDoubleValue( "Bglu" );
	aa = ConfigReader->GetDoubleValue( "aa" );
	W0 = ConfigReader->GetDoubleValue( "W0" );
        mres = ConfigReader->GetDoubleValue( "mres" );
        gres = ConfigReader->GetDoubleValue( "gres" );
	
	//Parse remaining data from generator
	nop = ConfigReader->GetIntValue( "nop" );
	tecm =  ConfigReader->GetDoubleValue( "tecm" );
	
	
	
	//Log configuration
	string ConfigLog = ConfigReader->GetStringValue( string("ConfigLogFile") );
	TLog * ConfigLogger = new TLog();
	ConfigLogger->setLogTime( false );
	ConfigLogger->openLogFile( ConfigLog );

	ConfigLogger->Write() << "# " << "-----------" << endl;
	ConfigLogger->Write() << "# " << "TWeightf0 " << endl;
	ConfigLogger->Write() << "# " << "-----------" << endl;
	
	ConfigLogger->Write() << "tecm = " << tecm << endl;
	ConfigLogger->Write() << "nop = " << nop << endl;
	
	
	for( int i = 0; i < 2; i++ )
	{
		ConfigLogger->Write() << "rintc" << i << " = " << rintc[i] << endl; 
	}
	for( int i = 0; i < 2; i++ )
	{
		ConfigLogger->Write() << "rslope" << i << " = " << rslope[i] << endl; 
	}
	for( int i = 0; i < 2; i++ )
	{
		ConfigLogger->Write() << "areal" << i << " = " << areal[i] << endl; 
	}
	for( int i = 0; i < 2; i++ )
	{
		ConfigLogger->Write() << "Cregge" << i << " = " << Cregge[i] << endl; 
	}
	for( int i = 0; i < 2; i++ )
	{
		ConfigLogger->Write() << "B0" << i << " = " << B0[i] << endl; 
	}
	ConfigLogger->Write() << "S0 = " << S0 << endl; 
	ConfigLogger->Write() << "aa = " << aa << endl;
	ConfigLogger->Write() << "W0 = " << W0 << endl;
	ConfigLogger->Write() << "Bglu= " << Bglu << endl;
	ConfigLogger->Write() << "mres= " << mres << endl;
	ConfigLogger->Write() << "gres= " << gres << endl;

	delete ConfigLogger;

	
	return;
};

////////////////////////////////////////////////////////////////////////
void TWeightf0::Kinematics2to3(TEvent * event)
{
	
	assert( nop >= 3 );
	
	//Transform momentum to the CM where Matrix Element should be calculated
	pb1 = event->pb[1];
	pb2 = event->pb[2];
	pf[1]=event->pf[1];
	pf[2]=event->pf[2];
	pf[3]=event->pf[3];

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
      t1=(pb1-pf[1]).M2();
      t2=(pb2-pf[2]).M2();
      s13=(pf[1]+pf[3]).M2();
      s23=(pf[2]+pf[3]).M2();
 	  dazim12 = ( azim[1] - azim[2] );
	  dazim12 = r2a*fdelta( dazim12 );
      coll12=1000.0*(PI-pf[1].Angle(pf[2].Vect()));//collinearity between protons
      dz1=pf[1].Pz()/inp1;
      dz2=pf[2].Pz()/inp1;

};
////////////////////////////////////////////////////////////////////////
double TWeightf0::getAcceptance( TEvent * event )
{      
       return( 0 );
};

////////////////////////////////////////////////////////////////////////
bool TWeightf0::getIsAccepted( TEvent * event )
{
	bool isAccepted = true;
	
      return isAccepted;
};

////////////////////////////////////////////////////////////////////////
double TWeightf0::fdelta(double dazim)
{
	double x = fmod( dazim + M_PI, 2.0*M_PI );
    if (x < 0)
        x += 2.0*M_PI;
    return x - M_PI;
    
};

////////////////////////////////////////////////////////////////////////
complex<double> TWeightf0::MatrixElement( void )
{
    complex<double> ME;
    //matrix element
    ME=pow(tecm,2)*reggeamp(t2,s13) * reggeamp(t1,s23)*pow(mres*gres,0.5)/pow(mres,2);
	return( ME );
	
};

////////////////////////////////////////////////////////////////////////
complex<double> TWeightf0::getMatrixElement( TEvent * event )
{
	SetEvent( event );

	
	return( MatrixElement() );
	
};


////////////////////////////////////////////////////////////////////////
complex<double> TWeightf0::reggeamp (double tt,double ss)
{
  //using namespace std;
   complex<double> Rgg;
   double alpha[2],Bslope[2];
   complex <double> phase[2];
   complex <double> iunit( 0.0, 1.0 );
   complex <double> temp( 0.0, 0.0 );
   for( int ii=0; ii<2; ii++ )
   {
      alpha[ii] = rintc[ii] + tt * rslope[ii];
      Bslope[ii] = B0[ii]+Bglu;//proton formfactor + gluon nonperturbative formfactor
      phase[ii] = areal[ii] + iunit;
      Rgg = phase[ii] * Cregge[ii] * exp( log( pow(tecm,2) /ss) * ( alpha[ii] - 1.0 ) + Bslope[ii] * tt);
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
void TWeightf0::SetEvent( TEvent * event, double eventWeight )
{
	
	Kinematics2to3( event );
	
};

////////////////////////////////////////////////////////////////////////
double TWeightf0::GetWeight( TEvent * event, double eventWeight )
{
	
	Kinematics2to3( event );
	double weight = GetWeight( eventWeight );
	
	return( weight );

};

////////////////////////////////////////////////////////////////////////
double TWeightf0::GetWeight( double eventWeight )
{

        complex<double> ME = MatrixElement();
	double wt = eventWeight;	
	double wtf = wt*norm(ME)*0.38935;// weight in mb
	
	return wtf;
	
};
