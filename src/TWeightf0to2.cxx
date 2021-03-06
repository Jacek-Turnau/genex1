/***********************************************************************
*                         TWeightf0to2                                 
*                                                                           
* Returns weight corresponding to matrix element for process pp-->ppf0(1500) 
* Machado        
*		       
*  									       
***********************************************************************/
#include"TWeightf0to2.h"
#include"TEventMaker2toN.h"

////////////////////////////////////////////////////////////////////////

const double TWeightf0to2::PI =  M_PI;

////////////////////////////////////////////////////////////////////////
TWeightf0to2::TWeightf0to2()
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

TWeightf0to2::~TWeightf0to2(void)
{
	 
	delete [] azim;
	delete [] theta;
	delete [] pt;
	delete [] prap;
  
   
	delete [] pf;
 	
  
};

////////////////////////////////////////////////////////////////////////
void TWeightf0to2::ReadConfigFile( const string & filename )
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
	mass_min = ConfigReader->GetDoubleValue( "TEventMaker2toN::mass_min" );
	mass_max = ConfigReader->GetDoubleValue( "TEventMaker2toN::mass_max" );
	
	//Parse remaining data from generator
	nop = ConfigReader->GetIntValue( "nop" );
	tecm =  ConfigReader->GetDoubleValue( "tecm" );
	
	
	
	//Log configuration
	string ConfigLog = ConfigReader->GetStringValue( string("ConfigLogFile") );
	TLog * ConfigLogger = new TLog();
	ConfigLogger->setLogTime( false );
	ConfigLogger->openLogFile( ConfigLog );

	ConfigLogger->Write() << "# " << "-----------" << endl;
	ConfigLogger->Write() << "# " << "TWeightf0to2 " << endl;
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
void TWeightf0to2::Kinematics2to4(TEvent * event)
{
	
	assert( nop >= 4 );
	
	//Extract particles
	pb1 = event->pb[1];
	pb2 = event->pb[2];
	pf[3]=event->pf[3];
	pf[4]=event->pf[4];
	pf[1]=event->pf[1];
	pf[2]=event->pf[2];

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
      //double r2a = 180.0 / M_PI;
      //fourbody kinematics kinematics
      t1=(pb1-pf[1]).M2();
      t2=(pb2-pf[2]).M2();
      s134=(pf[1]+pf[3]+pf[4]).M2();
      s234=(pf[2]+pf[3]+pf[4]).M2();
      cmass2=(pf[3]+pf[4]).M2();
      cmass=sqrt(cmass2);
      dz1=pf[1].Pz()/inp1;
      dz2=pf[2].Pz()/inp1;

};
////////////////////////////////////////////////////////////////////////
double TWeightf0to2::getAcceptance( TEvent * event )
{      
       return( 0 );
};

////////////////////////////////////////////////////////////////////////
bool TWeightf0to2::getIsAccepted( TEvent * event )
{
	bool isAccepted = true;
	
      return isAccepted;
};

////////////////////////////////////////////////////////////////////////
double TWeightf0to2::fdelta(double dazim)
{
	double x = fmod( dazim + M_PI, 2.0*M_PI );
    if (x < 0)
        x += 2.0*M_PI;
    return x - M_PI;
    
};

////////////////////////////////////////////////////////////////////////
complex<double> TWeightf0to2::MatrixElement( void )
{
    complex<double> ME;
    //matrix element
    //ME=pow(tecm,2)*reggeamp(t2,s134) * reggeamp(t1,s234)*BWamp()*pow(mres*gres,0.5)/pow(mres,2);
    ME=pow(tecm,2)*reggeamp(t2,s134) * reggeamp(t1,s234)*pow(mres*gres,0.5)/pow(mres,2);
	return( ME );
	
};

////////////////////////////////////////////////////////////////////////
complex<double> TWeightf0to2::getMatrixElement( TEvent * event )
{
	SetEvent( event );

	
	return( MatrixElement() );
	
};


////////////////////////////////////////////////////////////////////////
complex<double> TWeightf0to2::reggeamp (double tt,double ss)
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

///////////////////////////////////////////////////////
double TWeightf0to2::BW2()
{
  // modulus square of relativistic Breight Wigner for scalar f0
   double mres2=pow(mres,2.0);
   double gres2=pow(gres,2.0);
   double BW= mres*gres/(pow((mres2-cmass2),2)+gres2*mres2);
  return BW;
}
////////////////////////////////////////////////////////////////////////
complex<double> TWeightf0to2::BWamp()
{
  // modulus square of relativistic Breight Wigner for scalar f0 with integral normalized to 1
   complex <double> iunit( 0.0, 1.0 );
   double mres2=pow(mres,2.0);
   double gres2=pow(gres,2.0);
   complex<double> BW= pow(mres*gres,0.5)*((mres2-cmass2)+iunit*mres*gres)/(pow((mres2-cmass2),2)+gres2*mres2);
   BW=BW/pow(NormBW(),0.5);
   return BW;
}
////////////////////////////////////////////////////////////////////////
double TWeightf0to2::NormBW()
{
  // BW2 integral for range mass2_min - mass2_max
   double mres2=pow(mres,2.0);
   double s1 = pow(mass_min,2);
   double s2 = pow(mass_max,2);
   double Norm=atan((mres2-s1)/mres/gres)-atan((mres2-s2)/mres/gres) ;
  return Norm;
}
////////////////////////////////////////////////////////////////////////
void TWeightf0to2::SetEvent( TEvent * event, double eventWeight )
{
	
	Kinematics2to4( event );
	
};

////////////////////////////////////////////////////////////////////////
double TWeightf0to2::GetWeight( TEvent * event, double eventWeight )
{
	
	Kinematics2to4( event );
	double weight = GetWeight( eventWeight );
	
	return( weight );

};

////////////////////////////////////////////////////////////////////////
double TWeightf0to2::kcms(double s, double m1, double m2)
{
  // momentum in cms of systme of particles with masses m1, m2
  double m12=m1*m1;
  double m22=m2*m2;
  double temp= sqrt(pow((s-m12-m22),2)-4.0*m12*m22)/2.0/sqrt(s);
  return temp;

};
////////////////////////////////////////////////////////////////////////
double TWeightf0to2::GetWeight( double eventWeight )
{

        complex<double> ME = MatrixElement();
	double wt = eventWeight;	
	double wtf = wt*norm(ME);
	double wtdecay=kcms(cmass2,pf[3].M(),pf[4].M())/cmass/8.0/pow(M_PI,2);
	wtf=wtf/(pow(mass_max,2)-pow(mass_min,2))/wtdecay;
	/*	
        if(norm(ME) > 0.0){
	cout << " Phs_wt = " << wt << endl;	
	cout << " ME2  = " << norm(ME) << endl;
	cout << " wtf = " << wtf << endl;}*/	
	
	return wtf;
	
};
