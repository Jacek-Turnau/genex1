/***********************************************************************
*                TEventMaker2toN                                 
*                                                                           
* Generate event and weight of phase space for process 2 -> 2 + M, M->3+4+...+N                                            
*
* It bases on TGenPhaseSpace, which bases on GENBOD by F. James 
*            
*   									       
***********************************************************************/


#include"TEventMaker2toN.h"

////////////////////////////////////////////////////////////////////////



//Constants:
const double TEventMaker2toN::PI =  M_PI;



////////////////////////////////////////////////////////////////////////
TEventMaker2toN::TEventMaker2toN() 
{
  
	ConfigReader = ConfigReaderSubsystem::Instance();
	
	Logger = LoggerSubsystem::Instance();
 
	PDGDatabese = TDatabasePDG::Instance();

	ReadConfigFile( string("Generator.dat") );
   
	Decay = new TDecay();
   
   	mA = PDGDatabese->GetParticle(idA)->Mass();
	mB = PDGDatabese->GetParticle(idB)->Mass(); 
  
	///setup CM frame
	EA = (pow(mA,2) - pow(mB,2) + pow(tecm,2))/(2.0 * tecm);
	EB = (-pow(mA,2) + pow(mB,2) + pow(tecm,2))/(2.0 *tecm);
		
	double pzA = sqrt( EA*EA - mA*mA );
	double pzB = -sqrt( EB*EB - mB*mB );
		
	pb1.SetPxPyPzE( 0.0, 0.0, pzA, EA );
	pb2.SetPxPyPzE( 0.0, 0.0, pzB, EB );
  
	//initialize internal structure of the class
	Initialize();
  
};

////////////////////////////////////////////////////////////////////////
TEventMaker2toN::~TEventMaker2toN(void)
{
	
	delete [] mass;
	delete [] idOut;
	delete [] pf;
   
	delete Decay;
   
	if( bIDecay != 0)
	{
		delete DecayIntegral;
	}
  
};

////////////////////////////////////////////////////////////////////////
void TEventMaker2toN::ReadConfigFile( const string & filename )
{
	//Parse variables: setting up model variables	
	
	tecm =  ConfigReader->GetDoubleValue( "tecm" );
	p_min = ConfigReader->GetDoubleValue( "TEventMaker2toN::p_min" );	
	p_max = ConfigReader->GetDoubleValue( "TEventMaker2toN::p_max" );
	y_min = ConfigReader->GetDoubleValue( "TEventMaker2toN::y_min" ); 
	y_max = ConfigReader->GetDoubleValue( "TEventMaker2toN::y_max" ); 
	weightStrategy = ConfigReader->GetIntValue( "weightStrategy" );
	
	bIDecay = ConfigReader->GetIntValue( "TEventMaker2toN::bIDecay" );
 
	///Reading frame setup
	frameType = ConfigReader->GetIntValue( string("frameType") );
	idA = ConfigReader->GetIntValue( string("idA") );
	idB = ConfigReader->GetIntValue( string("idB") );
	tecm = ConfigReader->GetDoubleValue( string("tecm") );
	 
	//read output particles 
	int massIndication = ConfigReader->GetIntValue( string("TEventMaker2toN::massIndication") );

	nop = ConfigReader->GetIntValue( "TEventMaker2toN::nop" );

	//allocate tables
	mass  = new double [nop+1];
	idOut = new int[nop+1];
	pf = new TLorentzVector [nop+1];


	if( massIndication == 1)
	{
		ConfigReader->GetDoubleArray( mass, string("TEventMaker2toN::m"), nop+1, 1 );
		
		for( int i = 1; i < nop+1; i++ )
		{
			idOut[i] = 81; //unkonown particle if mass specified
			//cout << "id= " << idOut[i] << ", m = " << mass[i] << endl;
		}
		
	}
	else if( massIndication == 0 )
	{
		ConfigReader->GetIntArray( idOut, string("TEventMaker2toN::id"), nop+1, 1 );	
		
		for( int i = 1; i < nop+1; i++ )
		{
			mass[i] = PDGDatabese->GetParticle(idOut[i])->Mass();
			//cout << "id= " << idOut[i] << ", m = " << mass[i] << endl;
		}
	}
	else if( massIndication == 2)
	{
		ConfigReader->GetDoubleArray( mass, string("TEventMaker2toN::m"), nop+1, 1 );
		ConfigReader->GetIntArray( idOut, string("TEventMaker2toN::id"), nop+1, 1 );
	}
	

	//read central mass setup
	centralMassOption = ConfigReader->GetIntValue( "TEventMaker2toN::centralMassOption" );
	if(  centralMassOption == 0 )
	{
		double mmin = 0.0;
		for(int i = 3; i < nop+1; i++)
		{
			mmin += mass[i];
		}
	
		double eps = 0.001;
	
		mass_min=mmin+eps;
		//creation of particles from the range [mthr; tecm-mthr]
		mass_max=tecm - mass[1] - mass[2] - eps;
		
	}
	else if( centralMassOption == 1 )
	{
		mass_min = 	ConfigReader->GetDoubleValue( "TEventMaker2toN::mass_min" );
		mass_max = 	ConfigReader->GetDoubleValue( "TEventMaker2toN::mass_max" );
	}
	else if( centralMassOption == 2 )
	{
		double mmin = 0.0;
		for(int i = 3; i < nop+1; i++)
		{
			mmin += mass[i];
		}
	
		double eps = 0.001;
		mass_min=mmin+eps;
		
		mass_max = ConfigReader->GetDoubleValue( "TEventMaker2toN::mass_max" );
	}


	//Log configuration
	string ConfigLog = ConfigReader->GetStringValue( string("ConfigLogFile") );
	TLog * ConfigLogger = new TLog();
	ConfigLogger->setLogTime( false );
	ConfigLogger->openLogFile( ConfigLog );
 
	ConfigLogger->Write() << "# " << "-----------" << endl;
	ConfigLogger->Write() << "# " << "TEventMaker2toN " << endl;
	ConfigLogger->Write() << "# " << "-----------" << endl;
	ConfigLogger->Write() << "TEventMaker2toN::p_min = " << p_min << endl;
	ConfigLogger->Write() << "TEventMaker2toN::p_max = " << p_max << endl;
	ConfigLogger->Write() << "TEventMaker2toN::y_min = " << y_min << endl;
	ConfigLogger->Write() << "TEventMaker2toN::y_max = " << y_max << endl;
	ConfigLogger->Write() << "TEventMaker2toN::nop = " << nop << endl;
	ConfigLogger->Write() << "weightStrategy = " << weightStrategy << endl;
	ConfigLogger->Write() << "frameType = " << frameType << endl;
	ConfigLogger->Write() << "idA = " << idA << endl;
	ConfigLogger->Write() << "idB = " << idB << endl;
	ConfigLogger->Write() << "tecm = " << tecm << endl;
	
	ConfigLogger->Write() << "TEventMaker2toN::bIDecay = " << bIDecay << endl;
	
	ConfigLogger->Write() << "TEventMaker2toN::centralMassOption = " << centralMassOption << endl;
	if(  centralMassOption == 0 )
	{
	}
	else if( centralMassOption == 1 )
	{
		ConfigLogger->Write() << "TEventMaker2toN::mass_min = " << mass_min << endl;
		ConfigLogger->Write() << "TEventMaker2toN::mass_max = " << mass_max << endl;
	}
	else if( centralMassOption == 2 )
	{
		ConfigLogger->Write() << "TEventMaker2toN::mass_max = " << mass_max << endl;
	}
	
	
	ConfigLogger->Write() << "TEventMaker2toN::massIndication = " << massIndication << endl;
	for( int i = 1; i < 5; i++ )
	{
			ConfigLogger->Write() <<  "TEventMaker2toN::id" << i << " = " << idOut[i] << endl;
			ConfigLogger->Write() <<  "TEventMaker2toN::m" << i << " = " << mass[i] << endl;
	}
	
	
	delete ConfigLogger;

	
	return;
};

////////////////////////////////////////////////////////////////////////
void TEventMaker2toN::Initialize()
{ 
	
    //SETTING UP THE SYSTEM OF PARTICLES AND PROCESS
	
	mass2_min=pow(mass_min,2);
	mass2_max=pow(mass_max,2);
	
	//initialize decay integral calculation
	if( bIDecay != 0)
	{
		Logger->Write( logINFO ) << "Initialize TDecayIntegral..." << endl;
		DecayIntegral = new TDecayIntegral();
		DecayIntegral->Init( mass+3, idOut+3, nop-2);
		Logger->Write( logINFO ) << "..Initialize TDecayIntegral DONE" << endl;
	}
	else 
	{
		DecayIntegral = NULL;
	}
};

////////////////////////////////////////////////////////////////////////
void TEventMaker2toN::SetEvent(int nDim,double *Xarg)
{
	//GENERATE MOMENTA AND RAPIDITIES starting from nDim random numbers
	// nDim=3*(number of particles)-4

	//nomber of true degrees of freedom
	assert( nDim == 3 * nop - 4 );
	

	//set up  new generation
	generationFailed = false;

    double p1t = abs(p_max-p_min)*Xarg[0];
    double p2t = abs(p_max-p_min)*Xarg[1];
    double dphi1 = 2.0*PI*Xarg[2];
    double dphi2 = 2.0*PI*Xarg[3];
    double y = y_min+(y_max-y_min)*Xarg[4];
   //intermediate mass
    double mass2 = mass2_min + (mass2_max-mass2_min)*Xarg[5];
	
	////clear queue - in case there was error in previous run
	while (!rndQueue.empty())
	{
		rndQueue.pop();
	}

	//put rnd numbers into queue
	for( int i = 6; i < nDim; i++)
	{
		rndQueue.push( Xarg[i] );
	}  
 
    int isol=1;
    int r=0; 
     
    r = CalculateKinematics(p1t, p2t, dphi1, dphi2, y, mass2, isol);
    if(!r==0)generationFailed=true;
    
    return;


};


////////////////////////////////////////////////////////////////////////
int TEventMaker2toN::CalculateKinematics(double p1t, double p2t, double dphi1, double dphi2, double y, double vmass2, int isol)
{
  TLorentzVector p1,p2,p3;

    //   four-momenta of the initial nucleons;
    //   ------------------------------------;
    double s=tecm*tecm;

    TVector2 pp1;
    pp1.SetMagPhi(p1t,dphi1);
    TVector2 pp2;
    pp2.SetMagPhi(p2t,dphi2);
    TVector2 pp3t = (pp1 + pp2);
    TVector2 pp3=pp3t.Rotate(PI);
    //---------------------------------
    // cental object (variable mass)

    double am3t = sqrt( vmass2 + pp3.Mod2() );

    double p3_0 = am3t*cosh(y);
    double p3_z = am3t*sinh(y);
    if(p3_0>tecm) return 1;
    //set four momentum of intermediate state (or central particle)
    p3.SetXYZT(pp3.X(), pp3.Y(), p3_z, p3_0);


   // solving set of equations for p_{1z} and p_{2z};
   // a p1z^2 + b p1z + c = 0

    double Em = p3.E();
    double pmz = p3.Z();

	double m1 = mass[1];
	double m2 = mass[2];

    //Coefficients of the equation for p1z = x   in  ax^2+bx+c = 0:
    double a =-4*(pow(Em,2) - pow(pmz,2) - 2*Em*tecm + pow(tecm,2));

    double b = -4*pmz*(pow(Em,2) + pow(m1,2) - pow(m2,2) + pow(p1t,2) - pow(p2t,2) - pow(pmz,2) - 2*Em*tecm + pow(tecm,2));

    double c = pow(Em,4) + pow(m1,4) + pow(m2,4) - 2*pow(m2,2)*pow(p1t,2) + pow(p1t,4) + 2*pow(m2,2)*pow(p2t,2) - 2*pow(p1t,2)*pow(p2t,2) + pow(p2t,4) + 2*pow(m2,2)*pow(pmz,2) - 2*pow(p1t,2)*pow(pmz,2) + 2*pow(p2t,2)*pow(pmz,2) + pow(pmz,4) - 4*pow(Em,3)*tecm - 2*pow(m2,2)*pow(tecm,2) - 2*pow(p1t,2)*pow(tecm,2) - 2*pow(p2t,2)*pow(tecm,2) - 
		2*pow(pmz,2)*pow(tecm,2) + pow(tecm,4) - 2*pow(Em,2)*(pow(m1,2) + pow(m2,2) + pow(p1t,2) + pow(p2t,2) + pow(pmz,2) - 3*pow(tecm,2)) + 4*Em*tecm*(pow(m1,2) + pow(m2,2) + pow(p1t,2) + pow(p2t,2) + pow(pmz,2) - pow(tecm,2)) - 2*pow(m1,2)*(pow(m2,2) - pow(p1t,2) + pow(p2t,2) + pow(pmz,2) + pow(tecm,2));

    //Delta of the quadratic equation 
    double Delta = b*b-4.0*a*c;

	//fail of generation
    if( Delta < 0)
    {
    	return 2;
    };

    //two solutions
    double p1za = (-b + sqrt(Delta))/(2.0*a);
    double p1zb = (-b - sqrt(Delta))/(2.0*a);

    double pp1z[2],pp2z[2];
    double pp10[2],pp20[2];

    //construct full solution
    pp1z[0] = p1za;
    pp2z[0] = -(p3.Z()+p1za);

    pp1z[1] = p1zb;
    pp2z[1] = -(p3.Z()+p1zb);

    //transverse mass
    double am1t2 = sq(mass[1]) + sq(p1t);
    double am2t2 = sq(mass[2]) + sq(p2t);

    //construct energies
    pp10[0] = sqrt(am1t2+sq(pp1z[0]));
    pp10[1] = sqrt(am1t2+sq(pp1z[1]));

    pp20[0] = sqrt(am2t2+sq(pp2z[0]));
    pp20[1] = sqrt(am2t2+sq(pp2z[1]));

    assert( isol <2 && isol >=0 );

    double p1testz = pp1z[isol];		
    double p2testz = pp2z[isol];		

    double energy1 = sqrt(am1t2+sq(p1testz));
    double energy2 = sqrt(am2t2+sq(p2testz));

    double esum = energy1+energy2+p3_0;

    if(fabs(esum-tecm)>0.001)
		return 3;
    
    p1.SetXYZT(pp1.X(), pp1.Y(), pp1z[isol], pp10[isol]); 
    p2.SetXYZT(pp2.X(), pp2.Y(), pp2z[isol], pp20[isol]);


    TLorentzVector p13,p23;
    p13 = p1+p3;
    p23 = p2+p3;

    if( p13.M2() > s ) 
		return 4;
    
    if( p23.M2() > s ) 
		return 6;

    if( pb1.Z()*p1.Z() < 0||pb2.Z()*p2.Z()< 0 )return 7; 
   TLorentzVector p0a=pb1 - p1 - p3;
    TLorentzVector p0b=pb2 - p2 - p3;
    
    
    //decay third particle and boost  
  
  //internal switch between TDecay decay into N >=4 partilces (1)
  //and decay into 4 particles (for testing)
  int decayOption = 1;
  
  if( decayOption == 1 )
  {
    
    Decay->SetDecay(p3, nop-2, mass+3);
	
	wtdecay = Decay->Generate( rndQueue );
	
	pf[1]=p1;
	pf[2]=p2;
    
    for( int i = 3; i < nop+1; i++ )
    {
		pf[i] = *(Decay->GetDecay( i-3 ));
	}
    
    PM = p3;
   
 }
 else if( decayOption == 2 )
 {
	double dcos=-1.0+2.0*rndQueue.front(); rndQueue.pop();
    double decphi=2.0*PI*rndQueue.front(); rndQueue.pop();
	 
	double dmom=kcms(vmass2,mass[3],mass[4]);// cm momenum for variable intermediate mass
    double dsin = sqrt(1.0-dcos*dcos);
    double dpx=dmom*dsin*cos(decphi);
    double dpy=dmom*dsin*sin(decphi);
    double dpz=dmom*dcos;
    double edec1=sqrt(dmom*dmom+mass[3]*mass[3]);
    double edec2=sqrt(dmom*dmom+mass[4]*mass[4]);
    TLorentzVector pdec1,pdec2,p4,p5;
    //set four-momenta of particles from decay of intermediate state in their cms
    pdec1.SetXYZT(dpx,dpy,dpz,edec1);
    pdec2.SetXYZT(-dpx,-dpy,-dpz,edec2);
    //beta of boost to overall ref. frame
    TVector3 bst=p3.BoostVector();
    //boost to overall ref. frame
    pdec1.Boost(bst);
    pdec2.Boost(bst);
    pf[3]=pdec1;
    pf[4]=pdec2;
   
	TLorentzVector pf5 = pf[3]+pf[4];
	wtdecay = kcms(pf5.M2(), mass[3], mass[4])/pf5.M()*M_PI;
		
    
	 
 }  
   
   return 0;

};

////////////////////////////////////////////////////////////////////////
double TEventMaker2toN::kcms(double s, double m1, double m2)
{
  // momentum in cms of systme of particles with masses m1, m2
  double m12=m1*m1;
  double m22=m2*m2;
  double temp= sqrt(pow((s-m12-m22),2)-4.0*m12*m22)/2.0/sqrt(s);
  return temp;

};

////////////////////////////////////////////////////////////////////////
double TEventMaker2toN::SetEvent( int nDim, double *Xarg, TEvent * event )
{
	SetEvent( nDim, Xarg );

	//Set particles in corresponding frame
    event->pb[1] = pb1;
    event->pbInfo[1] = PDGDatabese->GetParticle( idA );
    event->pb[2] = pb2;
    event->pbInfo[2] = PDGDatabese->GetParticle( idB );
    for( int i = 1; i < nop+1; i++)
    {
		event->pf[i] = pf[i];
		event->pfInfo[i] = PDGDatabese->GetParticle( idOut[i] );
    }

	
	
	//generate weight 
	double norm3=pow( 2.0 * PI, -3*nop + 4 )*pow(2,-3);
  
	double wt = 0.0;
	
	double jacob_rndm_3 = 0.0;
	
	if ( false == generationFailed )
	{
         double jacob_rndm_3 = abs(p_max-p_min) * abs(p_max-p_min) * 2 *PI* 2*PI*abs(y_max-y_min)*(mass2_max-mass2_min);//PhaseSpace-->dXarg for a+b->1+2+M
         
         double jacobinv_phs = fabs( pf[1].Z()/pf[1].E() - pf[2].Z()/pf[2].E());
         
         double jacob3 = jacob_rndm_3/jacobinv_phs;

		 double flux = 2.0*tecm*tecm;
		 double gev2tomb = 0.3894;
		 // wt = pow(gev2tomb,-3)/flux;
		 wt = 1.0/flux;
				
	     double wt3 = wt * norm3 * jacob3 * pf[1].Pt() * pf[2].Pt()/pf[1].E()/pf[2].E();
	     
		 wt = wt3;
		 
		 wt *= wtdecay;
	 
	}
	else //generation failed
	{
		wt   = 0.0;
	}
	
	weight = wt/jacob_rndm_3;
	
	return( wt );
};

////////////////////////////////////////////////////////////////////////
bool TEventMaker2toN::isGenerationFailed( void )
{
	
	return( generationFailed );
};

////////////////////////////////////////////////////////////////////////
double TEventMaker2toN::getWeight( void )
{
	return( weight );
};
	
////////////////////////////////////////////////////////////////////////
double TEventMaker2toN::getDecayWeight( void )
{
	return( wtdecay );
};
	
////////////////////////////////////////////////////////////////////////
double TEventMaker2toN::getDecayPhaseSpaceIntegral( double M )
{
	if( bIDecay == 0 )
	{
		Logger->Write( logERROR ) << "TEventMaker2toN:: Set TEventMaker2toN::bIDecay to 1 in order to use getDecayPhaseSpaceIntegral( M ) - STOP" << endl;
		cerr << "TEventMaker2toN:: Set TEventMaker2toN::bIDecay to 1 in order to use getDecayPhaseSpaceIntegral( M ) - STOP" << endl;
		throw string( "TEventMaker2toN:: Set TEventMaker2toN::bIDecay to 1 in order to use getDecayPhaseSpaceIntegral( M ) - STOP" );
	}
	
	return( DecayIntegral->getDecayIntegral( M ) );
};	
