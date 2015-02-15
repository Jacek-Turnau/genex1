/***********************************************************************
*                       TEventMakerCylindricalPhaseSpace                                 
*                                                                           
* Generate particles and phase space weight according to 
* Cylindical phase space generator given by S. Jadach Comp.Phys.Commun 9(175)297         
*   									       
***********************************************************************/


#include"TEventMakerCylindricalPhaseSpace.h"

////////////////////////////////////////////////////////////////////////

const double TEventMakerCylindricalPhaseSpace::PI =  M_PI;

////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////
TEventMakerCylindricalPhaseSpace::TEventMakerCylindricalPhaseSpace()
{
  
	ConfigReader = ConfigReaderSubsystem::Instance();
	
	Logger = LoggerSubsystem::Instance();

	PDGDatabese = TDatabasePDG::Instance();
  
	ReadConfigFile( string("Generator.dat") );
  
	mA = PDGDatabese->GetParticle(idA)->Mass();
	mB = PDGDatabese->GetParticle(idB)->Mass(); 
  
	///setup CM frame
	EA = (pow(mA,2) - pow(mB,2) + pow(tecm,2))/(2.0 * tecm);
	EB = (-pow(mA,2) + pow(mB,2) + pow(tecm,2))/(2.0 *tecm);
		
	double pzA = sqrt( EA*EA - mA*mA );
	double pzB = -sqrt( EB*EB - mB*mB );
		
	pb1.SetPxPyPzE( 0.0, 0.0, pzA, EA );
	pb2.SetPxPyPzE( 0.0, 0.0, pzB, EB );
  
  
	//allocate tables:
  
	for( int i = 0; i < 6; i++ )
		pcm [i] = new double [nop+1];
  
	//mass  = new double [nop+1];
	tmass = new double [nop+1];
  
  
	//initialize internal structure of the class
	Initialize();
  
};

////////////////////////////////////////////////////////////////////////
TEventMakerCylindricalPhaseSpace::~TEventMakerCylindricalPhaseSpace(void)
{

   //delete pcm matrix
   for( int i = 0; i < 6; i++ )
     delete [] (pcm [i]);
   
   //delete [] mass;
   delete [] tmass;
   
   delete [] mass;
   delete [] idOut;
   
   delete [] pf;
  
};

////////////////////////////////////////////////////////////////////////
void TEventMakerCylindricalPhaseSpace::ReadConfigFile( const string & filename )
{	
	
	//Parse variables: setting up model variables	

	ptt = ConfigReader->GetDoubleValue( "TEventMakerCylindricalPhaseSpace::ptt" );
	nop = ConfigReader->GetIntValue( "TEventMakerCylindricalPhaseSpace::nop" );
	
	
	//Reading frame setup
	frameType = ConfigReader->GetIntValue( string("frameType") );
	idA = ConfigReader->GetIntValue( string("idA") );
	idB = ConfigReader->GetIntValue( string("idB") );
	tecm = ConfigReader->GetDoubleValue( string("tecm") );
 

	//allocate tables
	mass  = new double [nop+1];
	idOut = new int[nop+1];
	pf = new TLorentzVector [nop+1];


	//read output particles 
	int massIndication = ConfigReader->GetIntValue( string("TEventMakerCylindricalPhaseSpace::massIndication") );

	if( massIndication == 1)
	{
		ConfigReader->GetDoubleArray( mass, string("TEventMakerCylindricalPhaseSpace::m"), nop+1, 1 );
		
		for( int i = 1; i < nop+1; i++ )
		{
			idOut[i] = 81; //unkonown particle if mass specified
			//cout << "id= " << idOut[i] << ", m = " << mass[i] << endl;
		}
		
	}
	else if( massIndication == 0 )
	{
		ConfigReader->GetIntArray( idOut, string("TEventMakerCylindricalPhaseSpace::id"), nop+1, 1 );	
		
		for( int i = 1; i < nop+1; i++ )
		{
			mass[i] = PDGDatabese->GetParticle(idOut[i])->Mass();
			//cout << "id= " << idOut[i] << ", m = " << mass[i] << endl;
		}
	}
	else if( massIndication == 2)
	{
		ConfigReader->GetDoubleArray( mass, string("TEventMakerCylindricalPhaseSpace::m"), nop+1, 1 );
		ConfigReader->GetIntArray( idOut, string("TEventMakerCylindricalPhaseSpace::id"), nop+1, 1 );
	}
	

	//Log configuration
	string ConfigLog = ConfigReader->GetStringValue( string("ConfigLogFile") );
	TLog * ConfigLogger = new TLog();
	ConfigLogger->setLogTime( false );
	ConfigLogger->openLogFile( ConfigLog );
 
	ConfigLogger->Write() << "# " << "-----------" << endl;
	ConfigLogger->Write() << "# " << "TEventMakerCylindricalPhaseSpace " << endl;
	ConfigLogger->Write() << "# " << "-----------" << endl;
	ConfigLogger->Write() << "TEventMakerCylindricalPhaseSpace::ptt = " << ptt << endl;
	ConfigLogger->Write() << "TEventMakerCylindricalPhaseSpace::nop = " << nop << endl;
	ConfigLogger->Write() << "tecm = " << tecm << endl;
	ConfigLogger->Write() << "frameType = " << frameType << endl;
	ConfigLogger->Write() << "idA = " << idA << endl;
	ConfigLogger->Write() << "idB = " << idB << endl;
	ConfigLogger->Write() << "tecm = " << tecm << endl;
	
	
	ConfigLogger->Write() << "TEventMakerCylindricalPhaseSpace::massIndication = " << massIndication << endl;
	for( int i = 1; i < 5; i++ )
	{
			ConfigLogger->Write() <<  "TEventMakerCylindricalPhaseSpace::id" << i << " = " << idOut[i] << endl;
			ConfigLogger->Write() <<  "TEventMakerCylindricalPhaseSpace::m" << i << " = " << mass[i] << endl;
	}
	
	
	delete ConfigLogger;
	
	return;
};

////////////////////////////////////////////////////////////////////////
void TEventMakerCylindricalPhaseSpace::Initialize()
{ 
    R    = 0.0; 
    nNewtonbad = 0.0;
    generationFailed = false;
    sumMt = 0.0;
    sumPt2 = 0.0; 
    yNewt = 0.0;
    isAccepted = true;
    

	#ifdef DEBUG
	cout << "INITIAL CONDITIONS:" << endl;
	cout << "Total energy in CM frame = " << tecm << endl;

	cout << "Mass matrix of outgoing particles:" << endl;
		for( int i=1; i < nop+1; i++ )
			cout << "mass[" << i << "] = " << mass[i] << endl;
	#endif

	//dispersion of momenta
	double ptr = ptt;

	#ifdef CORRECT_MOMENTA_DISPERSION
	  //correct momenta dispersion for low energies
	  double massSum = 0.0;
	  for ( int i = 1; i < nop+1; i++)
		  massSum += mass[i];
	  double Q = 0.85 * sqrt( pow((tecm/long(nop)), 2) + pow( (massSum/long(nop)), 2) );
	  ptr = 1.0/sqrt( pow(1.0/ptt, 2) + pow(1.0/Q, 2) );
	#endif

	#ifdef DEBUG 
	cout << "Average momentum = " << ptt << endl;
	cout << "Corrected average momentum = " << ptr << endl;
	#endif
	// generate dispersion of transverse momenta
		R = (4.0/PI) * pow(ptr, 2) *( 1.0 + 1.0/long(nop-1) );  
		
	#ifdef DEBUG 
	  cout << "R = " << R << endl;
	#endif      
			
		
};

////////////////////////////////////////////////////////////////////////
void TEventMakerCylindricalPhaseSpace::SetEvent(int nDim,double *Xarg)
{
//GENERATE MOMENTA AND RAPIDITIES starting from nDim random numbers
// nDim=3*(number of particles)


	//set up  new generation
	generationFailed = false;

	//temp varaibles
	double a = 0.0;
	double b = 0.0;
	
	//sum of transverse momenta components
	double sumU = 0.0;
	double sumV = 0.0;

	// max and min rapidities
	double ymin = 1.0;
	int    imin = 1;
	double ymax = 0.0;
	int    imax = 1;

	//transversal momenta
	//double u [nop+1] = {0};
	//double v [nop+1] = {0};
	double * u = new double [nop+1];
	double * v = new double [nop+1];


	//rapidity
	//double rapidity [nop+1] = {0};
	double * rapidity = new double [nop+1];


	for ( int i = 1; i < nop+1; i++ )
	{
    // generate gaussian distribution from uniform one using BOX-MULLER algorithm
	  a = sqrt( - R * log( Xarg[3*(i-1)]) );
	  b = 2.0 * PI * Xarg[3*(i-1)+1];
		
		u[i] = a * sin(b);
		v[i] = a * cos(b);
	
		//additional summing in one loop - to be used later
		sumU += u[i];
		sumV += v[i];


		//rapidities
		rapidity[i] = Xarg[3*(i-1)+2];

		//searching for min and max rapidities in one loop - to be used later
		if ( rapidity[i]  > ymax )  
		{
			ymax = rapidity[i];
			imax = i;
		}		
		if ( rapidity[i] < ymin )
		{
			ymin = rapidity[i];
			imin = i;
		}

	}

	#ifdef DEBUGSetEvent
	  cout << "ymin = " << ymin << endl;
	  cout << "ymax = " << ymax << endl;
	#endif


	//scaling rapidities from (0;1) to [0;1] and calculating transverse momenta and masses

	sumMt = 0.0;	//sum of the transverse masses
	sumPt2 = 0.0;	//sum of the transverse momenta
	for ( int i = 1; i < nop+1; i++ )
	{
		rapidity[i] = ( rapidity[i] - ymin )/( ymax - ymin );

		pcm[2][i] = u[i] - sumU/long(nop);
		pcm[3][i] = v[i] - sumV/long(nop);

		tmass[i] = sqrt( pow( pcm[2][i], 2) + pow( pcm[3][i], 2 ) + pow( mass[i], 2 ) );
		
		sumMt += tmass[i];
		sumPt2 += pow( pcm[2][i], 2 ) +  pow( pcm[3][i], 2 );
	}

	//check if event is correctly generated
	if ( sumMt > tecm ) 
	{
		generationFailed = true;

		#ifdef DEBUG
		  cout << "BAD EVENT ( transverse mass  > total energi in CM) :(" << endl; 
		#endif
		
		delete [] u;
		delete [] v;
		delete [] rapidity;
		
		return; //generation failed
	}


	// iterative solution of equations for Y and Z (Newton iteration)
	const double eps = 0.0000001;	//precision of Newton method
	const int nMaxNewtonIter = 100;  //maximal number of Newton iterations

	//helper variables to calculating derivative
	double s1 = 0.0;
	double s2 = 0.0;
	double s3 = 0.0;
	double s4 = 0.0;

	yNewt	    = 0.0; //approximate Y 
	double dY   = 0.0; //change of Y from Newton method
	deriv       = 0.0; //derivative

	int j = 0 ;	//number of Newton iterations
	
	//derive initial value of Y for Newton method
	yNewt = 2.0 * log(tecm) - log( mass[imin] * mass[imax] );

	#ifdef DEBUGNewt
	  cout << "Start Newton iterations" << endl;
	  cout << "Initially Y = " << yNewt << endl;
	#endif

	do
	{
		s1 = 0.0;
		s2 = 0.0;
		s3 = 0.0;
		s4 = 0.0;		


		//calculate function and its derivative
		for ( int i = 1; i < nop+1; i++ )
		{
			s1 += tmass[i] * exp( rapidity[i] * yNewt );
			s2 += tmass[i] / exp( rapidity[i] * yNewt );
			s3 += tmass[i] * rapidity[i] * exp( rapidity[i] * yNewt );
			s4 += tmass[i] * rapidity[i] / exp( rapidity[i] * yNewt );
		}
		deriv = (s3/s1) - (s4/s2);

		dY = ( 2.0 * log(tecm) - log( s1 * s2 ) ) / deriv;
		yNewt = yNewt + dY; 

		j++;
		
		if ( j > nMaxNewtonIter ) 
		{
			generationFailed = true;
			nNewtonbad ++;
		
			delete [] u;
			delete [] v;
			delete [] rapidity;
		
			return; //generation failed
		}
	
	}
	while ( fabs( dY/yNewt )  > ( eps/long(nop) ) );

	#ifdef DEBUGNewt
	  cout << j << "Newton interations performed" << endl;
	  cout << "New Y = " << yNewt << endl;
	#endif



	//calculate longitudinal momenta, energies and rapidities 
	double y = 0.0;
	double z = 0.0; 
	s1 = 0.0;
	
	for ( int i = 1; i < nop+1; i++ )
	{
		s1 += tmass[i] * exp( rapidity[i] * yNewt );
	}

	z = log( tecm/s1 );	

	for ( int i = 1; i < nop+1; i++ )
	{
		y =  rapidity[i] * yNewt + z;
		
		pcm[1][i] = 0.5 * tmass[i] * ( exp( y ) - exp( -y ) );
		pcm[4][i] = 0.5 * tmass[i] * ( exp( y ) + exp( -y ) );
		pcm[5][i] = y;

	}
	
	delete [] u;
	delete [] v;
	delete [] rapidity;
		
	
	return;

  
};

////////////////////////////////////////////////////////////////////////
double TEventMakerCylindricalPhaseSpace::SetEvent( int nDim, double *Xarg, TEvent * event )
{
	SetEvent( nDim, Xarg );
	
	//Set particles in corresponding frame
    event->pb[1] = pb1;
    event->pbInfo[1] = PDGDatabese->GetParticle( idA );
    event->pb[2] = pb2;
    event->pbInfo[2] = PDGDatabese->GetParticle( idB );
     
    for(int i=1;i<nop+1;i++)
    {
		event->pf[i].SetPxPyPzE(pcm[2][i],pcm[3][i],pcm[1][i],pcm[4][i]);
		event->pfInfo[i] = PDGDatabese->GetParticle( idOut[i] );
    }
	
	// fix raction geometry 
	// if there is backscattering weight is 0
	if( ( event->pf[1].Pz() <= 0 ) && ( event->pf[2].Pz() >= 0 ) )
	{
		return( 0.0 );
	}	
	//else...
	
	
	
	//generate weight of Cylindrical Phase Space
	double flux = 2.0*tecm*tecm;
    double gev2tomb = 0.3894;
	double norm = pow( 2.0 * PI, -3*nop+4 );
	
	double wt = norm*pow(gev2tomb,-3)/flux;
	
	
	if ( false == generationFailed )
	{
		//transverse part
		wt = wt * pow( (PI*R), nop-1 ) * exp( sumPt2/R)/long(nop);
		//longitudinal part
		wt = wt * long(nop * (nop -1) ) * pow( yNewt, nop - 2 )* pow( 0.5, nop - 1)/(tecm * tecm * deriv);
	}
	else //generation failed
	{
		wt   = 0.0;
	}
	
	return( wt );
	
};

////////////////////////////////////////////////////////////////////////
bool TEventMakerCylindricalPhaseSpace::isGenerationFailed( void )
{
	return( generationFailed );
};
