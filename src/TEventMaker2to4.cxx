/***********************************************************************
*                         TEventMaker2to4                                 
* 
* Generates phase space for a+b -> a+b + C + D  event.              
*   									       
************************************************************************/
#include"TEventMaker2to4.h"

////////////////////////////////////////////////////////////////////////

const double TEventMaker2to4::PI =  M_PI;

////////////////////////////////////////////////////////////////////////
TEventMaker2to4::TEventMaker2to4()
{
  
	ConfigReader = ConfigReaderSubsystem::Instance();
	
	Logger = LoggerSubsystem::Instance();
  
	ReadConfigFile( string("Generator.dat") );
  
	PDGDatabese = TDatabasePDG::Instance();
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
TEventMaker2to4::~TEventMaker2to4(void)
{

	delete [] mass;
	delete [] idOut;
	delete [] pf;
   
};

////////////////////////////////////////////////////////////////////////
void TEventMaker2to4::ReadConfigFile( const string & filename )
{	
	
	//Parse variables: setting up model variables	

	nop = ConfigReader->GetIntValue( "nop" );
	p_min = ConfigReader->GetDoubleValue( "TEventMaker2to4::p_min" );
	p_max = ConfigReader->GetDoubleValue( "TEventMaker2to4::p_max" );
	y_min = ConfigReader->GetDoubleValue( "TEventMaker2to4::y_min" );
	y_max = ConfigReader->GetDoubleValue( "TEventMaker2to4::y_max" );
	pmt_max = ConfigReader->GetDoubleValue( "TEventMaker2to4::pmt_max" );
	weightStrategy = ConfigReader->GetIntValue( "weightStrategy" );

	///Reading frame setup
	frameType = ConfigReader->GetIntValue( string("frameType") );
	idA = ConfigReader->GetIntValue( string("idA") );
	idB = ConfigReader->GetIntValue( string("idB") );
	tecm = ConfigReader->GetDoubleValue( string("tecm") );

	assert( nop == 4 );

	//allocate tables
	mass  = new double [nop+1];
	idOut = new int[nop+1];
	pf = new TLorentzVector [nop+1];
	

	//read output particles 
	int massIndication = ConfigReader->GetIntValue( string("TEventMaker2to4::massIndication") );
	
	if( massIndication == 1)
	{
		ConfigReader->GetDoubleArray( mass, string("TEventMaker2to4::m"), nop+1, 1 );
		
		for( int i = 1; i < nop+1; i++ )
		{
			idOut[i] = 81; //unkonown particle if mass specified
			//cout << "id= " << idOut[i] << ", m = " << mass[i] << endl;
		}
		
	}
	else if( massIndication == 0 )
	{
		ConfigReader->GetIntArray( idOut, string("TEventMaker2to4::id"), nop+1, 1 );	
		
		for( int i = 1; i < nop+1; i++ )
		{
			mass[i] = PDGDatabese->GetParticle(idOut[i])->Mass();
			//cout << "id= " << idOut[i] << ", m = " << mass[i] << endl;
		}
	}
	else if( massIndication == 2)
	{
		ConfigReader->GetDoubleArray( mass, string("TEventMaker2to4::m"), nop+1, 1 );
		ConfigReader->GetIntArray( idOut, string("TEventMaker2to4::id"), nop+1, 1 );
	}

	
	//check if leading particles are the same
	assert( idOut[1] == idOut[2] );
	
	
	//Log configuration
	string ConfigLog = ConfigReader->GetStringValue( string("ConfigLogFile") );
	TLog * ConfigLogger = new TLog();
	ConfigLogger->setLogTime( false );
	ConfigLogger->openLogFile( ConfigLog );
 
	ConfigLogger->Write() << "# " << "-----------" << endl;
	ConfigLogger->Write() << "# " << "TEventMaker2to4 " << endl;
	ConfigLogger->Write() << "# " << "-----------" << endl;
	ConfigLogger->Write() << "nop = " << nop << endl;
	ConfigLogger->Write() << "TEventMaker2to4::p_min = " << p_min << endl;
	ConfigLogger->Write() << "TEventMaker2to4::p_max = " << p_max << endl;
	ConfigLogger->Write() << "TEventMaker2to4::y_min = " << y_min << endl;
	ConfigLogger->Write() << "TEventMaker2to4::y_max = " << y_max << endl;
	ConfigLogger->Write() << "TEventMaker2to4::pmt_max = " << pmt_max << endl;
	ConfigLogger->Write() << "weightStrategy = " << weightStrategy << endl;
	
	ConfigLogger->Write() << "frameType = " << frameType << endl;
	ConfigLogger->Write() << "idA = " << idA << endl;
	ConfigLogger->Write() << "idB = " << idB << endl;
	ConfigLogger->Write() << "tecm = " << tecm << endl;
	
	
	ConfigLogger->Write() << "TEventMaker2to4::massIndication = " << massIndication << endl;
	for( int i = 1; i < 5; i++ )
	{
			ConfigLogger->Write() <<  "TEventMaker2to4::id" << i << " = " << idOut[i] << endl;
			ConfigLogger->Write() <<  "TEventMaker2to4::m" << i << " = " << mass[i] << endl;
	}
	
	
	delete ConfigLogger;

	
	return;
};

////////////////////////////////////////////////////////////////////////
void TEventMaker2to4::Initialize()
{ 
		
};

////////////////////////////////////////////////////////////////////////
void TEventMaker2to4::SetEvent(int nDim,double *Xarg)
{
//GENERATE MOMENTA AND RAPIDITIES starting from nDim random numbers
// nDim=3*(number of particles)-4


	//set up  new generation
	generationFailed = false;

	double p1t = abs(p_max-p_min)*Xarg[0];
	double p2t = abs(p_max-p_min)*Xarg[1];
    double phi1 = 2.0*PI*Xarg[2];
    double phi2 = 2.0*PI*Xarg[3];
    double y3 = y_min+(y_max-y_min)*Xarg[4];
    double y4 = y_min+(y_max-y_min)*Xarg[5];
    double pmt = pmt_max*Xarg[6];
    double phipmt = 2.0*PI*Xarg[7];
    int isol;
    int r=0;

    //pick one of the roots
	isol = 1;

    r = CalculateKinematics(p1t,p2t,phi1,phi2,y3,y4,pmt,phipmt,isol);
    if(!r==0)
    	{
    	generationFailed=true;
    	//cout << "FAILED r=" << r << endl;
    	}
    return;

  
};

////////////////////////////////////////////////////////////////////////
int TEventMaker2to4::CalculateKinematics(double p1t, double p2t, double phi1, double phi2, double y3, double y4, double pmt, double phipmt, int isol)
{
	TLorentzVector p1,p2,p3,p4;

	double mout1 = mass[3];
	double mout2 = mass[4];


	//   four-momenta of the initial nucleons;
	//   ------------------------------------;
	double s=tecm*tecm;
   
   
    TVector2 pp1;
    pp1.SetMagPhi(p1t,phi1);
    TVector2 pp2;
    pp2.SetMagPhi(p2t,phi2);
    TVector2 pm;
    pm.SetMagPhi(pmt,phipmt);

    //four-momenta of all involved particles;

    TVector2 pp3 = (pm - pp1 - pp2)/2.;
    TVector2 pp4 = ((-1)*pm - pp1 - pp2)/2.;

    double am3t = sqrt( sq(mout1) + pp3.Mod2() );
    double am4t = sqrt( sq(mout2) + pp4.Mod2() );

    double p3_0 = am3t*cosh(y3);
    double p3_z = am3t*sinh(y3);

    double p4_0 = am4t*cosh(y4);
    double p4_z = am4t*sinh(y4);

    p3.SetXYZT(pp3.X(), pp3.Y(), p3_z, p3_0);
    p4.SetXYZT(pp4.X(), pp4.Y(), p4_z, p4_0);

    if((p3_0+p4_0)>tecm) return 1;


   // solving set of equations for p_{1z} and p_{2z};
   // a p1z^2 + b p1z + c = 0

    double Em = p3_0 + p4_0;
    double pmz = p3_z+p4_z;

    //Coefficients of the equation for p1z = x   in  ax^2+bx+c = 0:
    double a =-4*(pow(Em,2) - pow(pmz,2) - 2*Em*tecm + pow(tecm,2));

    double b = 4*pmz*(-pow(Em,2) - pow(p1t,2) + pow(p2t,2) + pow(pmz,2) + 2*Em*tecm - pow(tecm,2));

    double c = pow(Em,4) + pow(p1t,4) + pow(p2t,4) + 2*pow(p2t,2)*pow(pmz,2) + pow(pmz,4) -
    		   4*pow(Em,3)*tecm - 4*pow(mass[1],2)*pow(tecm,2) - 2*pow(p2t,2)*pow(tecm,2) -
    		   2*pow(pmz,2)*pow(tecm,2) + pow(tecm,4) -
    		   2*pow(Em,2)*(2*pow(mass[1],2) + pow(p1t,2) + pow(p2t,2) + pow(pmz,2) -
    		      3*pow(tecm,2)) + 4*Em*tecm*(2*pow(mass[1],2) + pow(p1t,2) + pow(p2t,2) +
    		      pow(pmz,2) - pow(tecm,2)) -
    		   2*pow(p1t,2)*(pow(p2t,2) + pow(pmz,2) + pow(tecm,2));

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
    pp2z[0] = -(p3_z+p4_z+p1za);

    pp1z[1] = p1zb;
    pp2z[1] = -(p3_z+p4_z+p1zb);


    //transverse mass
    double am1t2 = sq(mass[1]) + sq(p1t);
    double am2t2 = sq(mass[1]) + sq(p2t);

    //construct energies
    pp10[0] = sqrt(am1t2+sq(pp1z[0]));
    pp10[1] = sqrt(am1t2+sq(pp1z[1]));

    pp20[0] = sqrt(am2t2+sq(pp2z[0]));
    pp20[1] = sqrt(am2t2+sq(pp2z[1]));

    assert( isol <2 && isol >=0 );
    

    //pick one of the solutions
    double p1testz = pp1z[isol];		
    double p2testz = pp2z[isol];		

    double energy1 = sqrt(am1t2+sq(p1testz));
    double energy2 = sqrt(am2t2+sq(p2testz));

    double esum = energy1+energy2+p3_0+p4_0;

    if(fabs(esum-tecm)>0.001)return 3;

    //set 4-momentum of protons
    TLorentzVector p13,p23,p14,p24;
    p1.SetXYZT(pp1.X(), pp1.Y(), pp1z[isol], pp10[isol]); 
    p2.SetXYZT(pp2.X(), pp2.Y(), pp2z[isol], pp20[isol]);

    p13 = p1+p3;
    p14 = p1+p4;
    p23 = p2+p3;
    p24 = p2+p4;

    if( p13.M2() > s ) return 4;
    if( p14.M2() > s ) return 5;
    if( p23.M2() > s ) return 6;
    if( p24.M2() > s ) return 7;

    TLorentzVector p0a_1 = pb1 - p1 - p3;
    TLorentzVector p0a_2 = pb2 - p2 - p4;
    TLorentzVector p0b_1 = pb1 - p1 - p4;
    TLorentzVector p0b_2 = pb2 - p2 - p3;


    //set fourvectors - setting reflection of particles or going along initial directions

	pf[1]=p1;
	pf[2]=p2;

    pf[3]=p3;
    pf[4]=p4;

   return 0;
};


////////////////////////////////////////////////////////////////////////
double TEventMaker2to4::SetEvent( int nDim, double *Xarg, TEvent * event )
{
	
	SetEvent( nDim, Xarg );

	
	///Set particles in CM frame
    event->pb[1] = pb1;
    event->pbInfo[1] = PDGDatabese->GetParticle( idA );
    event->pb[2] = pb2;
    event->pbInfo[2] = PDGDatabese->GetParticle( idB );

	for( int i = 1; i < nop+1; i++)
    {
		event->pf[i] = pf[i];
		event->pfInfo[i] = PDGDatabese->GetParticle( idOut[i] );
    }
    
	
	//generate weight of event
	
	double flux = 2.0*tecm*tecm;
    double gev2tomb = 0.3894;
	double norm = pow( 2.0 * PI, -8 );
	
	double wt = norm*pow(gev2tomb,-3)/flux;
	
	
	if ( false == generationFailed )
	{


		//double jacobian = p_max*p_max * pmt_max * 8*PI*PI*PI * (y_max-y_min)*(y_max-y_min);//PhaseSpace-->dXarg
		double jacobian = abs(p_max-p_min)*abs(p_max-p_min) * pmt_max * 8*PI*PI*PI * (y_max-y_min)*(y_max-y_min);//PhaseSpace-->dXarg
		double ajacobinv = fabs( pf[1].Z()/pf[1].E() - pf[2].Z()/pf[2].E());
		double ajacob = 1./ajacobinv;
		wt = wt* ajacob*jacobian*pf[1].Pt()*pf[2].Pt()*((pf[3]-pf[4]).Pt())/pf[1].E()/pf[2].E();

		wt = wt * pow(2,-6);
		
	}
	else //generation failed
	{
		
		wt = 0.0;
		
		//cout << "Generation FAILED" << endl;
		
	}
	

	return( wt );
	
};


////////////////////////////////////////////////////////////////////////
bool TEventMaker2to4::isGenerationFailed( void )
{
		return( generationFailed );
};
