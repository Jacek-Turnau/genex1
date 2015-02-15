/***********************************************************************
		File defines cuts for exploration phase (if used )

***********************************************************************/ 	
	
	/*
		Fill below section using specific cuts.

		All 4-momenta are aviable
		
		Description of aviable variables:
		TLorenzVectors  pb1, pb2 - 4-momenta of the beam particles
		TLorenzVectors pf[1]...pf[nop] - 4-momenta of the outgoing particles

		double tecm - CM energy of reaction
	
		Example:
		
		//create useful kinematics variables
		double t1=(pb1-pf[1]).M2();
		double t2=(pb2-pf[2]).M2();
		//cuts
		if( ( abs( t1 ) < 1.0 ) && ( abs( t2 ) < 1.0 ) )
			Acc = 1.0;  //case passes acceptance test
		else
			Acc = 0.0; //otherwise do not pass
	
	*/


////////////////////////////////////////////////////////////////////////
double TExplorationAcceptance::cut0( void )
{
		return 1.0;
};

////////////////////////////////////////////////////////////////////////
double TExplorationAcceptance::cut1( void )
{
	double Acc = 0.0;
	
	double ptpip = pf[3].Pt();
	double ptpim = pf[4].Pt();
	
	double etapip = pf[3].PseudoRapidity();
	double etapim = pf[4].PseudoRapidity();

	
	Acc = 0.0;
	
	if( ( ptpip > 0.15 )&&( ptpim > 0.15 )&&( abs(etapim) < 1.0 )&&( abs(etapip) < 1.0 ))
	{
		Acc = 1.0;
	}
	else
	{
		Acc = 0.0;
	}
		
	return( Acc );	
		
};

////////////////////////////////////////////////////////////////////////
double TExplorationAcceptance::cut2( void )
{
	double Acc = 0.0;
	
	double t1=(pb1-pf[1]).M2();
    double t2=(pb2-pf[2]).M2();
	
	double ptpip = pf[3].Pt();
	double ptpim = pf[4].Pt();
	
	double etapip = pf[3].PseudoRapidity();
	double etapim = pf[4].PseudoRapidity();

	
	Acc = 0.0;
	
	if( ( ptpip > 0.15 )&&( ptpim > 0.15 )&&( abs(etapim) < 1.0 )&&( abs(etapip) < 1.0 )&&( abs(t1) > 0.003 )&&( abs(t1) < 0.035 )&&( abs(t2) > 0.003 )&&( abs(t2) < 0.035 ) )	//200GeV
	{
		Acc = 1.0;
	}
	else
	{
		Acc = 0.0;
	}
		
	return( Acc );		
	
};

////////////////////////////////////////////////////////////////////////
double TExplorationAcceptance::cut3( void )
{
	double Acc = 0.0;
	
	double t1=(pb1-pf[1]).M2();
    double t2=(pb2-pf[2]).M2();
	
	double ptpip = pf[3].Pt();
	double ptpim = pf[4].Pt();
	
	double etapip = pf[3].PseudoRapidity();
	double etapim = pf[4].PseudoRapidity();

	
	Acc = 0.0;
	
	
	if( ( ptpip > 0.15 )&&( ptpim > 0.15 )&&( abs(etapim) < 1.0 )&&( abs(etapip) < 1.0 )&&( abs(t1) > 0.1 )&&( abs(t1) < 1.5 )&&( abs(t2) > 0.1 )&&( abs(t2) < 1.5 ) ) // 500GeV
	{
		Acc = 1.0;
	}
	else
	{
		Acc = 0.0;
	}
		
	return( Acc );


};

////////////////////////////////////////////////////////////////////////
double TExplorationAcceptance::cut4( void )
{
	double Acc = 1.0;
	
	double ptpip = pf[3].Pt();
	double ptpim = pf[4].Pt();
	
	double etapip = pf[3].PseudoRapidity();
	double etapim = pf[4].PseudoRapidity();
	
	
	
	Acc = 0.0;
	
	if( ( ptpip > 0.1 )&&( ptpim > 0.1 )&&( abs(etapim) < 0.9 )&&( abs(etapip) < 0.9 ) )    //1st cut
	{
		Acc = 1.0;
	}
	else
	{
		Acc = 0.0;
	}
	

	return( Acc );

};

////////////////////////////////////////////////////////////////////////
double TExplorationAcceptance::cut5( void )
{
	double Acc = 1.0;
	
	double ptpip = pf[3].Pt();
	double ptpim = pf[4].Pt();
	
	double etapip = pf[3].PseudoRapidity();
	double etapim = pf[4].PseudoRapidity();
	
	
	
	Acc = 0.0;
	
	if( ( ptpip > 0.3 )&&( ptpim > 0.3 )&&( abs(etapim) < 0.9 )&&( abs(etapip) < 0.9 ) )	//2nd cut
	{
		Acc = 1.0;
	}
	else
	{
		Acc = 0.0;
	}
	

	return( Acc );

		
};

////////////////////////////////////////////////////////////////////////
double TExplorationAcceptance::cut6( void )
{
	double Acc = 0.0;
	
	double y3 = pf[3].Rapidity();
	double y4 = pf[4].Rapidity();
	double y5 = pf[5].Rapidity();
	
	if( (abs(y3) <2.0) && (abs(y4) <2.0) && (abs(y5) <2.0)  )
		Acc=1.0;
		
	return( Acc );	
	
};

////////////////////////////////////////////////////////////////////////
double TExplorationAcceptance::cut7( void )
{
	double Acc = 0.0;
	
	double t1=(pb1-pf[1]).M2();
    double t2=(pb2-pf[2]).M2();	

		
	const double tUPcut = 0.06;   //GeV
	const double tDOWNcut = 0.01; //GeV
	
	double y3 = pf[3].Rapidity();
	double y4 = pf[4].Rapidity();
	double y5 = pf[5].Rapidity();

	
	if( (abs(y3) <2.0) && (abs(y4) <2.0) && (abs(y5) <2.0)  && ( abs(t1) > tDOWNcut) &&  (abs(t1) < tUPcut) && ( abs(t2) > tDOWNcut) &&  (abs(t2) < tUPcut) )
		Acc=1.0;
		
	return( Acc );	
};

////////////////////////////////////////////////////////////////////////
double TExplorationAcceptance::cut8( void )
{
	double Acc = 0.0;
	
	double y3 = pf[3].Rapidity();
	double y4 = pf[4].Rapidity();
	double y5 = pf[5].Rapidity();
	
	if( (abs(y3) <2.0) && (abs(y4) <2.0) && (abs(y5) <2.0)  )
		Acc=1.0;
		
	return( Acc );	
};

////////////////////////////////////////////////////////////////////////
double TExplorationAcceptance::cut9( void )
{
	
	return 1.0;
	
};

////////////////////////////////////////////////////////////////////////
double TExplorationAcceptance::cut10( void )
{
		return 1.0;
};

////////////////////////////////////////////////////////////////////////
double TExplorationAcceptance::cut11( void )
{
	double Acc = 1.0;
	
	double ptpip = pf[3].Pt();
	double ptpim = pf[4].Pt();
	
	double etapip = pf[3].PseudoRapidity();
	double etapim = pf[4].PseudoRapidity();
	
	//double etapip = pf[3].Rapidity();
	//double etapim = pf[4].Rapidity();
	
	
	Acc = 0.0;

	
	if( ( ptpip > 0.1 )&&( ptpim > 0.1 )&&( abs(etapim) < 0.9 )&&( abs(etapip) < 0.9 ) )
	{
		Acc = 1.0;
	}
	else
	{
		Acc = 0.0;
	}
	

	
	return( Acc );
};
