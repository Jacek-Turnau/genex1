/***********************************************************************
			
				TInterpolationTable
  
 Class that realizes interpolation according to given tables of nodes of x and y(x).

***********************************************************************/


#include "TInterpolationTable.h"

////////////////////////////////////////////////////////////////////////
TInterpolationTable::TInterpolationTable() :   _bInitialized( false ), y_min(0.0), y_max(0.0)
{

	
	return;
};


////////////////////////////////////////////////////////////////////////
TInterpolationTable::~TInterpolationTable()
{
	return;
};


////////////////////////////////////////////////////////////////////////
void TInterpolationTable::SetNodes( double * x, double * y, int N )
{
	assert( N > 1 );
	
	//fill in table of nodes
	y_min = y[0];
	y_max = y[0];
	for( int i = 0; i < N; i++ )
	{
		_vNodes.push_back( node( x[i], y[i] ) );
		
		//set up min value
		if( y_min > y[i])
			y_min = y[i];
			
		if( y_max < y[i] )
			y_max = y[i];
	}
	
	//sort table according to x
	sort(_vNodes.begin(), _vNodes.end(), sort_pred_first());
	
	_bInitialized = true;
	
};


////////////////////////////////////////////////////////////////////////
void TInterpolationTable::DeleteNodes( void )
{

	_vNodes.clear();
	
	_bInitialized = false;	
	
};

////////////////////////////////////////////////////////////////////////
double TInterpolationTable::Interpolate( double x )
{
	if( _bInitialized == false )
	{
		cerr << "ERROR - set nodes first!" << endl;
		throw string( "ERROR - set nodes first!" );
	}
	
	
	// Assumes that _vNodes is sorted by the first elements
    // Check if x is out of bound
    if (x > _vNodes.back().first)
	{
		cerr << "ERROR - x not in  the interpolating interval!" << endl;
		throw string( "ERROR - x not in  the interpolating interval!" );
	}
    
    if (x < _vNodes[0].first)
    {
		cerr << "ERROR - x not in  the interpolating interval!" << endl;
		throw string( "ERROR - x not in  the interpolating interval!" );
	}
    
    vector< node >::iterator it1, it2;
  
    
    // extract iterator to smaller but not greather than x value
    it1 = lower_bound( _vNodes.begin(), _vNodes.end(), node(x, y_min), sort_pred_first() );
    
    // if iterator is at the boundary
    if (it1 == _vNodes.begin()) 
		return it1->second;
    
    // extract iterator to second element
    it2 = it1;
    --it2;
    
    //linear interpolation 
    double y_interp = it2->second + (it1->second - it2->second)*(x - it2->first)/(it1->first - it2->first);
	
	return( y_interp );	
	
};
