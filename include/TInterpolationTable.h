/***********************************************************************
			
			TInterpolationTable
  
  Class that realizes interpolation according to given tables of nodes of x and y(x).


***********************************************************************/

#ifndef TInterpolationTable_H
#define TInterpolationTable_H

/*!
  @class TInterpolationTable

  @brief Class that realizes interpolation according to given tables of nodes of x and y(x).

 */

#include <iostream> 

#include <utility>      
#include <vector>
#include <algorithm>
#include <string>
#include <assert.h>


using namespace std;


///Node of interpolation
typedef pair <double,double> node;

/*
///sorting predicate according to first element of node
static bool sort_pred_first(const node& left, const node& right)
{
    return left.first < right.first;
};
*/
/*
///sorting predicate according to second element of node
static bool sort_pred_second(const node& left, const node& right)
{
    return left.second < right.second;
};
*/

////////////////////////////////////////////////////////////////////////
class TInterpolationTable
{
private:

	/// Nodes vector
	vector< node > _vNodes;
	
	/// Initialization flag
	bool _bInitialized;
	
	/// Minimal y value in nodes
	double y_min;
	
	///Maximal y value in nodes
	double y_max;
	
	///sorting predicate according to first element of node
	struct sort_pred_first
    {
		bool operator () (const node& left, const node& right)
        {
			return left.first < right.first;
        }
    };
    
    ///sorting predicate according to second element of node
    struct sort_pred_second
    {
		bool operator () (const node& left, const node& right)
        {
			return left.second < right.second;
        }
     };
	

public:
	
	///constructor
	TInterpolationTable( void );
	
	///destructor
	virtual ~TInterpolationTable( void );

	/// Set interpolation tables/nodes y(x) of N elements
	void SetNodes( double * x, double * y, int N );
	
	/// Remove nodes and all interpolation information
	void DeleteNodes( void );
	
	/// @returns interpolated value of y(x)
	double Interpolate( double x );
	
};



#endif
