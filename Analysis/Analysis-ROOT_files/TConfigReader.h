/***********************************************************************

					TConfigReader.h 
  
  Simple class that reads config file into map (key, value).

	File format:
 
	key1 = value1
	key2 = value2
	...
 
	Lines with comments starts with # or ;


***********************************************************************/

#ifndef TCONFIGREADER_H
#define TCONFIGREADER_H

/*!
  @class TConfigReader

  @brief  Simple class that reads config file into map (key, value).

	File format:
 
	key1 = value1
	key2 = value2
	...
 
	Lines with comments starts with # or ;

 */


#include <iostream>
#include <map>
#include <string>
#include <fstream> 
#include <sstream>
#include <assert.h>
#include <cmath>
#include <math.h>
#include <complex>

////////////////////////////////////////////////////////////////////////
class TConfigReader: public std::map <std::string, std::string>
{
public:	
		
	/// Constructor
	TConfigReader() {};
	
	/// Destructor
	virtual ~TConfigReader() {};
		
		
    /// Check if key is in the map
	virtual bool IsKey( const std::string& s ) const;
	
	/// Reads configuration form given filename
	virtual void ReadConfFile( const std::string& filename );
      
	///Gets an integer value from a key. If the key does not exist, or if the value is not an integer, throws an int exception.
	virtual int GetIntValue( const std::string& key );
      

	///Gets a real value from a key. If the key does not exist, or if the value is not a real, throws an real exception.
	virtual float GetFloatValue( const std::string& key );  
      
	///Gets a double value from a key. If the key does not exist, or if the value is not a double, throws a double exception.
	virtual double GetDoubleValue( const std::string& key );        
      
	///Gets a long double value from a key. If the key does not exist, or if the value is not a long double, throws a long double exception.
	virtual long double GetLongDoubleValue( const std::string& key );       

	///Gets a complex<double> value from a key. If the key does not exist, or if the value is not a complex<double>, throws a complex<double> exception.
	virtual std::complex<double> GetComplexDoubleValue( const std::string& key ); 

	///Gets a boolean value from a key. If the key does not exist, or if the value is not a boolean, throws a boolean exception.
	virtual bool GetBoolValue( const std::string& key );        
      
	///Gets a string value from a key. If the key does not exist, or if the value is not a string, throws a string exception.
	virtual std::string GetStringValue( const std::string& key );
	
	///Gets array of int values from variables basename+number, where number runs form [start, dim-1]
	virtual void GetIntArray( int * array, std::string basename, int dim, int start = 0 );
	
	///Gets array of float values from variables basename+number, where number runs form [start, dim-1]
	virtual void GetFloatArray( float * array, std::string basename, int dim, int start = 0 );

	///Gets array of double values from variables basename+number, where number runs form [start, dim-1]
	virtual void GetDoubleArray( double * array, std::string basename, int dim, int start = 0 );
	
	///Gets array of complex<double> values from variables basename+number, where number runs form [start, dim-1]
	virtual void GetComplexDoubleArray( std::complex<double> * array, std::string basename, int dim, int start = 0 );

	///Gets array of strings values from variables basename+number, where number runs form [start, dim-1]
	virtual void GetStringArray( std::string * array, std::string basename, int dim, int start = 0 );
	  
};


/// The extraction operator reads configuration data until EOF. Invalid data is ignored.
std::istream& operator >> ( std::istream& ins, TConfigReader& d );


/// The insertion operator writes all configuration::data to stream.
std::ostream& operator << ( std::ostream& outs, const TConfigReader& d );

#endif
