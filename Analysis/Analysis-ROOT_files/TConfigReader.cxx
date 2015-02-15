/***********************************************************************

					TConfigReader.h 
  
Simple class that reads config file into map (key, value).

	File format:
 
	key1 = value1
	key2 = value2
	...
 
	Lines with comments starts with # or ;

***********************************************************************/

#include "TConfigReader.h"


////////////////////////////////////////////////////////////////////////
bool TConfigReader::IsKey( const std::string& s ) const
{
	return count( s ) != 0;
};
      
////////////////////////////////////////////////////////////////////////
int TConfigReader::GetIntValue( const std::string& key )
{
	if (!IsKey( key ))
	{
		throw 0;
	}	
	std::istringstream ss( this->operator [] ( key ) );
	int result;
	ss >> result;
	if (!ss.eof()) throw 1;
	return result;
};  
      

////////////////////////////////////////////////////////////////////////
float TConfigReader::GetFloatValue( const std::string& key )
{
	if (!IsKey( key )) 
	{
		throw 1.0f;
	}	
	std::istringstream ss( this->operator [] ( key ) );
	float result;
	ss >> result;
	if (!ss.eof()) throw 1;
	return result;
};  
      
////////////////////////////////////////////////////////////////////////
double TConfigReader::GetDoubleValue( const std::string& key )
{
	if (!IsKey( key )) 
	{
		throw 1.0;
	}	
	std::istringstream ss( this->operator [] ( key ) );
	double result;
	ss >> result;
	if (!ss.eof()) throw 1;
	return result;
};        
      
////////////////////////////////////////////////////////////////////////
long double TConfigReader::GetLongDoubleValue( const std::string& key )
{
	if (!IsKey( key )) 
	{
		throw 1.0l;
	}	
	std::istringstream ss( this->operator [] ( key ) );
	long double result;
	ss >> result;
	if (!ss.eof()) throw 1;
	return result;
};       

////////////////////////////////////////////////////////////////////////
std::complex<double> TConfigReader::GetComplexDoubleValue( const std::string& key )
{
	if (!IsKey( key )) 
	{
		throw std::complex<double>(0.0, 0.0);
	}	
	std::istringstream ss( this->operator [] ( key ) );
	std::complex<double> result;
	
	double re;
	double im;
	std::string bufer;
	ss >> bufer;
	sscanf(bufer.c_str(),"(%lf,%lf)",&re,&im);
	result = std::complex<double>(re,im);
	if (!ss.eof()) throw std::complex<double>(0.0, 0.0);
	
	return result;
};       


////////////////////////////////////////////////////////////////////////
bool TConfigReader::GetBoolValue( const std::string& key )
{
	if (!IsKey( key )) 
	{
		throw true;
	}	
	std::istringstream ss( this->operator [] ( key ) );
	bool result;
	ss >> result;
	if (!ss.eof()) throw 1;
	return result;
};        
      
////////////////////////////////////////////////////////////////////////
std::string TConfigReader::GetStringValue( const std::string& key )
{
	if (!IsKey( key )) 
	{
		throw std::string("a");
	}	
	std::istringstream ss( this->operator [] ( key ) );
	std::string result;
	ss >> result;
	if (!ss.eof()) throw 1;
	return result;
};


////////////////////////////////////////////////////////////////////////
void TConfigReader::GetIntArray( int * array, std::string basename, int dim, int start )
{
	assert( dim > 0 );
	assert( start <= dim );
	assert( array != 0 );
	
	std::string base;
	std::string result;
	std::stringstream buffer;
	
	base = basename;
	
	for (int i = start; i < dim ; i++ )
	{
			result.clear();
			buffer.clear();
			
			buffer << base << i;
			buffer >> result;
			array[i] = GetIntValue( result );
	}
	
	return;
	
};

////////////////////////////////////////////////////////////////////////
void TConfigReader::GetFloatArray( float * array, std::string basename, int dim, int start )
{
	assert( dim > 0 );
	assert( start <= dim );
	assert( array != 0 );
	
	std::string base;
	std::string result;
	std::stringstream buffer;
	
	base = basename;
	
	for (int i = start; i < dim ; i++ )
	{
			result.clear();
			buffer.clear();
			
			buffer << base << i;
			buffer >> result;
			array[i] = GetFloatValue( result );
	}
	
	return;
	
};

////////////////////////////////////////////////////////////////////////
void TConfigReader::GetDoubleArray( double * array, std::string basename, int dim, int start )
{
	assert( dim > 0 );
	assert( start <= dim );
	assert( array != 0 );
	
	std::string base;
	std::string result;
	std::stringstream buffer;
	
	base = basename;
	
	for (int i = start; i < dim ; i++ )
	{
			result.clear();
			buffer.clear();
			
			buffer << base << i;
			buffer >> result;
			array[i] = GetDoubleValue( result );
	}
	
	return;
	
};


////////////////////////////////////////////////////////////////////////
void TConfigReader::GetComplexDoubleArray( std::complex<double> * array, std::string basename, int dim, int start )
{
	assert( dim > 0 );
	assert( start <= dim );
	assert( array != 0 );
	
	std::string base;
	std::string result;
	std::stringstream buffer;
	
	base = basename;
	
	for (int i = start; i < dim ; i++ )
	{
			result.clear();
			buffer.clear();
			
			buffer << base << i;
			buffer >> result;
			array[i] = GetComplexDoubleValue( result );
	}
	
	return;
	
};


////////////////////////////////////////////////////////////////////////
void TConfigReader::GetStringArray( std::string * array, std::string basename, int dim, int start )
{
	assert( dim > 0 );
	assert( start <= dim );
	assert( array != 0 );
	
	std::string base;
	std::string result;
	std::stringstream buffer;
	
	base = basename;
	
	for (int i = start; i < dim ; i++ )
	{
			result.clear();
			buffer.clear();
			
			buffer << base << i;
			buffer >> result;
			array[i] = GetStringValue( result );
	}
	
	return;
	
};

////////////////////////////////////////////////////////////////////////
void TConfigReader::ReadConfFile( const std::string& filename )
{
	std::ifstream f( filename.c_str() );
	f >> *this;
	return;
};

////////////////////////////////////////////////////////////////////////
std::istream& operator >> ( std::istream& ins, TConfigReader& d )
{
	std::string s, key, value;

    // For each (key, value) pair in the file
    while (std::getline( ins, s ))
    {
		std::string::size_type begin = s.find_first_not_of( " \f\t\v" );

		// Skip blank lines
		if (begin == std::string::npos) continue;

		// Skip commentary
		if (std::string( "#;" ).find( s[ begin ] ) != std::string::npos) continue;

		// Extract the key value
		std::string::size_type end = s.find( '=', begin );
		key = s.substr( begin, end - begin );

		// (No leading or trailing whitespace allowed)
		key.erase( key.find_last_not_of( " \f\t\v" ) + 1 );

		// No blank keys allowed
		if (key.empty()) continue;

		// Extract the value (no leading or trailing whitespace allowed)
		begin = s.find_first_not_of( " \f\n\r\t\v", end + 1 );
		end   = s.find_first_of(  " \f\n\r\t\v#;", begin );

		value = s.substr( begin, end - begin );		

		// Checking if the value is already read
		if ( d.IsKey( key ) ) 
		{
			std::cout << "TConfigReader::WARNING: Replacing value of the variable " << key << std::endl;
		}	
		
		// Insert the properly extracted (key, value) pair into the map
		d[ key ] = value;
    }
	return ins;
};


////////////////////////////////////////////////////////////////////////
std::ostream& operator << ( std::ostream& outs, const TConfigReader& d )
{
	TConfigReader::const_iterator iter;
    for (iter = d.begin(); iter != d.end(); iter++)
		outs << iter->first << " = " << iter->second << std::endl;
	return outs;
};
