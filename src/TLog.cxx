/***********************************************************************
			TLog.cxx

	Simple class to log into file.

***********************************************************************/
#include "TLog.h"


////////////////////////////////////////////////////////////////////////
TLog :: TLog( ) : logTime(true)
{	
	//	
};

////////////////////////////////////////////////////////////////////////
std::ofstream& TLog::Write( TLogLevel level )
{
	if( logOpen )
	{
	
		if( logTime )
		{
			of << "#> " << getTime() << ": " << endl;
		}
		if( level != logCONFIG )
		{	
			of << "#>> " << levelToString( level ) << ": " << endl;;
		}
   }
   
   return of;
};

////////////////////////////////////////////////////////////////////////
TLog::~TLog()
{
	if( logOpen )
	{
		closeLogFile();
	}
};

////////////////////////////////////////////////////////////////////////
string TLog::getTime()
{
		std::string time_str;
		time_t raw_time;
		time( & raw_time );
		time_str = ctime( &raw_time );
		return time_str.substr( 0 , time_str.size() - 1 );
};


////////////////////////////////////////////////////////////////////////

string TLog::levelToString(TLogLevel level)
{
	string tmp;
	
	switch (level)
	{
		case(logERROR)  : tmp = "logERROR";
						break;
		case(logWARNING): tmp = "logWARNING";
						break;
		case(logINFO)   : tmp = "logINFO";
						break;
		case(logCONFIG) : tmp = "logCONFIG";
						break;
		case(logDEBUG)  : tmp = "logDEBUG";
						break;
		case(logDEBUG1) : tmp = "logDEBUG1";
						break;
		case(logDEBUG2) : tmp = "logDEBUG2";
						break;
		case(logDEBUG3) : tmp = "logDEBUG3";
						break;
		case(logDEBUG4) : tmp = "logDEBUG4";
						break;
	}
	
	return tmp;
}

////////////////////////////////////////////////////////////////////////
void TLog::openLogFile( string logfile )
{
	of.open (logfile.c_str(), ios::ate | ios::out | ios::app );
	
	if( logTime )
	{
		of << "#> " << getTime() << " - START OF LOG." << endl;
	}
	
	logOpen = true;
}
  
////////////////////////////////////////////////////////////////////////
void TLog::closeLogFile( void )
{
	if( logOpen )
	{
	
		if( logTime )
		{
			of << "#> " << getTime() << " - END OF LOG." << endl;
		}
	
	
		of.flush();
		of.close();
	}
	
}


