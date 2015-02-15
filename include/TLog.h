/***********************************************************************
				TLog

	Simple class to log into file.


***********************************************************************/

#ifndef TLog_H
#define TLog_H

/*!
  @class TLog

  @brief Simple class to log into file.

 */


#include <iostream>
#include <fstream>
#include <string>
#include <time.h>
#include <map>

using namespace std; 


///Log types
enum TLogLevel {logERROR, logWARNING, logINFO, logCONFIG, logDEBUG, logDEBUG1,logDEBUG2, logDEBUG3, logDEBUG4};
 

////////////////////////////////////////////////////////////////////////
class TLog
{
	
public:
	///	Constructor 
	///	create logger that saves logs into filename
	/// @param logTime  log the time or not
	TLog();
		
	/// Destructor
	virtual ~TLog();
  
	/// Open logfile for writing
	/// @param logfile  file with logs
	void openLogFile( string logfile = string( "Log.txt" ) );
  
	/// close logfile
	void closeLogFile( void );
  
	/// returns string that writes into logfile
	std::ofstream& Write(TLogLevel level = logCONFIG);

	///	returns current system time
	string getTime();
	
	/// set if save time of logs into file
	void setLogTime( bool logTime ){this->logTime = logTime;};
	
	bool isLogTimeSet( void ){ return logTime; };
	
	bool isLogOpen( void ){ return logOpen; };

	/// convert log level number into string name
	string levelToString(TLogLevel level);

protected:
	
	/// stream with logfile
	std::ofstream of;

private:

	/// log the time into file
	bool logTime;
	
	/// flag of opened log file
	bool logOpen;
   
};
#endif
