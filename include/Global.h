/***********************************************************************
	Global.h

Global variables for GenEx.


***********************************************************************/

#ifndef Global_H
#define Global_H

/*!

  @brief Global variables for GenEx.

*/

#include "TLog.h"
#include "TConfigReader.h"
#include "TPolicyReader.h"
#include "Singleton.h"

///Global logging system
typedef Singleton<TLog> LoggerSubsystem;

///Global configuration reader
typedef Singleton<TPolicyReader> ConfigReaderSubsystem;


#endif
