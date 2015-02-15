/***********************************************************************

	Singleton.h

Simple singelton wrapper.


Example use:
	
	//create logges subsystem in global scope
	typedef Singleton<Logger> LoggerSubsystem;
	
	//use logger system (specific for Logger class)
	LoggerSubsystem::Instance()->openLogFile("LogFile.txt");
	
	//delete log system
	LoggerSubsystem::Destroy();

***********************************************************************/


#ifndef SINGLETON_H_
#define SINGLETON_H_

/*!
  @class Singleton

  @brief Simple singelton wrapper.
  
  Example use:
	
	//create logges subsystem in global scope
	typedef Singleton<Logger> LoggerSubsystem;
	
	//use logger system (specific for Logger class)
	LoggerSubsystem::Instance()->openLogFile("LogFile.txt");
	
	//delete log system
	LoggerSubsystem::Destroy();

*/

 
#include <iostream>
 
template <class T>
class Singleton
{
	public:
		// Return existing or create new instance
		static T* Instance()
		{
			// Do we have an instance of this type? If so return it, otherwise create a new one.
			return _pInstance ? _pInstance : _pInstance = new T;
		}
 
		// Manually destroy an existing instance. Call at end of program to clean up.
		static void Destroy()
		{
			delete _pInstance;
			_pInstance = NULL;
		}
 
	private:
		Singleton();                            // Constructor                   (empty & cannot be called externally)
		~Singleton();                           // Destructor                    (empty & cannot be called externally)
		Singleton(Singleton const&);            // Copy constructor              (empty & cannot be called externally - no copies allowed)
		Singleton& operator=(Singleton const&); // Assignment operator           (empty & cannot be called externally - no assignment allowed)
		static T* _pInstance;                  // Static template-type instance
};
 
// Set static instance value to NULL
template <class T> T* Singleton<T>::_pInstance = NULL;
 
#endif
