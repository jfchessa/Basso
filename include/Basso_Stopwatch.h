#ifndef _BASSO_STOPWATCH_H_
#define _BASSO_STOPWATCH_H_

#include <iostream>
#include <time.h>


namespace Basso
{
    
using namespace std;

/**
\brief A simple timer class

This is a simpel class to emulate a wall clock timer/stopwatch

*/
class Basso_Stopwatch
{
public:

    typedef double time_type;
     
    Basso_Stopwatch() : running_(false) { }
  
    time_type Start( ) { startTime_ = clock(); running_=true; return 0.0; }
  
    time_type GetTime( ) const 
    { 
        if ( running_ )
            return ((time_type)( clock() - startTime_ ))/CLOCKS_PER_SEC;
        return 0.0;
    }
    
    time_type Stop( ) 
    { 
        if ( !running_ ) return 0.0;
        
        running_=false; 
        return ((time_type)( clock() - startTime_ ))/CLOCKS_PER_SEC;
    }
  
private:
    time_t startTime_;
    bool running_;
    
};

std::ostream &operator << ( std::ostream &out, const Basso_Stopwatch &A )
{
	out << A.GetTime();
	return out;
}

} // end of namespace
#endif

