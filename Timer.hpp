/* 
 * File:   Timer.hpp
 * Author: einstein
 *
 * Created on 7 de Julho de 2011, 11:03
 */

#ifndef TIMER_HPP
#define	TIMER_HPP
#include <ctime>

class MyTimer {
private:
    time_t sec0;
    clock_t clock0;
public:

    void start() {
        time(&sec0);
        clock0 = clock();
    }

    int getSec() {
        time_t t2;
        time(&t2);
        return difftime(t2, sec0);
    }
    
    double getMillSec() {        
        return ((float)clock() - (float)clock0)/CLOCKS_PER_SEC;
    }
};

MyTimer globalTimer;

#endif	/* TIMER_HPP */

