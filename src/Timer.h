
#include <sys/time.h>
#include <iostream>

class Timer
{
    double tot;
    std::string label;

    timeval tp;
    double t1;
    double t2;

public :
    Timer(std::string _label) : tot(0), label(_label) {}
    ~Timer() { std::cout<<"Timer "<<label<<" = "<<tot<<std::endl; }

    void Start()
    {
        gettimeofday(&tp,0);
        t1 = tp.tv_sec + tp.tv_usec/1.e6;
    }

    void Stop()
    {
        gettimeofday(&tp,0);
        t2 = tp.tv_sec + tp.tv_usec/1.e6;
        tot += t2-t1;
    }
};


