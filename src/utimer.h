#include <iostream>
#include <iomanip>
#include <chrono>

#ifndef utimer_H
#define utimer_H


#define START(timename) auto timename = std::chrono::system_clock::now();
#define STOP(timename,elapsed)  auto elapsed = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::system_clock::now() - timename).count();


class utimer {
  std::chrono::system_clock::time_point start;
  std::chrono::system_clock::time_point stop;
  std::string message; 
  using usecs = std::chrono::microseconds;
  using msecs = std::chrono::milliseconds;

private:
  long * us_elapsed;
  
public:

  utimer() : us_elapsed((long *)NULL) {
    start = std::chrono::system_clock::now();
  }
    
  utimer(long * us) : us_elapsed(us) {
    start = std::chrono::system_clock::now();
  }

  long getElapsed() {
    stop = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed = stop - start;
    return std::chrono::duration_cast<std::chrono::microseconds>(elapsed).count();
  }

  ~utimer() {
    long musec = getElapsed();
    
    if (us_elapsed != NULL)
      (*us_elapsed) = musec;
  }
};

#endif
