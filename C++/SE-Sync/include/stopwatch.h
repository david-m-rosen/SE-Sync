/** This file provides a convenient functional interface to the SESync algorithm
*
* jbriales 09 Sept 2017
*/

#ifndef _STOPWATCH_H_
#define _STOPWATCH_H_

#include <chrono>
#include <map>

#include <string>
#include <iostream>
#include <fstream>

namespace SESync
{
  namespace util
  {

    /** This class provides an intuitive and compact interface to measure
     * execution time of sections of code */
    class stopwatch
    {
    private:
      std::chrono::time_point<std::chrono::high_resolution_clock> start_time;
      std::chrono::duration<double> counter;
      double elapsed_time;

    public:
      /** Begin counting time */
      inline void start()
      {
        start_time = std::chrono::high_resolution_clock::now();
      }

      /** Store current elapsed time since start() */
      inline void stop()
      {
        counter = std::chrono::high_resolution_clock::now() - start_time;
        // elapsed_time = std::chrono::duration_cast<std::chrono::milliseconds>(counter).count() / 1e3;
        elapsed_time = counter.count();
      }

      /** Read elapsed time from start() until last stop() */
      inline double time() const {return elapsed_time;}

      /** Formated verbose output */
      inline void print()
      {
        std::cout << "elapsed computation time: "
                  << elapsed_time << " seconds"
                  << std::endl;
      }
    };

    /** This class provides the ability to define and handle multiple stopwatches
     * simultaneously in a neat way */
    class stopwatch_collection
    {
    private:
      /** Map of all name-stopwatch* pairs defined in the collection */
      std::map<std::string,stopwatch*> watches;

    public:
      stopwatch_collection(){}
      ~stopwatch_collection()
      {
        // Deallocate memory for all created watches
        for(auto ptr = watches.begin(); ptr != watches.end(); ptr++)
          delete ptr->second;
      }

      /** Create new stopwatch in the collection.
       * The pointer to the new watch is returned. */
      inline stopwatch* add(const std::string& name)
      {
        stopwatch* ptr = new stopwatch();
        watches[name] = ptr;
        return ptr;
      }

      /** Look for watch by name and return a pointer to it */
      inline stopwatch* operator[](const std::string& name)
      {
        return watches[name];
      }

      /** Returns a map which contains only the measured times */
      std::map<std::string,double> readTimes()
      {
        std::map<std::string,double> timesMap;
        for(auto ptr = this->watches.begin(); ptr != this->watches.end(); ptr++)
        {
          const std::string& name = ptr->first;
          const stopwatch* sw = ptr->second;
          // Allocate elements in the name-times map
          timesMap[name] = sw->time();
        }
        return timesMap;
      }
    };

  } /* namespace util */
} /* namespace SESync */

#endif /* _STOPWATCH_H_ */
