
#ifndef DEBUG_INCLUDE
#define DEBUG_INCLUDE

#include <fstream>
#include <iostream>
#include <vector>
#include "TrCluster.hh"
#include "vector2.h"
using namespace std;

namespace debug {

  //output log file name
  const char *filename;

  //output log file stream
  fstream out;

  void start_debug(const char *name = "debug.log") {
    out.open(name,ios_base::out);
    if(!out.is_open()) {
      cerr <<"[DEBUG] fatal error: cannot open output file";
      exit(1);
    }
    filename = name;
  }

  inline void end_debug()
  {
    out.close();
    std::cout <<"\n[DEBUG]: output written in " <<filename <<"\n\n";
  }

  void print_cl(const TrCluster *cl) {
    if(!out.is_open())
      return;

    out <<"layer:\t" <<cl->layer <<endl
      <<"ladder:\t" <<cl->ladder <<endl
      <<"strip:\t" <<cl->strip <<endl
      <<"position:\t(" <<cl->pos[0] <<"," <<cl->pos[1] <<","
      <<cl->pos[2] <<")\n"
      <<"clust:\t(" <<cl->clust[0] <<"," <<cl->clust[1] <<")\n";
  }

  void print_vec2(const char *name, const vector2<double> &v) {
    if(!out.is_open())
      return;

    out <<endl <<name <<" {\n";

    for(int i = 0; i < (int) v.size(); ++i)
      for(int j = 0; j < (int) v[i].size(); ++j)
        if(v[i][j] != 0)
          out <<i <<":" <<j <<"\t" <<v[i][j] <<endl;

    out <<"}\n\n";
  }

  void print_vecpair(const char *name, const vector_pair<double> &v) {
    if(!out.is_open())
      return;

    out <<endl <<name <<" {\n";

    for(int i = 0; i < (int) v.size(); ++i)
      out <<i <<":\t(" <<v[i].first <<"," <<v[i].second <<")\n";

    out <<"}\n\n";
  }
}

#endif
