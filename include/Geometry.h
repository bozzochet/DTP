
// read geometry parameters set in macros/geo.mac

#ifndef GEOMETRY_INCLUDE
#define GEOMETRY_INCLUDE


#include <string>
#include <cstring>
#include <map>
#include <fstream>
#include <iostream>


class Geometry
{
  const std::string filename_ = "macros/geo.mac";

  int Nlayers_ = -1;
  int Nsquares_ = -1;
  int Nrows_ = -1;
  int Nstrips_ = -1;
  double pitch_ = -1;
  double thickness_ = -1;

  //name and references to class variables
  const std::map <std::string, void*> map_name_ = {
    {"Nlayers", (void*) &Nlayers_},
    {"Nstrips", (void*) &Nstrips_},
    {"Nsquares", (void*) &Nsquares_},
    {"Nrows", (void*) &Nrows_},
    {"pitch", (void*) &pitch_},
    {"thickness", (void*) &thickness_}
  };


  //return line with a value
  std::string InterestingLine(std::fstream&);

  //get value name
  std::string GetName(const std::string&);

  //get value as a string
  std::string GetValueStr(const std::string&);

  //get reference to class variable corresponding to variable read
  void* GetValueRef(const std::string&);

public:
  Geometry();

  inline int GetNlayers()
  { return Nlayers_; }

  inline int GetNsquares()
  { return Nsquares_; }

  inline int GetNrows()
  { return Nrows_; }

  inline int GetNstrips()
  { return Nstrips_; }

  inline double GetPitch()
  { return pitch_; }

  inline double GetThickness()
  { return thickness_; }

  inline double GetSquareSide()
  { return pitch_ * Nstrips_; }

  inline int GetNladders()
  { return Nsquares_ * Nrows_ * Nlayers_; }

  void Print();
};


#endif //include guard
