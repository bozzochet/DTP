
#include "Geometry.h"


Geometry::Geometry()
{
  std::fstream file;
  file.open(filename_.c_str(), std::ios::in);

  if( !file.is_open())
  {
    std::cerr <<"[GEOMETRY] fatal error: unable to open " <<filename_
      <<std::endl;
    exit(-1);
  }

  int found = 0; //variables found

  while(found < 5)
  {
    //search line with values and get values

    std::string line = InterestingLine(file);
    std::string name = GetName(line);
    std::string str_value = GetValueStr
    (
      line.substr
        (11 + name.length(), line.length() - 11 - name.length())
    );

    void *value = GetValueRef(name);

    //translate void* to right type

    if( std::strcmp(name.c_str(), "pitch") == 0)
    {
      double *pitch = (double*) value;
      *pitch = std::stod(str_value);
    }
    else
    {
      int *val_int = (int*) value;
      *val_int = std::stoi(str_value);
    }

    ++found;
  }

  file.close();
}


std::string Geometry::InterestingLine(std::fstream &file)
{
  while(true)
  {
    std::string line;
    getline(file, line);

    if
    (
      line.length() > 10
      && std::strcmp(line.substr(1,8).c_str(), "Detector") == 0
    )
      return line;
  }
}


std::string Geometry::GetName(const std::string &line)
{
  std::string name;

  for(int i = 10; i < (int) line.length(); ++i) //read name
  {
    const char c = line[i];

    if(c == ' ')
      break;

    name += c;
  }

  return name;
}


void* Geometry::GetValueRef(const std::string &name)
{
  for(auto it = map_name_.begin(); it != map_name_.end(); ++it)
    if( std::strcmp(name.c_str(), it->first.c_str()) == 0)
      return it->second;

  std::cerr <<"[GEOMETRY] fatal error: unable to correspond " <<name
    <<std::endl;
  exit(-1);
}


std::string Geometry::GetValueStr(const std::string &line)
{
  std::string str_value;

  for(int i = 0; i < (int) line.length(); ++i)
  {
    const char c = line[i];

    if(c == ' ' && str_value.length() == 0)
      continue;
    else if (c == ' ')
      break;
    else
      str_value += c;
  }

  return str_value;
}


void Geometry::Print()
{
  for(auto it = map_name_.begin(); it != map_name_.end(); ++it)
  {
    std::cout << it->first << ": ";

    if( std::strcmp(it->first.c_str(), "pitch") == 0 )
      std::cout << *((double*) it->second);
    else
      std::cout << *((int*) it->second);

    std::cout <<std::endl;
  }
}
