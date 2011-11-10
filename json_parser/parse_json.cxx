#include <json_spirit_reader.h>
#include <json_spirit_error_position.h>
#include <iostream>
#include <string>
#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>
#include <sstream>

#include "print_xml.h"

void add_defaults(json_spirit::Value &root);

std::string parse_json(const boost::filesystem::path &filename)
{
  std::stringstream result;
  try
    {
      json_spirit::Value root;
      boost::filesystem::ifstream infile(filename);
      json_spirit::read_or_throw(infile,root);
      add_defaults(root);

      result << "<?xml version=\"1.0\"?>\n";
      print_xml(result,"StGermainData",root);
    }
  catch (json_spirit::Error_position &e)
    {
      std::cerr << "Error on line "
                << e.line_ << " at column "
                << e.column_ << "\n\t"
                << e.reason_
                << "\n";
    }
  catch (std::exception &e)
    {
      std::cerr << e.what()
                << "\n";
    }
  return result.str();
}
