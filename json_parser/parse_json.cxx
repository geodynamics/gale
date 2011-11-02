#include <json_spirit_reader.h>
#include <json_spirit_error_position.h>
#include <iostream>
#include <fstream>
#include <string>
#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>

#include "print_xml.h"

void add_defaults(json_spirit::Value &root);

boost::filesystem::path parse_json(const boost::filesystem::path &filename)
{
  boost::filesystem::path xml_name(filename.stem().string()+"_json.xml");
  try
    {
      json_spirit::Value root;
      boost::filesystem::ifstream infile(filename);
      json_spirit::read_or_throw(infile,root);

      add_defaults(root);

      json_spirit::Object &o(root.get_obj());
      for(json_spirit::Object::iterator i=o.begin(); i!=o.end(); ++i)
        {
          if(i->name_=="outputPath")
            {
              if(i->value_.type()!=json_spirit::str_type)
                throw
                  std::runtime_error("The field 'outputPath' must be a string");
              xml_name=i->value_.get_str() / xml_name;
              break;
            }
        }

      std::ofstream outfile(xml_name.c_str());
      outfile << "<?xml version=\"1.0\"?>\n";
      print_xml(outfile,"StGermainData",root);
    }
  catch (json_spirit::Error_position &e)
    {
      std::cerr << "Error on line "
                << e.line_ << " at column "
                << e.column_ << "\n\t"
                << e.reason_
                << "\n";
      xml_name.clear();
    }
  catch (std::exception &e)
    {
      std::cerr << e.what()
                << "\n";
      xml_name.clear();
    }
  return xml_name;
}
