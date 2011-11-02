#ifndef PRINT_XML_H
#define PRINT_XML_H

#include <json_spirit_value.h>
#include <iostream>
#include <string>

void print_xml(std::ostream &os, const std::string &name,
               const json_spirit::Value &o);

inline void print_xml(std::ostream &os, const json_spirit::Value &o)
{
  print_xml(os,"",o);
}

#endif
