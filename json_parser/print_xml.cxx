#include <json_spirit_value.h>
#include <iostream>
#include <string>

#include "print_xml.h"
#include "fix_start_end.h"

std::string fix_comparisons(const std::string &s)
{
  std::string result(s);
  size_t i(result.find('<'));
  while(i!=std::string::npos)
    {
      result.replace(i,1,"&lt;");
      i=result.find('<',i+1);
    }
  i=result.find('>');
  while(i!=std::string::npos)
    {
      result.replace(i,1,"&gt;");
      i=result.find('>',i+1);
    }
  i=result.find('&');
  while(i!=std::string::npos)
    {
      result.replace(i,1,"&amp;");
      i=result.find('&',i+1);
    }
  return result;
}

void print_xml(std::ostream &os, const std::string &name,
               const json_spirit::Value &o)
{
  std::string start, end;
  switch(o.type())
    {
    case json_spirit::str_type:
      fix_start_end("param",name,start,end);
      os << start << fix_comparisons(o.get_str()) << end;
      break;
    case json_spirit::bool_type:
      fix_start_end("param",name,start,end);
      os << start << o.get_bool() << end;
      break;
    case json_spirit::int_type:
      fix_start_end("param",name,start,end);
      os << start << o.get_int() << end;
      break;
    case json_spirit::real_type:
      fix_start_end("param",name,start,end);
      os << start << o.get_real() << end;
      break;
    case json_spirit::null_type:
      fix_start_end("param",name,start,end);
      os << start << end;
      break;
    case json_spirit::obj_type:
      fix_start_end("struct",name,start,end);
      os << start;
      for(json_spirit::Object::const_iterator i=o.get_obj().begin();
          i!=o.get_obj().end(); ++i)
        {
          print_xml(os,i->name_,i->value_);
        }
      os << end;
      break;
    case json_spirit::array_type:
      fix_start_end("list",name,start,end);
      os << start;
      json_spirit::Array::const_iterator i=o.get_array().begin();
      if(o.get_array().size()>1 && i->type()==json_spirit::str_type
         && i->get_str()=="asciidata")
        {
          os << "<asciidata>\n";
          ++i;
          for(json_spirit::Array::const_iterator j=i->get_array().begin();
              j!=i->get_array().end(); ++j)
            {
              os << "<columnDefinition name=\"" << j->get_str()
                 << "\" type=\"double\"/>\n";
            }
          ++i;
          for(;i!=o.get_array().end(); ++i)
            {
              if(i->type()==json_spirit::int_type)
                {
                  os << i->get_int() << " ";
                }
              else
                {
                  os << i->get_real() << " ";
                }
            }
          os << "</asciidata>\n";
        }
      else
        {
          for(;i!=o.get_array().end(); ++i)
            {
              print_xml(os,*i);
            }
        }
      os << end;
      break;
    }
}

