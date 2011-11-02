#include <json_spirit_value.h>
#include <iostream>
#include <string>

#include "print_xml.h"
#include "fix_start_end.h"

void print_xml(std::ostream &os, const std::string &name,
               const json_spirit::Value &o)
{
  std::string start, end;
  switch(o.type())
    {
    case json_spirit::str_type:
      fix_start_end("param",name,start,end);
      os << start << o.get_str() << end;
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
      for(json_spirit::Array::const_iterator i=o.get_array().begin();
          i!=o.get_array().end(); ++i)
        {
          print_xml(os,*i);
        }
      os << end;
      break;
    }
}

