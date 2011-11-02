#include <algorithm>
#include <list>
#include "fix_start_end.h"

void fix_start_end(const std::string &default_name,
                   const std::string &name,
                   std::string &start, std::string &end)
{
  /* Horribly inefficient */
  std::list<std::string> reserved_words;
  reserved_words.push_back("struct");
  reserved_words.push_back("import");
  reserved_words.push_back("plugins");
  reserved_words.push_back("asciidata");
  reserved_words.push_back("columnDefinition");
  reserved_words.push_back("toolbox");

  if(std::find(reserved_words.begin(),
               reserved_words.end(), name)!=reserved_words.end())
    {
      start="<" + name + ">";
      end="</" + name + ">\n";
    }
  else if(name=="")
    {
      start="<" + default_name + ">";
      end="</" + default_name + ">\n";
    }
  else if(name=="StGermainData")
    {
      start="<StGermainData xmlns=\"http://www.vpac.org/StGermain/XML_IO_Handler/Jun2003\">\n";
      end="</StGermainData>\n";
    }
  else
    {
      start="<" + default_name + " name=\"" + name + "\">";
      end="</" + default_name + ">\n";
    }
}
