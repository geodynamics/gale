#include <json_spirit_reader.h>
#include <stdexcept>

class Compare_Name
{
public:
  std::string name;
  Compare_Name(const std::string &n): name(n) {}
  bool operator()(json_spirit::Pair p)
  {
    return p.name_==name;
  }
};

void add_components(json_spirit::Object &o,
                    const json_spirit::Object &components,
                    const std::string &enable_key, const bool &default_enable)
{
  bool enable(default_enable);
  for(json_spirit::Object::const_iterator i=o.begin(); i!=o.end(); ++i)
    {
      if(i->name_==enable_key)
        {
          if(i->value_.type()!=json_spirit::bool_type)
            throw std::runtime_error(("The field " + enable_key +
                                      "must be either true or false (no quotes)").c_str());
          enable=i->value_.get_bool();
          break;
        }
    }
  json_spirit::Object temp;
  if(enable)
    {
      for(json_spirit::Object::const_iterator i=components.begin();
          i!=components.end(); ++i)
        {
          if(std::find_if(o.begin(),o.end(),Compare_Name(i->name_))==o.end())
            {
              temp.push_back(*i);
            }
        }
      temp.insert(temp.end(),o.begin(),o.end());
      temp.swap(o);
    }
}
