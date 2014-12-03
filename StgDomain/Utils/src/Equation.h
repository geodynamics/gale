#ifndef __StGermain_Utils_Equation_h__
#define __StGermain_Utils_Equation_h__

#include <string>

double Equation_eval(const double *coord, DomainContext *context,
                     const std::string &equation, const bool &add_dt);

inline double Equation_eval(const double *coord, DomainContext *context,
                            const std::string &equation)
{
  return Equation_eval(coord,context,equation,false);
}

#endif
