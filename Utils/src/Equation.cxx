/* We have to include mpParser first because StGermin #defines
   GetParent, which is a function used in mpParser */

#include "mpParser.h"
#include "StGermain/StGermain.h"
#include "StgDomain/StgDomain.h"

double Equation_eval(const double *coord, DomainContext *context,
                     const std::string &equation, const bool &add_dt)
{
  mup::Value result;
  try
    {
      mup::ParserX p(mup::pckALL_NON_COMPLEX);

      p.DefineConst("x", coord[0]);
      p.DefineConst("y", coord[1]); 
      if(context->dim==3)
        p.DefineConst("z", coord[2]); 
      p.DefineConst("t", context->currentTime+(add_dt ? context->dt : 0));
      p.EnableAutoCreateVar(true);
      p.SetExpr(equation);

      result=p.Eval();

      if(context->dim==2)
        Journal_PrintfL(Journal_Register( Info_Type,"Equation"),
                        2, "Equation %s:  x=%g y=%g t=%g result=%g\n",
                        equation.c_str(),coord[0],coord[1],
                        context->currentTime+(add_dt ? context->dt : 0),
                        result.GetFloat());
      else
        Journal_PrintfL(Journal_Register( Info_Type,"Equation"),
                        2, "Equation %s:  x=%g y=%g z=%g t=%g result=%g\n",
                        equation.c_str(),coord[0],coord[1],coord[2],
                        context->currentTime+(add_dt ? context->dt : 0),
                        result.GetFloat());
    }
  catch (mup::ParserError &e)
    {
      Journal_Firewall(false,
                       Journal_Register( Error_Type,"Equation"),
                       "Error when parsing equation: (%s)\n\t%s\n",
                       equation.c_str(),e.GetMsg().c_str());
    }
  return result.GetFloat();
}
