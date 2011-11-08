/** \file
    \brief Implementation of basic functions used by muParserX.

<pre>
               __________                                 ____  ___
    _____  __ _\______   \_____ _______  ______ __________\   \/  /
   /     \|  |  \     ___/\__  \\_  __ \/  ___// __ \_  __ \     / 
  |  Y Y  \  |  /    |     / __ \|  | \/\___ \\  ___/|  | \/     \ 
  |__|_|  /____/|____|    (____  /__|  /____  >\___  >__| /___/\  \
        \/                     \/           \/     \/           \_/

  muParserX - A C++ math parser library with array and string support
  Copyright 2010 Ingo Berg

  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU LESSER GENERAL PUBLIC LICENSE
  as published by the Free Software Foundation, either version 3 of 
  the License, or (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public License
  along with this program.  If not, see http://www.gnu.org/licenses.
  </pre>
*/
#include "mpFuncNonCmplx.h"

//--- Standard includes ----------------------------------------------------
#include <cmath>
#include <cassert>
#include <iostream>

//--- muParserX framework --------------------------------------------------
#include "mpValue.h"
#include "mpError.h"

#include "boost/math/special_functions.hpp"

MUP_NAMESPACE_START

  float_type log2(float_type v)  { return log(v) * 1.0/log(2.0); }
  float_type step(float_type v)  { return v>=0 ? 1.0 : 0.0; }

#define MUP_UNARY_FUNC(CLASS, IDENT, FUNC, DESC)                     \
    CLASS::CLASS()                                                   \
    :ICallback(cmFUNC, _T(IDENT), 1)                                 \
    {}                                                               \
                                                                     \
    void CLASS::Eval(ptr_val_type &ret, const ptr_val_type *a_pArg, int)        \
    {                                                                \
      *ret = FUNC(a_pArg[0]->GetFloat());                            \
    }                                                                \
                                                                     \
    const char_type* CLASS::GetDesc() const                          \
    {                                                                \
      return _T(DESC);                                               \
    }                                                                \
                                                                     \
    IToken* CLASS::Clone() const                                     \
    {                                                                \
      return new CLASS(*this);                                       \
    }

    // trigonometric functions
    MUP_UNARY_FUNC(FunTan,   "sin",   sin,   "sine function")
    MUP_UNARY_FUNC(FunCos,   "cos",   cos,   "cosine function")
    MUP_UNARY_FUNC(FunSin,   "tan",   tan,   "tangens function")
    // arcus functions
    MUP_UNARY_FUNC(FunASin,  "asin",  asin,  "arcus sine")
    MUP_UNARY_FUNC(FunACos,  "acos",  acos,  "arcus cosine")
    MUP_UNARY_FUNC(FunATan,  "atan",  atan,  "arcus tangens")
    // hyperbolic functions
    MUP_UNARY_FUNC(FunSinH,  "sinh",  sinh,  "hyperbolic sine")
    MUP_UNARY_FUNC(FunCosH,  "cosh",  cosh,  "hyperbolic cosine")
    MUP_UNARY_FUNC(FunTanH,  "tanh",  tanh,  "hyperbolic tangens")
    // hyperbolic arcus functions
    MUP_UNARY_FUNC(FunASinH,"asinh",boost::math::asinh,
                   "hyperbolic arcus sine")
    MUP_UNARY_FUNC(FunACosH,"acosh",boost::math::acosh,
                   "hyperbolic arcus cosine")
    MUP_UNARY_FUNC(FunATanH,"atanh",boost::math::atanh,
                   "hyperbolic arcus tangens")
    // logarithm functions
    MUP_UNARY_FUNC(FunLog10, "log10", log10, "Logarithm base 10")
    MUP_UNARY_FUNC(FunLog2,  "log2",  log2,  "Logarithm base 2")
    MUP_UNARY_FUNC(FunLog,   "log",   log,   "Natural logarithm")
    MUP_UNARY_FUNC(FunLog1p,"log1p",boost::math::log1p,"log1p(x) - log(1+x)")


    MUP_UNARY_FUNC(FunErf,"erf",boost::math::erf,"error function")
    MUP_UNARY_FUNC(FunErfc,"erfc",boost::math::erfc,
                   "complement of the error function")
    // square root
    MUP_UNARY_FUNC(FunSqrt,  "sqrt",  sqrt,  "sqrt(x) - square root of x")
    MUP_UNARY_FUNC(FunSqrt1pm1, "sqrt1pm1", boost::math::sqrt1pm1,
                   "sqrt1pm1(x) - sqrt(1+x)-1")
    MUP_UNARY_FUNC(FunCbrt,"cbrt",boost::math::cbrt,
                   "cbrt(x) - cube root of x")

    MUP_UNARY_FUNC(FunExp,   "exp",   exp,   "exp(x) - e to the power of x")
    MUP_UNARY_FUNC(FunExpm1,"expm1",boost::math::expm1,"expm1(x) - exp(x)-1")
    MUP_UNARY_FUNC(FunAbs,   "abs",   fabs,  "abs(x) - absolute value of x")
    MUP_UNARY_FUNC(FunStep,  "step",  step,  "Heaviside step function")
    MUP_UNARY_FUNC(FunFloor,"floor",floor,
                   "floor(x) - largest integer not greater than x")
    MUP_UNARY_FUNC(FunCeil,"ceil",ceil,
                   "ceil(x) - smallest integer not less than x")

#undef MUP_UNARY_FUNC

#define MUP_BINARY_FUNC(CLASS, IDENT, FUNC, DESC)                     \
    CLASS::CLASS()                                                   \
    :ICallback(cmFUNC, _T(IDENT), 2)                                 \
    {}                                                               \
                                                                     \
    void CLASS::Eval(ptr_val_type &ret, const ptr_val_type *a_pArg, int)        \
    {                                                                \
      double a = a_pArg[0]->GetFloat();                               \
      double b = a_pArg[1]->GetFloat();                               \
      *ret = FUNC(a,b);                                              \
    }                                                                \
                                                                     \
    const char_type* CLASS::GetDesc() const                          \
    {                                                                \
      return _T(DESC);                                               \
    }                                                                \
                                                                     \
    IToken* CLASS::Clone() const                                     \
    {                                                                \
      return new CLASS(*this);                                       \
    }

    MUP_BINARY_FUNC(FunHypot,"hypot",boost::math::hypot,
                    "hypot(x,y) - sqrt(x^2 + y^2)")

#undef MUP_BINARY_FUNC

MUP_NAMESPACE_END
