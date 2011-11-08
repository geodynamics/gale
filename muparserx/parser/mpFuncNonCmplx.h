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
#ifndef MUP_FUNC_NON_CMPLX_H
#define MUP_FUNC_NON_CMPLX_H

#include "mpICallback.h"

/** \defgroup functions Function callback objects.

  This group lists the objects representing parser functions.
*/


MUP_NAMESPACE_START

#define MUP_UNARY_FUNC_DEF(CLASS)                                          \
    class CLASS : public ICallback                                         \
    {                                                                      \
    public:                                                                \
      CLASS();                                                             \
      virtual void Eval(ptr_val_type &ret, const ptr_val_type *a_pArg, int a_iArgc);  \
      virtual const char_type* GetDesc() const;                            \
      virtual IToken* Clone() const;                                       \
    }; 

    MUP_UNARY_FUNC_DEF(FunSin)
    MUP_UNARY_FUNC_DEF(FunCos)
    MUP_UNARY_FUNC_DEF(FunTan)
    // arcus functions
    MUP_UNARY_FUNC_DEF(FunASin)
    MUP_UNARY_FUNC_DEF(FunACos)
    MUP_UNARY_FUNC_DEF(FunATan)
    // hyperbolic functions
    MUP_UNARY_FUNC_DEF(FunSinH)
    MUP_UNARY_FUNC_DEF(FunCosH)
    MUP_UNARY_FUNC_DEF(FunTanH)
    // hyperbolic arcus functions
    MUP_UNARY_FUNC_DEF(FunASinH)
    MUP_UNARY_FUNC_DEF(FunACosH)
    MUP_UNARY_FUNC_DEF(FunATanH)
    // logarithm functions
    MUP_UNARY_FUNC_DEF(FunLog10)
    MUP_UNARY_FUNC_DEF(FunLog2)
    MUP_UNARY_FUNC_DEF(FunLog)
    MUP_UNARY_FUNC_DEF(FunLog1p)
    // error functions
    MUP_UNARY_FUNC_DEF(FunErf)
    MUP_UNARY_FUNC_DEF(FunErfc)
    // square root
    MUP_UNARY_FUNC_DEF(FunSqrt)
    MUP_UNARY_FUNC_DEF(FunSqrt1pm1)
    MUP_UNARY_FUNC_DEF(FunCbrt)

    MUP_UNARY_FUNC_DEF(FunExp)
    MUP_UNARY_FUNC_DEF(FunExpm1)
    MUP_UNARY_FUNC_DEF(FunAbs)
    MUP_UNARY_FUNC_DEF(FunStep)
    MUP_UNARY_FUNC_DEF(FunFloor)
    MUP_UNARY_FUNC_DEF(FunCeil)

    // This is really a binary function, but the signature is the
    // same.
    MUP_UNARY_FUNC_DEF(FunHypot)

#undef MUP_UNARY_FUNC_DEF

}  // namespace mu

#endif
