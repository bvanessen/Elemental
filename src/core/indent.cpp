/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/

#include <string>

#include "El/core/environment/decl.hpp"
#include "El/Typedefs.hpp"

namespace
{

El::Int indentLevel=0;
El::Int spacesPerIndent=2;

}

namespace El
{

Int PushIndent() { return ::indentLevel++; }
Int PopIndent() { return ::indentLevel--; }
void SetIndent( Int indent ) { ::indentLevel = indent; }
void ClearIndent() { ::indentLevel = 0; }
Int IndentLevel() { return ::indentLevel; }

std::string Indent()
{
    std::string ind;
    for( Int i=0; i < ::spacesPerIndent * ::indentLevel; ++i )
        ind = ind + " ";
    return ind;
}

} // namespace El
