#ifndef EL_CORE_DISTMATRIX_DECL_HPP_
#define EL_CORE_DISTMATRIX_DECL_HPP_

#include "El/core/Grid.hpp"

namespace El
{

enum class DistWrap
{
    ELEMENT,
    BLOCK
};

template<typename T=double,
         Dist U=Dist::MC, Dist V=Dist::MR, DistWrap wrap=DistWrap::ELEMENT>
class DistMatrix;

}// namespace El
#endif // EL_CORE_DIST_MATRIX_DECL_HPP_
