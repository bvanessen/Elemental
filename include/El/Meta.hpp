#ifndef EL_META_HPP_
#define EL_META_HPP_

#include <type_traits>

namespace El
{

template<typename Condition,class T=void>
using EnableIf = typename std::enable_if<Condition::value,T>::type;
template<typename Condition,class T=void>
using DisableIf = typename std::enable_if<!Condition::value,T>::type;

}//namespace El
#endif // EL_META_HPP_
