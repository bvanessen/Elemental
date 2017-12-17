/*
  Copyright (c) 2009-2016, Jack Poulson
  All rights reserved.

  This file is part of Elemental and is under the BSD 2-Clause License,
  which can be found in the LICENSE file in the root directory, or at
  http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_ENVIRONMENT_DECL_HPP
#define EL_ENVIRONMENT_DECL_HPP

// FIXME (trb 12/15/17): THIS FILE IS TERRRRRRRIBLE AND, FRANKLY, A
// GOOD EXAMPLE OF SOMETHING THAT LEADS TO ALL THOSE REASONS WHY
// PEOPLE RELIGIOUSLY HATE C++

#include <exception>
#include <functional>
#include <iostream>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <vector>

#include "El/core/imports/mpi.hpp"
#include "El/core/imports/mpi_choice.hpp"

#include "El/macros.h"

namespace El
{

void PrintVersion( std::ostream& os=std::cout );
void PrintConfig( std::ostream& os=std::cout );
void PrintCCompilerInfo( std::ostream& os=std::cout );
void PrintCxxCompilerInfo( std::ostream& os=std::cout );
bool Using64BitInt();
bool Using64BitBlasInt();

// For manually initializing and finalizing Elemental; their direct usage
// in C++ programs now deprecated.
void Initialize();
void Initialize( int& argc, char**& argv );
void Finalize();
bool Initialized();

// For initializing/finalizing Elemental using RAII
class Environment
{
public:
    Environment() { Initialize(); }
    Environment( int& argc, char**& argv ) { Initialize( argc, argv ); }
    ~Environment() { Finalize(); }
};

// For getting the MPI argument instance (for internal usage)
class Args : public choice::MpiArgs
{
public:
    Args
    ( int argc, char** argv,
      mpi::Comm comm=mpi::COMM_WORLD, std::ostream& error=std::cerr )
        : choice::MpiArgs(argc,argv,comm,error)
    { }
    virtual ~Args() { }
protected:
    virtual void HandleVersion( std::ostream& os=std::cout ) const;
    virtual void HandleBuild( std::ostream& os=std::cout ) const;
};
Args& GetArgs();

// For processing command-line arguments
template<typename T>
T Input( std::string name, std::string desc );
template<typename T>
T Input( std::string name, std::string desc, T defaultVal );
void ProcessInput();
void PrintInputReport();

// For getting and setting the algorithmic blocksize
Int Blocksize();
void SetBlocksize( Int blocksize );

// For manipulating the algorithmic blocksize as a stack
void PushBlocksizeStack( Int blocksize );
void PopBlocksizeStack();
void EmptyBlocksizeStack();

template<typename T,
         typename=EnableIf<IsScalar<T>>>
const T& Max( const T& m, const T& n ) EL_NO_EXCEPT;
inline const Int& Max( const Int& m, const Int& n ) EL_NO_EXCEPT;

template<typename T,
         typename=EnableIf<IsScalar<T>>>
const T& Min( const T& m, const T& n ) EL_NO_EXCEPT;
inline const Int& Min( const Int& m, const Int& n ) EL_NO_EXCEPT;

// Replacement for std::memcpy, which is known to often be suboptimal.
// Notice the sizeof(T) is no longer required.
template<typename T,
         typename=EnableIf<IsPacked<T>>>
void MemCopy
(       T* dest,
        const T* source,
        size_t numEntries );
template<typename T,
         typename=DisableIf<IsPacked<T>>,
         typename=void>
void MemCopy
(       T* dest,
        const T* source,
        size_t numEntries );

template<typename T,
         typename=EnableIf<IsPacked<T>>>
void MemSwap
( T* a,
  T* b,
  T* temp,
  size_t numEntries );
template<typename T,
         typename=DisableIf<IsPacked<T>>,
         typename=void>
void MemSwap
( T* a,
  T* b,
  T* temp,
  size_t numEntries );

// Generalization of std::memcpy so that unit strides are not required
template<typename T,
         typename=EnableIf<IsPacked<T>>>
void StridedMemCopy
(       T* dest,   Int destStride,
        const T* source, Int sourceStride, Int numEntries );
template<typename T,
         typename=DisableIf<IsPacked<T>>,
         typename=void>
void StridedMemCopy
(       T* dest,   Int destStride,
        const T* source, Int sourceStride, Int numEntries );

// A thin wrapper around std::copy
template<typename S,typename T>
void CopySTL( const S& a, T& b );

// Replacement for std::memset, which is likely suboptimal and hard to extend
// to non-POD datatypes. Notice that sizeof(T) is no longer required.
template<typename T,
         typename=EnableIf<IsPacked<T>>>
void MemZero( T* buffer, size_t numEntries );
template<typename T,
         typename=DisableIf<IsPacked<T>>,
         typename=void>
void MemZero( T* buffer, size_t numEntries );

// Clear the contents of x by swapping with an empty object of the same type
template<typename T>
void SwapClear( T& x );

// Reserve memory in a vector without zero-initializing the variables unless
// valgrind is currently running or the datatype *requires* construction.
template<typename T,
         typename=EnableIf<IsPacked<T>>>
void FastResize( std::vector<T>& v, Int numEntries );
template<typename T,
         typename=DisableIf<IsPacked<T>>,
         typename=void>
void FastResize( std::vector<T>& v, Int numEntries );

inline void BuildStream( std::ostringstream& ) { }

template<typename T,typename... ArgPack>
void BuildStream( std::ostringstream& os, const T& item,
                  const ArgPack& ... args );
template<typename... ArgPack>
std::string BuildString( const ArgPack& ... args );

class UnrecoverableException : public std::runtime_error
{
public:
    UnrecoverableException( const char* msg="Unrecoverable exception" )
        : std::runtime_error( msg ) { }
};

template<typename... ArgPack>
void UnrecoverableError( const ArgPack& ... args );
template<typename... ArgPack>
void LogicError( const ArgPack& ... args );
template<typename... ArgPack>
void RuntimeError( const ArgPack& ... args );

// This is the only place that Elemental is currently using duck-typing.
// I'm not sure if it's a good idea to use it more often.
template<class MatType>
std::string DimsString( const MatType& A, std::string label="Matrix" );

// This is defined in choice.hpp
class ArgException;

// An exception which signifies that a matrix was unexpectedly singular.
class SingularMatrixException : public std::runtime_error
{
public:
    SingularMatrixException( const char* msg="Matrix was singular" )
        : std::runtime_error( msg ) { }
};

// An exception which signifies a zero pivot was chosen, though the matrix
// may not actually be singular
class ZeroPivotException : public std::runtime_error
{
public:
    ZeroPivotException( const char* msg="Zero pivot was chosen" )
        : std::runtime_error( msg ) { }
};

// An exception which signifies that a matrix was unexpectedly non-HPD
class NonHPDMatrixException  : public std::runtime_error
{
public:
    NonHPDMatrixException( const char* msg="Matrix was not HPD" )
        : std::runtime_error( msg ) { }
};

// An exception which signifies that a matrix was unexpectedly non-HPSD
class NonHPSDMatrixException  : public std::runtime_error
{
public:
    NonHPSDMatrixException( const char* msg="Matrix was not HPSD" )
        : std::runtime_error( msg ) { }
};

EL_DEBUG_ONLY(
    void EnableTracing();
    void DisableTracing();

    void PushCallStack( std::string s );
    void PopCallStack();
    void DumpCallStack( std::ostream& os=std::cerr );

    class CallStackEntry
    {
    public:
        CallStackEntry( std::string s )
        {
            if( !std::uncaught_exception() )
                PushCallStack(s);
        }
        ~CallStackEntry()
        {
            if( !std::uncaught_exception() )
                PopCallStack();
        }
    };
    typedef CallStackEntry CSE;
    )

void OpenLog( const char* filename );

std::ostream & LogOS();

template<typename... ArgPack>
void Log( const ArgPack& ... args );

void CloseLog();

void ReportException( const std::exception& e, std::ostream& os=std::cerr );

void ComplainIfDebug();

Int PushIndent();
Int PopIndent();
void SetIndent( Int level );
void ClearIndent();
Int IndentLevel();
std::string Indent();

template<typename... ArgPack>
void Output( const ArgPack& ... args );

template<typename... ArgPack>
void OutputFromRoot( mpi::Comm comm, const ArgPack& ... args );

template<typename T>
void EnsureConsistent( T alpha, mpi::Comm comm, std::string name="" );

// This will be guaranteed by C++14 via std::make_unique
#if __cplusplus < 201402L
template<typename T,typename ...ArgPack>
std::unique_ptr<T> MakeUnique( ArgPack&& ...args )
{ return std::unique_ptr<T>( new T( std::forward<ArgPack>(args)... ) ); }
#else
template<typename T,typename ...ArgPack>
std::unique_ptr<T> MakeUnique( ArgPack&& ...args )
{ return std::make_unique<T>(std::forward<ArgPack>(args)...); }
#endif

// MakeFunction allows for the automatic conversion of a lambda to an
// std::function and is similar to http://stackoverflow.com/a/24068396/1119818.

/*
  namespace make_function {

  template<typename T>
  struct Helper { using type = void; };
  template<typename Ret,typename Class,typename... Args>
  struct Helper<Ret(Class::*)(Args...) const>
  { using type = std::function<Ret(Args...)>; };

  } // namespace make_function

  template<typename Function>
  typename make_function::Helper<decltype(&Function::operator())>::type
  MakeFunction(Function const& func) { return func; }
*/

// Handles generic types that are functors, delegate to its 'operator()'
template<typename T>
struct function_traits : public function_traits<decltype(&T::operator())>
{};

// Handles pointers to member functions
template<typename ClassType,typename ReturnType,typename... Args>
struct function_traits<ReturnType(ClassType::*)(Args...) const>
{
    enum { arity = sizeof...(Args) };
    typedef std::function<ReturnType (Args...)> f_type;
};

// Handles pointers to member functions
template<typename ClassType,typename ReturnType,typename... Args>
struct function_traits<ReturnType(ClassType::*)(Args...)>
{
    enum { arity = sizeof...(Args) };
    typedef std::function<ReturnType (Args...)> f_type;
};

// Handles function pointers
template<typename ReturnType,typename... Args>
struct function_traits<ReturnType (*)(Args...)>
{
    enum { arity = sizeof...(Args) };
    typedef std::function<ReturnType (Args...)> f_type;
};

template<typename L>
static typename function_traits<L>::f_type MakeFunction(L l)
{ return (typename function_traits<L>::f_type)(l); }

// Handles std::bind & multiple function call operator()'s
template<typename ReturnType,typename... Args,class T>
auto MakeFunction(T&& t)
    -> std::function<decltype(ReturnType(t(std::declval<Args>()...)))(Args...)>
{ return {std::forward<T>(t)}; }

// For explicit overloads
template<typename ReturnType,typename... Args>
auto MakeFunction(ReturnType(*p)(Args...))
    -> std::function<ReturnType(Args...)>
{ return {p}; }

// For explicit overloads
template<typename ReturnType,typename... Args,typename ClassType>
auto MakeFunction(ReturnType(ClassType::*p)(Args...))
    -> std::function<ReturnType(Args...)>
{ return {p}; }

template<typename T>
T Scan( const std::vector<T>& counts, std::vector<T>& offsets );

template<typename T>
bool IsSorted( const std::vector<T>& x );
// While is_strictly_sorted exists in Boost, it does not exist in the STL (yet)
template<typename T>
bool IsStrictlySorted( const std::vector<T>& x );

void Union
( std::vector<Int>& both,
  const std::vector<Int>& first, const std::vector<Int>& second );
std::vector<Int>
Union( const std::vector<Int>& first, const std::vector<Int>& second );

void RelativeIndices(std::vector<Int>& relInds, const std::vector<Int>& sub,
                     const std::vector<Int>& full);
std::vector<Int> RelativeIndices(
    const std::vector<Int>& sub, const std::vector<Int>& full );

// Insists that the index can be found
Int Find( const std::vector<Int>& sortedInds, Int index );

#ifdef EL_HAVE_PRETTY_FUNCTION
# define EL_FUNCTION __PRETTY_FUNCTION__
#else
# define EL_FUNCTION __func__
#endif

#define EL_LOGIC_ERROR(...)                                             \
    El::LogicError(EL_FUNCTION," in ",__FILE__,"@",__LINE__,": ",__VA_ARGS__);
#define EL_RUNTIME_ERROR(...)                                           \
    El::RuntimeError(EL_FUNCTION," in ",__FILE__,"@",__LINE__,": ",__VA_ARGS__);
#define EL_DEBUG_CSE EL_DEBUG_ONLY(El::CSE cse(EL_FUNCTION))

} // namespace El

#endif // ifndef EL_ENVIRONMENT_DECL_HPP
