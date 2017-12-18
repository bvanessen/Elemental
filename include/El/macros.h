#ifndef EL_MACROS_H_
#define EL_MACROS_H_

#include "El/config.h"

// MACROS TO CONTROL RELEASE/DEBUG-SPECIFIC COMPILATION
// FIXME (trb 12/16/17): Not sure how I feel about these wrapping huge
// code blocks.
#ifdef EL_RELEASE
# define EL_DEBUG_ONLY(cmd)
# define EL_RELEASE_ONLY(cmd) cmd;
#else
# define EL_DEBUG_ONLY(cmd) cmd;
# define EL_RELEASE_ONLY(cmd)
#endif

// NOEXCEPT SPECIFIERS.
#ifdef EL_HAVE_NO_EXCEPT
# define EL_NO_EXCEPT noexcept
#else
# define EL_NO_EXCEPT
#endif

#ifdef EL_RELEASE
# define EL_NO_RELEASE_EXCEPT EL_NO_EXCEPT
#else
# define EL_NO_RELEASE_EXCEPT
#endif

#endif // EL_MACROS_H_
