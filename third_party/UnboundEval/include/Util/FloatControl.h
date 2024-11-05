
#pragma once
#include <float.h>

#ifndef TARGET_CPU_ARM64
#include <xmmintrin.h>
#endif
#ifdef __clang__
#include <fenv.h>
#endif

namespace Util
{
	class SetAndRestoreFloatControlDownward
	{
	public:
		SetAndRestoreFloatControlDownward()
		{
            
			#if defined(_WIN32)
			// Save current rounding mode.
			m_previousMmCsr = _MM_GET_ROUNDING_MODE();
			m_previousControl87 = _control87(0, 0);

			// Set Denormal mode to "Flush to zero"
			_MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);
		    _control87( _DN_FLUSH, _MCW_DN	);

			// Set rounding mode to down.
			_MM_SET_ROUNDING_MODE(_MM_ROUND_DOWN);
			_control87(_RC_DOWN, _MCW_RC);
			#if !defined(_M_X64)
			_control87(_PC_53, _MCW_PC);
			#endif
			#endif

            #ifdef __clang__
                fegetenv(&m_prevFenv);
                fesetround( FE_DOWNWARD );
                feclearexcept( FE_ALL_EXCEPT  );
             
            #ifndef TARGET_CPU_ARM64
                m_previousMmCsr = _MM_GET_ROUNDING_MODE();
			    _MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);
		        _MM_SET_ROUNDING_MODE(_MM_ROUND_DOWN);
            #endif
			#endif
		}

		~SetAndRestoreFloatControlDownward()
		{
			// Restore previous rounding mode.
			#if defined(_WIN32)
			_MM_SET_ROUNDING_MODE(m_previousMmCsr);
			_control87(m_previousControl87, _MCW_RC);
			#if !defined(_M_X64)
			_control87(m_previousControl87, _MCW_PC);
			#endif
			#endif

            
            #ifdef __clang__
            #ifndef TARGET_CPU_ARM64
            _MM_SET_ROUNDING_MODE(m_previousMmCsr);
            #endif
            fesetenv(&m_prevFenv);
            #endif
		}

	private:
        #ifdef __clang__
            fenv_t m_prevFenv;
        #endif
		unsigned int m_previousMmCsr;
		unsigned int m_previousControl87;
	};
}
