#ifndef CONFIG_H_
#define CONFIG_H_

#include <chrono>

#ifdef LOG_OUTPUT

#define lprintf(...) printf(__VA_ARGS__)
#define lputs(...) puts(__VA_ARGS__)

#else

#define lprintf(...) void(0)
#define lputs(...) void(0)

#endif


namespace qflow {

// simulation of Windows GetTickCount()
    unsigned long long inline GetCurrentTime64() {
        using namespace std::chrono;
        return duration_cast<milliseconds>(steady_clock::now().time_since_epoch()).count();
    }
} // namespace qflow

#endif
