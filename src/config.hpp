#ifndef CONFIG_H_
#define CONFIG_H_

#ifdef LOG_OUTPUT

#define lprintf(...) printf(__VA_ARGS__)
#define lputs(...) puts(__VA_ARGS__)

#else

#define lprintf(...) void(0)
#define lputs(...) void(0)

#endif


#endif
