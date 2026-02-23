#ifndef __THROW_WRAPPER
#define __THROW_WRAPPER

// Enable OpenMP debugging by defining VERBOSE_THROW
// #define VERBOSE_THROW

#ifdef VERBOSE_THROW

    #include <iostream>

    // Make sure that we capture error strings even in OpenMP blocks
    #define THROW(X) \
        std::cerr << X << std::endl;\
        throw X;
#else

    #define THROW(X) \
        throw X;
#endif // VERBOSE_THROW 

#endif // __THROW_WRAPPER