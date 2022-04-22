/* Copyright 2022 Eetu Karppinen. Subject to the MIT license. */
#ifndef GLOBAL_H
#define GLOBAL_H

#include <random>

class GF2n;
class Extension;

namespace util
{
    class rand64bit
    {
    private:
        std::mt19937_64 gen;
        std::uniform_int_distribution<uint64_t> dist;
    public:
        rand64bit() {}
        void init(uint64_t seed) { gen = std::mt19937_64(seed); }
        uint64_t operator() () { return this->dist(gen); }
    };
}

namespace global
{
    /* these are defined in main.cc */
    extern util::rand64bit randgen;
    extern GF2n F;
    extern Extension E;
    extern bool output;
}


#endif
