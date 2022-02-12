#ifndef GLOBAL_H
#define GLOBAL_H

#include "util.hh"
#include "gf.hh"
#include "extension.hh"

namespace global
{
    /* these are defined in main.cc */
    extern util::rand64bit randgen;
    extern GF2n F;
    extern Extension E;
}


#endif
