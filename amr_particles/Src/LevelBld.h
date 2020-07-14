
#ifndef LEVELBLD_H_
#define LEVELBLD_H_

#include "AmrLevel.h"

using namespace amrex;

class LevelBld
{
    static LevelBld bld;
    
    LevelBld* getLevelBld () { return &bld; }

    void variableSetUp () { AmrLevel::variableSetUp(); }
    void variableCleanUp () { AmrLevel::variableCleanUp(); }

    AmrLevel *operator() () { return new AmrLevel; }
    AmrLevel *operator() (Amr&            papa,
                          int             lev,
                          const Geometry& level_geom,
                          const BoxArray& ba,
                          const DistributionMapping& dm,
                          Real            time)
    { return new AmrLevel(papa, lev, level_geom, ba, dm, time); }
};

extern AMREX_ATTRIBUTE_WEAK LevelBld* getLevelBld ();

#endif