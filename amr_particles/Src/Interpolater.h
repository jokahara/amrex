#ifndef FillPatchUtil_I_H_
#define FillPatchUtil_I_H_

#include <AMReX_Box.H>
#include <AMReX_BCRec.H>
#include <AMReX_REAL.H>
#include <AMReX_GpuControl.H>
#include <AMReX_Geometry.H>

using namespace amrex;

bool ProperlyNested (const IntVect& ratio, const IntVect& blocking_factor, int ngrow,
                     const IndexType& boxType)
{
    int ratio_max = ratio[0];
#if (AMREX_SPACEDIM > 1)
    ratio_max = std::max(ratio_max, ratio[1]);
#endif
#if (AMREX_SPACEDIM == 3)
    ratio_max = std::max(ratio_max, ratio[2]);
#endif
    // There are at least this many coarse cells outside fine grids
    // (except at physical boundaries).
    const IntVect& nbuf = blocking_factor / ratio_max;

    Box crse_box(IntVect(AMREX_D_DECL(0 ,0 ,0 )), IntVect(AMREX_D_DECL(4*nbuf[0]-1,
                                                                       4*nbuf[1]-1,
                                                                       4*nbuf[2]-1)));
    crse_box.convert(boxType);
    Box fine_box(nbuf, IntVect(AMREX_D_DECL(3*nbuf[0]-1,3*nbuf[1]-1,3*nbuf[2]-1)));
    fine_box.convert(boxType);
    fine_box.refine(ratio_max);
    fine_box.grow(ngrow);

    Box& fine_box_coarsened = amrex::coarsen(fine_box, ratio_max).grow(1);

    return crse_box.contains(fine_box_coarsened);
}

#endif