#ifndef Extrapolater_H_
#define Extrapolater_H_

#include <AMReX_iMultiFab.H>
#include <AMReX_Geometry.H>
#include "CellFabArray.h"

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace amrex;

namespace Extrapolater
{
    // finebnd: boundary cells covered by fine cells (including periodically shifted fine cells)
    // crsebnd: boundary cells not covered by fine cells
    // physbnd: boundary cells outside the domain (excluding periodic boundaries)
    // interior: interior cells
    const int finebnd = 1;
    const int crsebnd = 0;
    const int physbnd = 0;
    const int interior = 1;

    //! It is expected that FillBoundary (w/ periodicity) has been called on mf.
    void FirstOrderExtrap (CellFabArray& mf, const Geometry& geom, int scomp, int ncomp)
    {
        Gpu::LaunchSafeGuard lsg(false); // xxxxx TODO gpu

		BL_ASSERT(mf.nGrow() == 1);
		BL_ASSERT(scomp >= 0);
		BL_ASSERT(ncomp <= mf.nComp());

		iMultiFab mask(mf.boxArray(), mf.DistributionMap(), 1, 1, MFInfo(),
					   DefaultFabFactory<IArrayBox>());
		mask.BuildMask(geom.Domain(), geom.periodicity(),
				       finebnd, crsebnd, physbnd, interior);

		int N = mf.nComp();

#ifdef _OPENMP
#pragma omp parallel
#endif
		for (MFIter mfi(mf); mfi.isValid(); ++mfi)
		{
			const Box& bx = mfi.validbox();
			const IArrayBox& maskfab = mask[mfi];
			const Box& maskbox = maskfab.box();
			CellFab& datafab = mf[mfi];
			const Box& databox = datafab.box();

			//amrex_first_order_extrap()
		}
    }
}

#endif
