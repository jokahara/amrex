#ifndef MyInterpolater_H_
#define MyInterpolater_H_

#include <AMReX_Box.H>
#include <AMReX_REAL.H>
#include <AMReX_Geometry.H>
#include <AMReX_BaseFab.H>
#include <AMReX_Interpolater.H>
#include "Cell.h"

using namespace amrex;
using CellFab = BaseFab< std::shared_ptr<Cell> >;

/**
* Specifies interpolater interface for coarse-to-fine interpolation in space.
* See AMReX_Interpolater.H
*/

class MyInterpolater : public Interpolater
{
public:

    MyInterpolater() {}

    ~MyInterpolater () {};

    /**
    * Returns coarsened box given fine box and refinement ratio.
    */
    inline Box CoarseBox (const Box& fine, int ratio) 
    {
        return CoarseBox(fine, ratio*IntVect::TheUnitVector());
    }

    /**
    * Returns coarsened box given fine box and refinement ratio.
    */
    inline Box CoarseBox (const Box& fine, const IntVect& ratio)
    {
        Box coarse = amrex::coarsen(fine, ratio);
        return coarse;
    }

    // Interpolates from coarse level to fine.
    // fine_region specifies the region to fill
    void interp (const CellFab&   crse,
                 int              crse_comp,
                 CellFab&         fine,
                 int              fine_comp,
                 int              ncomp,
                 const Box&       fine_region,
                 const IntVect&   ratio,
                 const Geometry&  crse_geom,
                 const Geometry&  fine_geom,
                 RunOn            runon) 
    {
        auto const crse_arr = crse.const_array();
        auto fine_arr = fine.array();
        
        Box fine_domain = fine_geom.Domain();
        Box crse_domain = crse_geom.Domain();
        Box coarse_region = CoarseBox(fine_region, ratio);

        // clear receiving cells
        ParallelFor(fine_region, ncomp, [&] (int i, int j, int k, int n) 
        {
            fine_arr(i,j,k,n)->clear();
        });

        // moving particles to fine cells

        //#pragma omp parallel
        ParallelFor(coarse_region, ncomp, [&] (int i, int j, int k, int n) 
        {
            IntVect iv{i,j,k};
            Cell* parent = crse_arr(iv, n).get();
            
            for (uint p = 0; p < parent->number_of_particles; p++)
            {
                auto& particle = parent->particles[p];
                // get the new location corresponding to fine geometry
                // (does not check for periodicity)
                IntVect new_iv = fine_geom.CellIndex(particle.data());
                
                if (fine_region.contains(new_iv)) 
                {
                    Cell* child = fine_arr(new_iv, n).get();
                    // add particles to new locations
                    child->particles.push_back(particle);
                    child->number_of_particles++;
                    continue;
                }
            }
        }); 
    }


    // Averages from fine level to coarse.
    void average(const CellFab&   crse,
                 int              crse_comp,
                 CellFab&         fine,
                 int              fine_comp,
                 int              ncomp,
                 const Box&       fine_region,
                 const IntVect&   ratio,
                 const Geometry&  crse_geom,
                 const Geometry&  fine_geom,
                 RunOn            runon) 
    {
        auto const crse_arr = crse.const_array();
        auto fine_arr = fine.array();
        
        Box coarse_region = CoarseBox(fine_region, ratio);
        
        ParallelFor(coarse_region, ncomp, [&] (int i, int j, int k, int n) 
        {
            IntVect cvect{i,j,k};
            Cell* parent = crse_arr(i,j,k,n).get();
            parent->clear();
            
            // make refined unit box
            Box fbox(cvect, cvect);
            fbox.refine(ratio);
            
            // copy particles from children to parent cell 
            For(fbox, [&] (int ii, int jj, int kk) 
            {
                Cell* child = fine_arr(ii,jj,kk,n).get();
                parent->particles.insert(parent->particles.begin(), 
                            child->particles.begin(), child->particles.end());
            });

            parent->number_of_particles = parent->particles.size();
        }); 
    }

    // This will not be used
    virtual void interp (const FArrayBox& crse,
                         int              crse_comp,
                         FArrayBox&       fine,
                         int              fine_comp,
                         int              ncomp,
                         const Box&       fine_region,
                         const IntVect&   ratio,
                         const Geometry&  crse_geom,
                         const Geometry&  fine_geom,
                         Vector<BCRec> const & bcr,
                         int              actual_comp,
                         int              actual_state,
                         RunOn            gpu_or_cpu) override {};
};

#endif