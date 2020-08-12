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

    /**
    * Coarse to fine interpolation in space.
    * This is a pure virtual function and hence MUST
    * be implemented by derived classes.
    */
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
        Box coarse_region = crse.box();
        
        // clear receiving cells
        #pragma omp parallel
        ParallelFor(fine_region, ncomp, [&] (int i, int j, int k, int n) 
        {
            fine_arr(i,j,k,n)->clear();

            // Alternative way of iterating:
            /*IntVect fv{i,j,k};
            Cell* child = fine_arr(fv, n).get();
            IntVect cv = amrex::coarsen(fv, ratio);
            Cell* parent = crse_arr(cv, n).get();*/
        });

        // moving particles to fine cells

                //Print() << coarse_region << "->" << fine_region << "\n";
        ParallelFor(coarse_region, ncomp, [&] (int i, int j, int k, int n) 
        {
            IntVect iv{i,j,k};
            Cell* parent = crse_arr(iv, n).get();
            
            for (uint p = 0; p < parent->number_of_particles; p++)
            {
                auto& particle = parent->particles[p];
                // get the new location corresponding to fine geometry
                IntVect new_iv = fine_geom.CellIndex(particle.data());

                if (fine_arr.contains(new_iv[0], new_iv[1], new_iv[2])) {
                    Cell* child = fine(new_iv, n).get();
                    // add particles to new locations
                    //Print() << "push: " << particle[0] << "\n";
                    child->particles.push_back(particle);
                    child->number_of_particles++;
                    continue;
                }

                // periodicity checks
                if (i < 0 || i >= crse_domain.length(0)) {
                    if (fine_geom.isPeriodic(0)) {
                        new_iv[0] += fine_domain.length(0);
                    }
                }
                if (j < 0 || j >= crse_domain.length(1)) {
                    if (fine_geom.isPeriodic(1)) {
                        new_iv[0] += fine_domain.length(1);
                    }
                }
                //Print() << new_iv << "\n";
                if (fine_arr.contains(new_iv[0], new_iv[1], new_iv[2])) {
                    Cell* child = fine(new_iv, n).get();
                    child->particles.push_back(particle);
                    child->number_of_particles++;
                    continue;
                }
            }
        }); 
    }


    /**
    * Fine to coarse interpolation in space.
    * This is a pure virtual function and hence MUST
    * be implemented by derived classes.
    */
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
        
        Box coarse_region = crse.box();

        //const int volume = ratio[0]*ratio[1]*ratio[2];
        
        ParallelFor(coarse_region, ncomp, [&] (int i, int j, int k, int n) 
        {
            IntVect cvect{i,j,k};
            Cell* parent = crse_arr(i,j,k,n).get();
            parent->clear();
            
            Box fbox(cvect, cvect);
            fbox.refine(ratio);
            // add particles from fine cells to parent cell 
            For(fbox, [&] (int ii, int jj, int kk) 
            {
                Cell* child = fine_arr(ii,jj,kk,n).get();
                parent->particles.insert(parent->particles.begin(), 
                            child->particles.begin(), child->particles.end());
            });

            parent->number_of_particles = parent->particles.size();
            //Print() << "added " << parent->number_of_particles << " to " << cvect << "\n";
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