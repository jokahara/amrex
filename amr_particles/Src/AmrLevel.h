
#ifndef AmrLevel_H_
#define AmrLevel_H_

#include <AMReX_REAL.H>
#include <AMReX_Geometry.H>
#include <AMReX_LayoutData.H>
#include <AMReX_DistributionMapping.H>

#include "Amr.h"
#include "CellFabArray.h"

#include <vector>
#include <memory>
#include <map>

using namespace amrex;

class Amr;

/**
* \brief base class for managing individual levels.
* AmrLevel functions both as a container for state data on a level
* and also manages the advancement of data in time.
*/
class AmrLevel
{
    friend class AmrIter;

public:
    //
    //Default constructor.  Builds invalid object.
    //
    AmrLevel () noexcept;
    //
    //The basic constructor.
    //
    AmrLevel (Amr&                       papa,
              int                        lev,
              const Geometry&            level_geom,
              const BoxArray&            bl,
              const DistributionMapping& dm);

    AmrLevel (const AmrLevel&) = delete;
    AmrLevel& operator = (const AmrLevel&) = delete;

    //! The destructor.
    ~AmrLevel ();

    /**
    * Init data on this level after regridding if old AmrLevel
    * AmrLevel* old is given, level is filled with it.
    * Must be called initializing!
    */
    void init (AmrLevel* old = nullptr);

    //! Operations to be done after regridding
    void post_regrid (int lbase, int new_finest) {};
    
    //! Returns this AmrLevel.
    int Level () const noexcept { return level; }
    //! List of propagated grids at this level.
    const BoxArray& boxArray () const noexcept { return grids; }
    //! Distribution of propagated grids among processes
    const DistributionMapping& DistributionMap () const noexcept { return dmap; }
    //
    const FabFactory<CellFab>& Factory () const noexcept { return *m_factory; }
    //! Number of grids at this level.
    int numGrids () const noexcept { return grids.size(); }
    //! Returns the indices defining physical domain.
    const Box& Domain () const noexcept { return geom.Domain(); }
    //! Returns the geometry object.
    const Geometry& Geom () const noexcept { return geom; }
    //! Refinement ratio to finer level.
    const IntVect& fineRatio () const noexcept { return fine_ratio; }
    //! Returns number of cells on level.
    Long countCells () const noexcept { return grids.numPts(); }

    // Builds coarse and fine boundaries, and sets which areas will be propagated
    void constructCrseFineBdry(AmrLevel* fine);

    //! Get the area not to tag.
    const BoxArray& getAreaNotToTag() noexcept { return m_AreaNotToTag; }
    //! Get the area to tag.
    const Box& getAreaToTag() noexcept { return m_AreaToTag; }
    //! Constuct the area not to tag.
    void constructAreaNotToTag();
    //! Set the area not to tag.
    void setAreaNotToTag(BoxArray& ba) noexcept { m_AreaNotToTag = ba; }
    
    //! Error estimation for regridding.
    void errorEst (TagBoxArray& tb,
                   int          clearval,
                   int          tagval,
                   Real         time,
                   int          n_error_buf = 0,
                   int          ngrow = 0);

    //! Interpolate from coarse level to the valid area in dest.
    void FillFromCoarsePatch (CellFabArray& dest,
                              int       icomp,
                              int       ncomp,
                              int       nghost = 0);
                                        
    //! Data container.
    CellFabArray& getCells () noexcept 
    { return *state; }

    CellFab& operator[](MFIter& mfi) noexcept 
    { return (*state)[mfi]; }

    /** 
    * Called in grid_places after other tagging routines to modify
    * the list of tagged points.  Default implementation does nothing.
    */
    void manual_tags_placement (TagBoxArray&           tags,
                                const Vector<IntVect>& bf_lev) {}

    // Estimate the amount of work required to advance just this level
    void estimateWork(MultiFab& mf) { state->EstimateWork(mf); }

    void FillPatch (AmrLevel& AmrLevel,
                    int       boxGrow,
                    int       icomp,
                    int       ncomp);

    /*void FillPatchAdd (AmrLevel& AmrLevel,
                        int       boxGrow,
                        Real      time,
                        int       icomp,
                        int       ncomp);*/
    
    // container for overlapping coarse and fine data. 
    struct BoundaryContainer
    {
        std::unique_ptr<CellFabArray> coarse; // coarse data
        std::unique_ptr<CellFabArray> fine;   // fine data
        AmrLevel* src;

        BoundaryContainer() {};

        // reset boundaries and connect cells to their parents
        void reset (CellFabArray* new_coarse = nullptr, AmrLevel* coarse_parent = nullptr, 
                    CellFabArray* new_fine = nullptr, AmrLevel* fine_parent = nullptr,
                    AmrLevel* source_level = nullptr) noexcept;

    private:
        void connect_cells(CellFabArray& data, AmrLevel* parent);
    };

    BoundaryContainer& getFineBoudary() noexcept { return fine_boundary; };
    BoundaryContainer& getCoarseBoudary() noexcept { return coarse_boundary; };

    // fill fine cells from coarse data
    void FillCoarseToFine();
    // fill coarse cells from fine data
    void FillFineToCoarse();

protected:
    //
    // The Data.
    //
    int                   level;        // AMR level (0 is coarsest).
    Geometry              geom;         // Geom at this level.
    BoxArray              grids;        // Cell-centered locations of grids.
    DistributionMapping   dmap;         // Distribution of grids among processes
    Amr*                  parent;       // Pointer to parent AMR structure.
    IntVect               crse_ratio;   // Refinement ratio to coarser level.
    IntVect               fine_ratio;   // Refinement ratio to finer level.
    
    // All cells
    std::unique_ptr<CellFabArray> state;
    // std::unique_ptr<CellFabArray> new_state;
    
    // local Cellfabs that overlap with finer levels are not propagated.
    std::vector<int> propagate_fab;

    // Ghost cells on the coarse/fine boundaries;
    BoundaryContainer coarse_boundary;  // boundary with coarser level
    BoundaryContainer fine_boundary;    // boundary with finer level

    BoxArray              m_AreaNotToTag; //Area which shouldn't be tagged on this level.
    Box                   m_AreaToTag;    //Area which is allowed to be tagged on this level.
    
    std::unique_ptr<FabFactory<CellFab> > m_factory;
};

class AmrIter : public MFIter 
{
public:
    friend class AmrLevel;

    AmrIter (AmrLevel& AmrLevel)
    : MFIter(*AmrLevel.state)
    {
        propagate = AmrLevel.propagate_fab.data();
    }

    inline void operator++() noexcept 
    {
        propagate++;
        ++currentIndex;
    }

    inline bool isPropagated() noexcept
    { return *propagate; }
    
    ~AmrIter () {};
    
private:
    // Disallowed.
    AmrIter ();
    AmrIter (const AmrIter& rhs);
    AmrIter& operator= (const AmrIter& rhs);

    int* propagate;
};

#endif /*_AmrLevel_H_*/
