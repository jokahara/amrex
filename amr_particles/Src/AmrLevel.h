
#ifndef AmrLevel_H_
#define AmrLevel_H_

#include <AMReX_REAL.H>
#include <AMReX_Geometry.H>
#include <AMReX_LayoutData.H>
#include <AMReX_DistributionMapping.H>

#include "Amr.h"
#include "CellFabArray.h"
#include "CopyDescriptor.h"

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
    friend class FillPatchIterator;
    friend class FillPatchIteratorHelper;

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
              const DistributionMapping& dm,
              Real                       time=0);

    AmrLevel (const AmrLevel&) = delete;
    AmrLevel& operator = (const AmrLevel&) = delete;

    //! The destructor.
    ~AmrLevel ();

    //! Do an integration step on this level. Returns maximum safe time step.
    //virtual Real advance (Real time, Real dt, int  iteration, int  ncycle) { return 0; }

    //! Operations to be done after regridding
    void post_regrid (int lbase, int new_finest) { }
    
    //! Returns this AmrLevel.
    int Level () const noexcept { return level; }
    //! List of grids at this level.
    const BoxArray& boxArray () const noexcept { return grids; }
    const BoxArray& getEdgeBoxArray (int dir) const noexcept;
    const BoxArray& getNodalBoxArray () const noexcept;
    //
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
    Long countCells () const noexcept;

    //! Get the area not to tag.
    const BoxArray& getAreaNotToTag() noexcept;
    //! Get the area to tag.
    const Box& getAreaToTag() noexcept;
    //! Constuct the area not to tag.
    void constructAreaNotToTag();
    //! Set the area not to tag.
    void setAreaNotToTag(BoxArray& ba) noexcept;

    //! Error estimation for regridding.
    void errorEst (TagBoxArray& tb,
                   int          clearval,
                   int          tagval,
                   Real         time,
                   int          n_error_buf = 0,
                   int          ngrow = 0) { }

    //! Interpolate from coarse level to the valid area in dest.
    void FillCoarsePatch (CellFabArray& dest,
                          Real      time,
                          int       icomp,
                          int       ncomp,
			              int       nghost = 0);
    //! Function to set physical boundary conditions.
    /*void setPhysBoundaryValues (CellFab& dest, Real     time,
                                int      dest_comp,
                                int      src_comp,
                                int      num_comp);*/
                                        
    //! Data container.
    CellFabArray& getData (Real time=0) noexcept { return state; }

    //! Boundary condition access function.
    /*Vector<int> getBCArray (int gridno,
                            int icomp,
                            int ncomp);*/
    //! Hack to allow override of (non-fine-fine) fillpatched boundary data
    /*void set_preferred_boundary_values (CellFabArray& S,
                                        int       icomp,
                                        int       dcomp,
                                        int       ncomp,
                                        Real      time) const;*/
    /** 
    * \brief called in grid_places after other tagging routines to modify
    * the list of tagged points.  Default implementation does nothing.
    */
    void manual_tags_placement (TagBoxArray&           tags,
                                const Vector<IntVect>& bf_lev) {}
    /**
    * \brief Estimate the amount of work required to advance Just this level
    * based on the number of cells.
    * This estimate can be overwritten with different methods
    */
    virtual Real estimateWork();

    static void FillPatch (AmrLevel& AmrLevel,
                           CellFabArray& leveldata,
                           int       boxGrow,
                           Real      time,
                           int       icomp,
                           int       ncomp);

    static void FillPatchAdd (AmrLevel& AmrLevel,
                              CellFabArray& leveldata,
                              int       boxGrow,
                              Real      time,
                              int       icomp,
                              int       ncomp);

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
    CellFabArray          state;        // state data.

    BoxArray              m_AreaNotToTag; //Area which shouldn't be tagged on this level.
    Box                   m_AreaToTag;    //Area which is allowed to be tagged on this level.

    std::unique_ptr<FabFactory<CellFab> > m_factory;

private:

    mutable BoxArray      edge_grids[AMREX_SPACEDIM];  // face-centered grids
    mutable BoxArray      nodal_grids;              // all nodal grids
};

//
// Forward declaration.
//
class FillPatchIteratorHelper;

class FillPatchIterator : public MFIter
{
  public:

    friend class AmrLevel;

    FillPatchIterator (AmrLevel&  AmrLevel,
                       CellFabArray& leveldata);

    FillPatchIterator (AmrLevel& AmrLevel,
                       CellFabArray& leveldata,
                       int       boxGrow,
                       Real      time,
                       int       icomp,
                       int       ncomp);

    void Initialize (int  boxGrow,
                     Real time,
                     int  icomp,
                     int  ncomp);

    ~FillPatchIterator () {};

    CellFab& operator() () noexcept { return m_fabs[MFIter::index()]; }

    Box UngrownBox () const noexcept { return MFIter::validbox(); }

    CellFabArray& get_mf() noexcept { return m_fabs; }
    
  private:
    //
    // Disallowed.
    //
    FillPatchIterator ();
    FillPatchIterator (const FillPatchIterator& rhs);
    FillPatchIterator& operator= (const FillPatchIterator& rhs);

    void FillFromLevel0 (Real time, int icomp, int ncomp);
    void FillFromTwoLevels (Real time, int icomp, int ncomp);

    //
    // The data.
    //
    AmrLevel&                         m_Amrlevel;
    CellFabArray&                     m_leveldata;
    std::vector< std::pair<int,int> > m_range;
    CellFabArray                      m_fabs;
    int                               m_ncomp;
};

class FillPatchIteratorHelper
{
public:

    friend class FillPatchIterator;

    FillPatchIteratorHelper (AmrLevel&     AmrLevel,
                             CellFabArray& leveldata);

    FillPatchIteratorHelper (AmrLevel&     AmrLevel,
                             CellFabArray& leveldata,
                             int           boxGrow,
                             Real          time,
                             int           icomp,
                             int           ncomp);

    void Initialize (int           boxGrow,
                     Real          time,
                     int           icomp,
                     int           ncomp);

    ~FillPatchIteratorHelper () {};

    void fill (CellFab& fab, int idx);

private:
    //
    // Disallowed.
    //
    FillPatchIteratorHelper ();
    FillPatchIteratorHelper (const FillPatchIteratorHelper& rhs);
    FillPatchIteratorHelper& operator= (const FillPatchIteratorHelper& rhs);
    //
    // The data.
    //
    AmrLevel&                   m_Amrlevel;
    CellFabArray&               m_leveldata;
    CopyDescriptor              m_mfcd;
    Vector<Vector<FabArrayId> > m_mfid;
    std::map<int,Box>           m_ba;
    Real                        m_time;
    int                         m_growsize;
    int                         m_icomp;
    int                         m_ncomp;

    std::map< int,Vector< Vector<Box> > >                 m_fbox; // [grid][level][validregion]
    std::map< int,Vector< Vector<Box> > >                 m_cbox; // [grid][level][fillablesubbox]
    std::map< int,Vector< Vector< Vector<FillBoxId> > > > m_fbid; // [grid][level][fillablesubbox][oldnew]
};

#endif /*_AmrLevel_H_*/
