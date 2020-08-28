
#ifndef AMREX_Amr_H_
#define AMREX_Amr_H_

#include <fstream>
#include <memory>
#include <list>

#include <AMReX_AmrCore.H>
#include <AMReX_Box.H>
#include <AMReX_Geometry.H>
#include <AMReX_BoxArray.H>
#include <AMReX_Array.H>
#include <AMReX_Vector.H>

#include "CellFabArray.h"
#include "AmrLevel.h"

/**
* \brief Manage hierarchy of levels for time-dependent AMR computations.
*
* The Amr class is designed to manage parts of the computation  which do
* not belong on a single level, like establishing and updating the hierarchy
* of levels, global timestepping, and managing the different AmrLevels
*/
using namespace amrex;

//TODO: namespace Amr {

class AmrLevel;

class Amr : public AmrCore
{

public:
    //! The constructor.
    Amr ();

    Amr (const RealBox* rb, int max_level_in, const Vector<int>& n_cell_in, int coord);

    Amr (const Amr& rhs) = delete;
    Amr& operator= (const Amr& rhs) = delete;

    static void Initialize ();
    static void Finalize ();

    void InitAmr ();
    
    //! Initialize coarsest level.
    void defBaseLevel (const BoxArray* lev0_grids = 0, const Vector<int>* pmap = 0);

    //! Build refined levels. Must be called after initBaseLevel.
    void bldFineLevels ();

    //! The destructor.
    virtual ~Amr ();

	//! Fills ghost cells with data from other levels and remote neighbours.
    void FillAllBoundaries();

    //! Rebuild grid hierarchy finer than lbase and load balance.
    void regrid (int  lbase, bool initial = false);

    //! AmrLevel lev.
    AmrLevel& getLevel (int lev) noexcept { return *amr_level[lev]; }
    AmrLevel& operator[] (int lev) { return *amr_level[lev]; }
    
    //! Number of components.
    int nComp () noexcept { return n_comp; }
    //! Number of ghost cells.
    IntVect nGrow () noexcept { return n_grow; }
    //! Total number of cells.
    Long cellCount () noexcept;
    //! Number of cells at given level.
    Long cellCount (int lev) noexcept;
    //! Total number of grids.
    int numGrids () noexcept;
    //! Number of grids at given level.
    int numGrids (int lev) noexcept;
    //! Should we regrid this level?
    bool okToRegrid (int level) noexcept;
    //! Array of BoxArrays read in to initially define grid hierarchy
    static const BoxArray& initialBa (int level) noexcept
        { BL_ASSERT(level-1 < initial_ba.size()); return initial_ba[level-1]; }
    //! Number of levels at which the grids are initially specified
    static int initialBaLevels () noexcept { return initial_ba.size(); }

    const Vector<BoxArray>& getInitialBA() noexcept;


    void printGridInfo (std::ostream& os,
                        int           min_lev,
                        int           max_lev);

protected:

    //! Initialize grid hierarchy -- called by Amr::init.
    void initialInit (const BoxArray* lev0_grids = 0, const Vector<int>* pmap = 0);
    //! Check for valid input.
    void checkInput ();
    //! Define new grid locations (called from regrid) and put into new_grids.
    void MakeGrids (int               lbase,
                    int&              new_finest,
                    Vector<BoxArray>& new_grids);

    // Makes distribution maps from given grids
    void MakeDistributionMaps (int new_finest, 
                               Vector<BoxArray>& new_grids,
                               Vector<DistributionMapping>& new_dmap);

    void InstallNewLevels (int new_finest, 
                           const Vector<BoxArray>& grids,
                           const Vector<DistributionMapping>& dmap);

    // Make a new level using provided BoxArray and DistributionMapping and
    // fill with interpolated coarse level data.
    // overrides the pure virtual function in AmrCore
    virtual void MakeNewLevelFromCoarse (int lev, Real time, const BoxArray& ba,
					 const DistributionMapping& dm) override
	{ amrex::Abort("How did we get her!"); }

    // Remake an existing level using provided BoxArray and DistributionMapping and
    // fill with existing fine and coarse data.
    // overrides the pure virtual function in AmrCore
    virtual void RemakeLevel (int lev, Real time, const BoxArray& ba,
			      const DistributionMapping& dm) override
	{ amrex::Abort("How did we get her!"); }

    // Delete level data
    // overrides the pure virtual function in AmrCore
    void ClearLevel (int lev) override;

    // Make a new level from scratch using provided BoxArray and DistributionMapping.
    // Only used during initialization.
    // overrides the pure virtual function in AmrCore
    virtual void MakeNewLevelFromScratch (int lev, Real time, const BoxArray& ba,
					  const DistributionMapping& dm) override
	{ amrex::Abort("How did we get her!"); }

    // tag all cells for refinement
    // overrides the pure virtual function in AmrCore
    void ErrorEst (int lev, TagBoxArray& tags, Real time, int ngrow);
    // contains cells we do not want to refine (for instance outer boundary)
    BoxArray GetAreaNotToTag (int lev);

    // Remove or add tagged points which violate/satisfy additional
    // user-specified criteria. (called after ErrorEst)
    void ManualTagsPlacement (int lev, TagBoxArray& tags, const Vector<IntVect>& bf_lev);

    //
    // The data ...
    //
    std::string     regrid_grids_file;   //!< Grids file that will bypass regridding.
    std::string     initial_grids_file;  //!< Grids file that will bypass regridding only at initialization.
    Vector<std::unique_ptr<AmrLevel> > amr_level;    //!< Vector of levels

    bool            isPeriodic[AMREX_SPACEDIM];  //!< Domain periodic?
    Vector<int>     regrid_int;      //!< Interval between regridding.

    bool            abort_on_stream_retry_failure;
    int             stream_max_tries;
    int             loadbalance_with_workestimates;
    int             loadbalance_level0_int;
    Real            loadbalance_max_fac;

    int             n_comp;
    IntVect         n_grow;

    //
    // Static class members.  Set defaults in Initialize()!!!
    //
    static bool initialized;
    //! Array of BoxArrays read in to initially define grid hierarchy
    static Vector<BoxArray> initial_ba;
    //! Array of BoxArrays read in to externally define grid hierarchy at each regrid
    static Vector<BoxArray> regrid_ba;

};

#endif /*_Amr_H_*/
