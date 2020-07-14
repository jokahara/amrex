
#ifndef AmrLevel_H_
#define AmrLevel_H_

#include <AMReX_REAL.H>
#include <AMReX_Geometry.H>
#include <AMReX_LayoutData.H>
#include <AMReX_DistributionMapping.H>

#include "Amr.h"
#include "CellFabArray.h"
#include "Derive.h"
#include "StateDescriptor.h"
#include "StateData.h"

#include <memory>
#include <map>

using namespace amrex;

/**
* \brief Virtual base class for managing individual levels.
* AmrLevel functions both as a container for state data on a level
* and also manages the advancement of data in time.
*/
class AmrLevel
{
    friend class Amr;
    friend class FillPatchIterator;
    friend class FillPatchIteratorHelper;

public:
    //! What time are we at?
    enum TimeLevel { AmrOldTime,
                     AmrHalfTime,
                     AmrNewTime,
                     Amr1QtrTime,
                     Amr3QtrTime,
                     AmrOtherTime };
    //! The destructor.
    virtual ~AmrLevel ();
    //! Get the level directory names
    void LevelDirectoryNames (const std::string &dir,
                              std::string &LevelDir,
                              std::string &FullPath);
    //! Create the Level_ directory for checkpoint and plot files
    virtual void CreateLevelDirectory (const std::string &dir);
    /**
    * \brief Set if the Level_ directory was created or to clear the value.
    * CreateLevelDirectory sets levelDirectoryCreated = true
    */
    void SetLevelDirectoryCreated(bool ldc) noexcept { levelDirectoryCreated = ldc; }
    /**
    * \brief A string written as the first item in writePlotFile() at
    * level zero. 
    * It is so we can distinguish between different types of
    * plot files.
    * This default "HyperCLaw-V1.1" is for VisIt software and some of our
    * internal postprocessing routines
    */
    virtual std::string thePlotFileType () const
    {
        static const std::string the_plot_file_type("HyperCLaw-V1.1");
        return the_plot_file_type;
    }
    /**
    * \brief Write plot file stuff to specified directory.  This is a
    * pure virtual function and hence MUST be implemented by
    * derived classes.
    */
    virtual void writePlotFile (const std::string& dir,
                                std::ostream&      os);

    //! Do pre-plotfile work to avoid synchronizations while writing the amr hierarchy
    virtual void writePlotFilePre (const std::string& dir,
                                   std::ostream&      os);

    //! Do post-plotfile work to avoid synchronizations while writing the amr hierarchy
    virtual void writePlotFilePost (const std::string& dir,
                                    std::ostream&      os);

    /**
    * \brief Write small plot file stuff to specified directory.  
    * Unlike writePlotFile, this is NOT a pure virtual function 
    * so implementation by derived classes is optional.
    */
    virtual void writeSmallPlotFile (const std::string& dir,
                                     std::ostream&      os) {};
    //! Write current state to checkpoint file.
    virtual void checkPoint (const std::string& dir,
                             std::ostream&      os,
                             bool               dump_old = true);
    //! Do pre-checkpoint work to avoid synchronizations while writing the amr hierarchy
    virtual void checkPointPre (const std::string& dir,
                                std::ostream&      os);
    //! Do post-checkpoint work to avoid synchronizations while writing the amr hierarchy
    virtual void checkPointPost (const std::string& dir,
                                 std::ostream&      os);
    //! Restart from a checkpoint file.
    virtual void restart (Amr&          papa,
                          std::istream& is,
			  bool          bReadSpecial = false);

    //! Old checkpoint may have different number of states than the new source code.
    virtual void set_state_in_checkpoint (Vector<int>& state_in_checkpoint);

    //! Is name a state variable?
    static bool isStateVariable (const std::string& name,
                                 int&               state_indx,
                                 int&               ncomp);

    static void FlushFPIcache ();
    /**
    * \brief Compute the initial time step.  This is a pure virtual function
    * and hence MUST be implemented by derived classes.
    */
    virtual void computeInitialDt (int                    finest_level,
                                   int                    sub_cycle,
                                   Vector<int>&           n_cycle,
                                   const Vector<IntVect>& ref_ratio,
                                   Vector<Real>&          dt_level,
                                   Real                   stop_time) = 0;
    /**
    * \brief Compute the next time step.  This is a pure virtual function
    * and hence MUST be implemented by derived classes.
    */
    virtual void computeNewDt (int                    finest_level,
                               int                    sub_cycle,
                               Vector<int>&           n_cycle,
                               const Vector<IntVect>& ref_ratio,
                               Vector<Real>&          dt_min,
                               Vector<Real>&          dt_level,
                               Real                   stop_time,
                               int                    post_regrid_flag) = 0;
    /**
    * \brief Do an integration step on this level.  Returns maximum safe
    * time step.  This is a pure virtual function and hence MUST
    * be implemented by derived classes.
    */
    virtual Real advance (Real time,
                          Real dt,
                          int  iteration,
                          int  ncycle) = 0;

    /**
    * \brief Contains operations to be done after a timestep.  This is a
    * pure virtual function and hence MUST be implemented by derived
    * classes.
    */
    virtual  void post_timestep (int iteration) = 0;
    /**
    * \brief Contains operations to be done only after a full coarse
    * timestep.  The default implementation does nothing.
    */
    virtual void postCoarseTimeStep (Real time);
    /**
    * \brief Operations to be done after restart.
    */
    virtual void post_restart () {};
    /**
    * \brief Operations to be done after regridding
    * This is a pure virtual function and hence MUST be
    * implemented by derived classes.
    */
    virtual  void post_regrid (int lbase,
                               int new_finest) = 0;
    /**
    * \brief Operations to be done after initialization.
    * This is a pure virtual function and hence MUST be
    * implemented by derived classes.
    */
    virtual  void post_init (Real stop_time) = 0;
    /**
    * \brief Is it ok to continue the calculation?
    */
    virtual  int okToContinue () { return 1; }
    /**
    * \brief Should I regrid with this level as base level?
    * This test is only evaluated if regrid_int > 0 and 
    * level_count >= regrid_int as well. Defaults to true.
    */
    virtual  int okToRegrid ();
    /**
    * \brief Init grid data at problem start-up.
    * This is a pure virtual function and hence MUST be
    * implemented by derived classes.
    */
    virtual void initData () = 0;
    //! Set the time levels of state data.
    virtual void setTimeLevel (Real time,
                               Real dt_old,
                               Real dt_new);
    //! Alloc space for old time data.
    virtual void allocOldData ();
    //! Delete old-time data.
    virtual void removeOldData ();
    /**
    * \brief Init data on this level from another AmrLevel (during regrid).
    * This is a pure virtual function and hence MUST be
    * implemented by derived classes.
    */
    virtual void init (AmrLevel &old) = 0;
    /**
    * Init data on this level after regridding if old AmrLevel
    * did not previously exist. This is a pure virtual function
    * and hence MUST be implemented by derived classes.
    */
    virtual void init () = 0;
    //! Reset data to initial time by swapping new and old time data.
    void reset ();
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
    //! Number of states at this level.
    int numStates () const noexcept { return state.size(); }
    //! Returns the indices defining physical domain.
    const Box& Domain () const noexcept { return geom.Domain(); }
    //! Timestep n at this level.
    int nStep () const noexcept { return parent->levelSteps(level); }
    //! Returns the geometry object.
    const Geometry& Geom () const noexcept { return geom; }
    //
    const IntVect& fineRatio () const noexcept { return fine_ratio; }
    //! Returns number of cells on level.
    Long countCells () const noexcept;

    //! Get the area not to tag.
    const BoxArray& getAreaNotToTag() noexcept;
    const Box& getAreaToTag() noexcept;
    //! Constuct the area not to tag.
    void constructAreaNotToTag();
    //! Set the area not to tag.
    void setAreaNotToTag(BoxArray& ba) noexcept;

    /**
    * \brief Error estimation for regridding. This is a pure virtual
    * function and hence MUST be implemented by derived classes.
    */
    virtual void errorEst (TagBoxArray& tb,
                           int          clearval,
                           int          tagval,
                           Real         time,
			               int          n_error_buf = 0,
                           int          ngrow = 0) = 0;
    //! Interpolate from coarse level to the valid area in dest.
    void FillCoarsePatch (CellFabArray& dest,
                          int       dcomp,
                          Real      time,
                          int       state_idx,
                          int       scomp,
                          int       ncomp,
			              int       nghost = 0);
    //! Function to set physical boundary conditions.
    virtual void setPhysBoundaryValues (CellFab& dest,
                                        int      state_indx,
                                        Real     time,
                                        int      dest_comp,
                                        int      src_comp,
                                        int      num_comp);
    /**
    * \brief Returns a CellFabArray containing the derived data for this level.
    * The user is responsible for deleting this pointer when done
    * with it.  If ngrow>0 the CellFabArray is built on the appropriately
    * grown BoxArray.
    */
    virtual std::unique_ptr<CellFabArray> derive (const std::string& name,
					                           Real               time,
					                           int                ngrow);
    /**
    * \brief This version of derive() fills the dcomp'th component of mf
    * with the derived quantity.
    */
    virtual void derive (const std::string& name,
                         Real               time,
                         CellFabArray&         mf,
                         int                dcomp);
    //! State data object.
    StateData& get_state_data (int state_indx) noexcept { return state[state_indx]; }
    //! State data at old time.
    CellFabArray& get_old_data (int state_indx) noexcept { return state[state_indx].oldData(); }
    //! State data at old time.
    const CellFabArray& get_old_data (int state_indx) const noexcept { return state[state_indx].oldData(); }
    //! State data at new time.
    CellFabArray& get_new_data (int state_indx) noexcept { return state[state_indx].newData(); }
    //! State data at new time.
    const CellFabArray& get_new_data (int state_indx) const noexcept { return state[state_indx].newData(); }
    //! Returns list of Descriptors.
    static const DescriptorList& get_desc_lst () noexcept { return desc_lst; }
    //! Returns list of derived variables.
    static DeriveList& get_derive_lst () noexcept;
    //! Returns whether or not we want a post-timestep regrid.
    int postStepRegrid () noexcept { return post_step_regrid; }
    //! Sets a new value for the post-timestep regrid trigger.
    void setPostStepRegrid (int new_val) noexcept { post_step_regrid = new_val; }

    //! Update the distribution maps in StateData based on the size of the map
    void UpdateDistributionMaps ( DistributionMapping& dmap );

    //! Boundary condition access function.
    Vector<int> getBCArray (int State_Type,
                            int gridno,
                            int scomp,
                            int ncomp);
    //! Get state data at specified index and time.
    CellFabArray& get_data (int  state_indx, Real time) noexcept;
    //! Hack to allow override of (non-fine-fine) fillpatched boundary data
    virtual void set_preferred_boundary_values (CellFabArray& S,
                                                int       state_index,
                                                int       scomp,
                                                int       dcomp,
                                                int       ncomp,
                                                Real      time) const;
    /** 
    * \brief called in grid_places after other tagging routines to modify
    * the list of tagged points.  Default implementation does nothing.
    */
    virtual void manual_tags_placement (TagBoxArray&           tags,
                                        const Vector<IntVect>& bf_lev);
    //! Modify list of variables to be plotted
    virtual void setPlotVariables ();
    //! Modify list of variables to be plotted
    virtual void setSmallPlotVariables ();
    /**
    * \brief Estimate the amount of work required to advance Just this level
    * based on the number of cells.
    * This estimate can be overwritten with different methods
    */
    virtual Real estimateWork();

    //! Which state data type is for work estimates? -1 means none
    virtual int WorkEstType () { return -1; }

    /**
    * \brief Returns one the TimeLevel enums.
    * Asserts that time is between AmrOldTime and AmrNewTime.
    */ 
    TimeLevel which_time (int  state_indx, Real time) const noexcept;

    //! Does the AmrLevel want Amr to write a plotfile now?
    virtual bool writePlotNow ();

    //! Does the AmrLevel want Amr to write a small plotfile now?
    virtual bool writeSmallPlotNow ();

#ifdef AMREX_PARTICLES
    //! This function can be called from the parent 
    virtual void particle_redistribute (int lbase = 0, bool a_init = false) {;}
#endif

    static void FillPatch (AmrLevel& amrlevel,
                           CellFabArray& leveldata,
                           int       boxGrow,
                           Real      time,
                           int       index,
                           int       scomp,
                           int       ncomp,
                           int       dcomp=0);

    static void FillPatchAdd (AmrLevel& amrlevel,
                              CellFabArray& leveldata,
                              int       boxGrow,
                              Real      time,
                              int       index,
                              int       scomp,
                              int       ncomp,
                              int       dcomp=0);

protected:
    //! The constructors -- for derived classes.
    AmrLevel () noexcept;

    AmrLevel (Amr&            papa,
              int             lev,
              const Geometry& level_geom,
              const BoxArray& bl,
	          const DistributionMapping& dm,
              Real            time);

    AmrLevel (const AmrLevel&) = delete;
    AmrLevel& operator = (const AmrLevel&) = delete;

    //! Common code used by all constructors.
    void finishConstructor (); 

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
    static DeriveList     derive_lst;   // List of derived quantities.
    static DescriptorList desc_lst;     // List of state variables.
    Vector<StateData>     state;        // Array of state data.

    BoxArray              m_AreaNotToTag; //Area which shouldn't be tagged on this level.
    Box                   m_AreaToTag;    //Area which is allowed to be tagged on this level.

    int                   post_step_regrid; // Whether or not to do a regrid after the timestep.

    bool                  levelDirectoryCreated;    // for checkpoints and plotfiles

    std::unique_ptr<FabFactory<CellFab> > m_factory;

private:

    mutable BoxArray      edge_grids[AMREX_SPACEDIM];  // face-centered grids
    mutable BoxArray      nodal_grids;              // all nodal grids
};

//
// Forward declaration.
//
class FillPatchIteratorHelper;

class FillPatchIterator
    :
    public MFIter
{
  public:

    friend class AmrLevel;

    FillPatchIterator (AmrLevel&  amrlevel,
                       CellFabArray& leveldata);

    FillPatchIterator (AmrLevel& amrlevel,
                       CellFabArray& leveldata,
                       int       boxGrow,
                       Real      time,
                       int       state_indx,
                       int       scomp,
                       int       ncomp);

    void Initialize (int  boxGrow,
                     Real time,
                     int  state_indx,
                     int  scomp,
                     int  ncomp);

    ~FillPatchIterator ();

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

    void FillFromLevel0 (Real time, int index, int scomp, int dcomp, int ncomp);
    void FillFromTwoLevels (Real time, int index, int scomp, int dcomp, int ncomp);

    //
    // The data.
    //
    AmrLevel&                         m_amrlevel;
    CellFabArray&                        m_leveldata;
    std::vector< std::pair<int,int> > m_range;
    CellFabArray                         m_fabs;
    int                               m_ncomp;
};

class FillPatchIteratorHelper
{
public:

    friend class FillPatchIterator;

    FillPatchIteratorHelper (AmrLevel&  amrlevel,
                             CellFabArray& leveldata);

    FillPatchIteratorHelper (AmrLevel&     amrlevel,
                             CellFabArray&    leveldata,
                             int           boxGrow,
                             Real          time,
                             int           state_indx,
                             int           scomp,
                             int           ncomp,
                             Interpolater* mapper);

    void Initialize (int           boxGrow,
                     Real          time,
                     int           state_indx,
                     int           scomp,
                     int           ncomp,
                     Interpolater* mapper);

    ~FillPatchIteratorHelper ();

    void fill (CellFab& fab, int dcomp, int idx);

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
    AmrLevel&                   m_amrlevel;
    CellFabArray&                  m_leveldata;
    CellFabArrayCopyDescriptor     m_mfcd;
    Vector< Vector<CellFabArrayId> > m_mfid;     // [level][oldnew]
    Interpolater*               m_map;
    std::map<int,Box>           m_ba;
    Real                        m_time;
    int                         m_growsize;
    int                         m_index;
    int                         m_scomp;
    int                         m_ncomp;
    bool                        m_FixUpCorners;

    std::map< int,Vector< Vector<Box> > >                 m_fbox; // [grid][level][validregion]
    std::map< int,Vector< Vector<Box> > >                 m_cbox; // [grid][level][fillablesubbox]
    std::map< int,Vector< Vector< Vector<FillBoxId> > > > m_fbid; // [grid][level][fillablesubbox][oldnew]
};

#endif /*_AmrLevel_H_*/
