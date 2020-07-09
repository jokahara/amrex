
#include <string>
#include <limits>
#include <memory>

#ifdef _OPENMP
#include <omp.h>
#endif

#include <AMReX_AmrCore.H>
#include "CellFabArray.hpp"
//#include <AMReX_FluxRegister.H>
#include <AMReX_BCRec.H>

using namespace amrex;

class AmrGrid : public AmrCore
{
public:

    ////////////////
    // public member functions

    // constructor - reads in parameters from inputs file
    //             - sizes multilevel arrays and data structures
    AmrGrid ();
    virtual ~AmrGrid();

    // initializes multilevel data
    void InitData ();

    // advance a level by dt
    // includes a recursive call for finer levels
    void TimeStep (int lev, Real time, int iteration);

    // The required virtual function from AmrCore:

    // Make a new level using provided BoxArray and DistributionMapping and
    // fill with interpolated coarse level data.
    // overrides the pure virtual function in AmrCore
    virtual void MakeNewLevelFromCoarse (int lev, Real time, const BoxArray& ba,
					 const DistributionMapping& dm) override;

    // Remake an existing level using provided BoxArray and DistributionMapping and
    // fill with existing fine and coarse data.
    // overrides the pure virtual function in AmrCore
    virtual void RemakeLevel (int lev, Real time, const BoxArray& ba,
			      const DistributionMapping& dm) override;

    // Delete level data
    // overrides the pure virtual function in AmrCore
    virtual void ClearLevel (int lev) override;

    // Make a new level from scratch using provided BoxArray and DistributionMapping.
    // Only used during initialization.
    // overrides the pure virtual function in AmrCore
    virtual void MakeNewLevelFromScratch (int lev, Real time, const BoxArray& ba,
					  const DistributionMapping& dm) override;

    // tag all cells for refinement
    // overrides the pure virtual function in AmrCore
    virtual void ErrorEst (int lev, TagBoxArray& tags, Real time, int ngrow) override;

private:

    ////////////////
    // private member functions

    // read in some parameters from inputs file
    void ReadParameters();

    // compute a new multifab by coping in phi from valid region and filling ghost cells
    // works for single level and 2-level cases (fill fine grid ghost by interpolating from coarse)
    void FillPatch (int lev, Real time, CellFabArray& mf, int icomp, int ncomp);

    // fill an entire multifab by interpolating from the coarser level
    // this comes into play when a new level of refinement appears
    void FillCoarsePatch (int lev, Real time, CellFabArray& mf, int icomp, int ncomp);

    // utility to copy in data into another multifab
    void GetData (int lev, Real time, Vector<CellFabArray*>& data,
                  Vector<Real>& datatime);

    void ComputeDt ();

    // get plotfile name
    std::string PlotFileName (int lev) const;

    // put together an array of multifabs for writing
    Vector<const CellFabArray*> PlotFileMF () const;

    // set plotfile variables names
    Vector<std::string> PlotFileVarNames () const;

    // write plotfile to disk
    void WritePlotFile () const;

    ////////////////
    // private data members

    Vector<int> istep;      // which step?
    Vector<int> nsubsteps;  // how many substeps on each level?

    // keep track time, and time step at each level
    Vector<Real> t;
    Vector<Real> dt;

    // array of multifabs to store the solution at each level of refinement
    // after advancing a level we use "swap".
    Vector<CellFabArray> data;

    // BCRec is essentially a 2*DIM integer array storing the physical boundary
    // condition types at the lo/hi walls in each direction
    Vector<BCRec> bcs;  // 1-component

    ////////////////
    // runtime parameters

    // maximum number of steps and stop time
    int max_step = std::numeric_limits<int>::max();
    Real stop_time = std::numeric_limits<Real>::max();

    // if >= 0 we restart from a checkpoint
    std::string restart_chkfile = "";

    // how often each level regrids the higher levels of refinement
    // (after a level advances that many time steps)
    int regrid_int = 2;
    
    // plotfile prefix and frequency
    std::string plot_file {"plt"};
    int plot_int = -1;

    // checkpoint prefix and frequency
    std::string chk_file {"chk"};
    int chk_int = -1;
};