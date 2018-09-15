#ifndef TASK_LIST_TASK_LIST_HPP_
#define TASK_LIST_TASK_LIST_HPP_
//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//!   \file task_list.hpp
//    \brief provides functionality to control dynamic execution using tasks

#include <stdint.h>
#include <string>

// Athena++ headers
#include "../athena.hpp"

// use 128-bit integers for this task list 
typedef unsigned __int128 uint128_t; 

// forward declarations
class Mesh;
class MeshBlock;
class TaskList;
class GravitySolverTaskList;

// return codes for functions working on individual Tasks and TaskList
enum TaskStatus {TASK_FAIL, TASK_SUCCESS, TASK_NEXT};
enum TaskListStatus {TL_RUNNING, TL_STUCK, TL_COMPLETE, TL_NOTHING_TO_DO};

//----------------------------------------------------------------------------------------
//! \struct IntegratorWeight
//  \brief weights used in time integrator tasks

struct IntegratorWeight {
  // 2S or 3S* low-storage RK coefficients, Ketchenson (2010)
  Real delta; // low-storage coefficients to avoid double F() evaluation per substage
  Real gamma_1, gamma_2, gamma_3; // low-storage coeff for weighted ave of registers
  Real beta; // Coefficients from bidiagonal Shu-Osher form Beta matrix, -1 diagonal terms
};

//----------------------------------------------------------------------------------------
//! \struct Task
//  \brief data and function pointer for an individual Task

struct Task {
  uint128_t task_id;      // encodes step & task using bit positions in HydroTasks
  uint128_t dependency;   // encodes dependencies to other tasks using " " " "
  enum TaskStatus (TaskList::*TaskFunc)(MeshBlock*, int);  // ptr to member function
};


//---------------------------------------------------------------------------------------
//! \class TaskState
//  \brief container for task states

class TaskState {
  public:
  uint128_t finished_tasks;
  int indx_first_task, num_tasks_left;
  void Reset(int ntasks) {
    indx_first_task = 0;
    num_tasks_left = ntasks;
    finished_tasks = (uint128_t)0LL;
  }
};


//----------------------------------------------------------------------------------------
//! \class TaskList
//  \brief data and function definitions for task list base class

class TaskList {
friend class TimeIntegratorTaskList;
friend class GravitySolverTaskList;
//friend class ClessTimeIntegratorTaskList;
public:
  explicit TaskList(Mesh *pm);
  virtual ~TaskList();

  // data
  int ntasks;     // number of tasks in this list
  int nsub_steps; // number of times task list should be repeated per full time step

  // functions
  enum TaskListStatus DoAllAvailableTasks(MeshBlock *pmb, int step, TaskState &ts);
  void DoTaskListOneSubstep(Mesh *pmesh, int step);
	// cless functions 
	//enum TaskListStatus DoAllAvailableTasksCL(MeshBlock *pmb, int step, TaskState &ts);
  //void DoTaskListOneSubstepCL(Mesh *pmesh, int step);

private:
  Mesh* pmy_mesh_;
  struct Task task_list_[128];
	//struct Task task_listCL_[64];
};

//----------------------------------------------------------------------------------------
//! \class TimeIntegratorTaskList
//  \brief data and function definitions for TimeIntegratorTaskList derived class

class TimeIntegratorTaskList : public TaskList {
public:
  TimeIntegratorTaskList(ParameterInput *pin, Mesh *pm);
  ~TimeIntegratorTaskList() {}

  // data
  std::string integrator;
  Real cfl_limit; // dt stability limit for the particular time integrator + spatial order
  struct IntegratorWeight step_wghts[MAX_NSTEP];

  void AddTimeIntegratorTask(uint128_t id, uint128_t dep);

  // functions
  enum TaskStatus StartAllReceive(MeshBlock *pmb, int step);
  enum TaskStatus ClearAllBoundary(MeshBlock *pmb, int step);

  enum TaskStatus CalculateFluxes(MeshBlock *pmb, int step);
  enum TaskStatus CalculateEMF(MeshBlock *pmb, int step);

  enum TaskStatus FluxCorrectSend(MeshBlock *pmb, int step);
  enum TaskStatus EMFCorrectSend(MeshBlock *pmb, int step);

  enum TaskStatus FluxCorrectReceive(MeshBlock *pmb, int step);
  enum TaskStatus EMFCorrectReceive(MeshBlock *pmb, int step);

  enum TaskStatus HydroIntegrate(MeshBlock *pmb, int step);
  enum TaskStatus FieldIntegrate(MeshBlock *pmb, int step);

  enum TaskStatus HydroSourceTerms(MeshBlock *pmb, int step);

  enum TaskStatus HydroDiffusion(MeshBlock *pmb, int step);
  enum TaskStatus FieldDiffusion(MeshBlock *pmb, int step);
  enum TaskStatus CalcDiffusivity(MeshBlock *pmb, int step);

	enum TaskStatus InternalEnergySync(MeshBlock *pmb, int step); 
	enum TaskStatus InternalEnergyCheck(MeshBlock *pmb, int step); 

  enum TaskStatus HydroSend(MeshBlock *pmb, int step);
  enum TaskStatus FieldSend(MeshBlock *pmb, int step);

  enum TaskStatus HydroReceive(MeshBlock *pmb, int step);
  enum TaskStatus FieldReceive(MeshBlock *pmb, int step);

  enum TaskStatus HydroShearSend(MeshBlock *pmb, int step);
  enum TaskStatus HydroShearReceive(MeshBlock *pmb, int step);
  enum TaskStatus FieldShearSend(MeshBlock *pmb, int step);
  enum TaskStatus FieldShearReceive(MeshBlock *pmb, int step);
  enum TaskStatus EMFShearSend(MeshBlock *pmb, int step);
  enum TaskStatus EMFShearReceive(MeshBlock *pmb, int step);
  enum TaskStatus EMFShearRemap(MeshBlock *pmb, int step);

  enum TaskStatus Prolongation(MeshBlock *pmb, int step);
  enum TaskStatus Primitives(MeshBlock *pmb, int step);
  enum TaskStatus PhysicalBoundary(MeshBlock *pmb, int step);
  enum TaskStatus UserWork(MeshBlock *pmb, int step);
  enum TaskStatus NewBlockTimeStep(MeshBlock *pmb, int step);
  enum TaskStatus CheckRefinement(MeshBlock *pmb, int step);

  enum TaskStatus GravSend(MeshBlock *pmb, int step);
  enum TaskStatus GravReceive(MeshBlock *pmb, int step);
  enum TaskStatus GravSolve(MeshBlock *pmb, int step);
  enum TaskStatus GravFluxCorrection(MeshBlock *pmb, int step);

  enum TaskStatus StartupIntegrator(MeshBlock *pmb, int step);
  enum TaskStatus UpdateTimeStep(MeshBlock *pmb, int step);

	// cless tasks  
	enum TaskStatus ClessCalculateFluxes(MeshBlock *pmb, int step); 
	enum TaskStatus ClessIntegrate(MeshBlock *pmb, int step); 
	enum TaskStatus ClessSourceTerms(MeshBlock *pmb, int step); 
	enum TaskStatus ClessSend(MeshBlock *pmb, int step);
	enum TaskStatus ClessReceive(MeshBlock *pmb, int step);
	enum TaskStatus ClessProlongation(MeshBlock *pmb, int step);
  enum TaskStatus ClessPrimitives(MeshBlock *pmb, int step);
  enum TaskStatus ClessPhysicalBoundary(MeshBlock *pmb, int step);
	enum TaskStatus ClessFluxCorrectSend(MeshBlock *pmb, int step);
	enum TaskStatus ClessFluxCorrectReceive(MeshBlock *pmb, int step);
	enum TaskStatus ClessStartupIntegrator(MeshBlock *pmb, int step); 
	enum TaskStatus ClessNewBlockTimeStep(MeshBlock *pmb, int step);
};


//----------------------------------------------------------------------------------------
// 128-bit integers with "1" in different bit positions used to ID each hydro task.

namespace HydroIntegratorTaskNames {
  const uint128_t NONE=0;
  const uint128_t START_ALLRECV=1LL<<0;
  const uint128_t CLEAR_ALLBND=1LL<<1;

  const uint128_t CALC_HYDFLX=1LL<<2;
  const uint128_t CALC_FLDFLX=1LL<<3;
  const uint128_t CALC_RADFLX=1LL<<4;
  const uint128_t CALC_CHMFLX=1LL<<5;

  const uint128_t ADD_VISCFLX=1LL<<6;
  const uint128_t ADD_HEATFLX=1LL<<7;
  const uint128_t ADD_OHMFLX=1LL<<8;
  const uint128_t ADD_ADFLX=1LL<<9;
  const uint128_t ADD_HALLFLX=1LL<<10;

  const uint128_t SEND_HYDFLX=1LL<<11;
  const uint128_t SEND_FLDFLX=1LL<<12;
  const uint128_t SEND_RADFLX=1LL<<13;
  const uint128_t SEND_CHMFLX=1LL<<14;

  const uint128_t RECV_HYDFLX=1LL<<15;
  const uint128_t RECV_FLDFLX=1LL<<16;
  const uint128_t RECV_RADFLX=1LL<<17;
  const uint128_t RECV_CHMFLX=1LL<<18;

  const uint128_t SRCTERM_HYD=1LL<<19;
  const uint128_t SRCTERM_FLD=1LL<<20;
  const uint128_t SRCTERM_RAD=1LL<<21;
  const uint128_t SRCTERM_CHM=1LL<<22;

  const uint128_t INT_HYD=1LL<<23;
  const uint128_t INT_FLD=1LL<<24;
  const uint128_t INT_RAD=1LL<<25;
  const uint128_t INT_CHM=1LL<<26;

  const uint128_t SEND_HYD=1LL<<27;
  const uint128_t SEND_FLD=1LL<<28;
  const uint128_t SEND_RAD=1LL<<29;
  const uint128_t SEND_CHM=1LL<<30;

  const uint128_t RECV_HYD=1LL<<31;
  const uint128_t RECV_FLD=1LL<<32;
  const uint128_t RECV_RAD=1LL<<33;
  const uint128_t RECV_CHM=1LL<<34;

  const uint128_t PROLONG =1LL<<35;
  const uint128_t CON2PRIM=1LL<<36;
  const uint128_t PHY_BVAL=1LL<<37;
  const uint128_t USERWORK=1LL<<38;
  const uint128_t NEW_DT  =1LL<<39;
  const uint128_t AMR_FLAG=1LL<<40;

  const uint128_t SOLV_GRAV=1LL<<41;
  const uint128_t SEND_GRAV=1LL<<42;
  const uint128_t RECV_GRAV=1LL<<43;
  const uint128_t CORR_GFLX=1LL<<44;

  const uint128_t STARTUP_INT=1LL<<45;
  const uint128_t UPDATE_DT  =1LL<<46;

  const uint128_t SEND_HYDSH=1LL<<47;
  const uint128_t SEND_EMFSH=1LL<<48;
  const uint128_t SEND_FLDSH=1LL<<49;
  const uint128_t RECV_HYDSH=1LL<<50;
  const uint128_t RECV_EMFSH=1LL<<51;
  const uint128_t RECV_FLDSH=1LL<<52;
  const uint128_t RMAP_EMFSH=1LL<<53;

  const uint128_t DIFFUSE_HYD=1LL<<54;
  const uint128_t DIFFUSE_FLD=1LL<<55;
  const uint128_t CALC_DIFFUSIVITY=1LL<<56;

	const uint128_t SYNC_IE=1LL<<57;
	const uint128_t CHECK_IE=1LL<<58;

	// cless tasks 
	const uint128_t CL_CALC_FLX =(uint128_t)1LL<<59;
	const uint128_t CL_INT      =(uint128_t)1LL<<60;
	const uint128_t CL_SRCTERM  =(uint128_t)1LL<<61;
	const uint128_t CL_SEND     =(uint128_t)1LL<<62;
	const uint128_t CL_RECV     =(uint128_t)1LL<<63;
	const uint128_t CL_PROLONG  =(uint128_t)1LL<<64; 
	const uint128_t CL_CON2PRIM =(uint128_t)1LL<<65;
	const uint128_t CL_PHY_BVAL =(uint128_t)1LL<<66;
	const uint128_t CL_SEND_FLX =(uint128_t)1LL<<67;
	const uint128_t CL_RECV_FLX =(uint128_t)1LL<<68; 
	const uint128_t CL_START_INT=(uint128_t)1LL<<69;
	const uint128_t CL_NEW_DT   =(uint128_t)1LL<<70; 

}; // namespace HydroIntegratorTaskNames


#endif // TASK_LIST_TASK_LIST_HPP_
