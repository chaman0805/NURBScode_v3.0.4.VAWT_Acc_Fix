# NURBScode_v3.0.4.VAWT_Acc_Fix
Source Code Changes (.f90 files)
code/driver.f90
Three targeted fixes for the VAWT acceleration problem:

Enabled shell/blade solving — solshell = .false. and NBlade = 0 were replaced with solshell = myid.eq.0 and NBlade = 1, re-activating the structural blade solver on rank 0.

Fixed MPI broadcast scope — MPI_BARRIER + MPI_BCAST for xrotOld were moved outside the if(solshell) block (lines ~161–170), so all MPI ranks receive the broadcast regardless of whether they are the shell-solving rank. Previously, non-shell ranks would never receive the updated rotation, causing a hang or data divergence.

Fixed mismatched end if — an end if was added after call get_omega(...) and the if(solshell) guard was moved down to wrap only the rigid-body matrix section (Rfix, etc.), correcting a structural error in the time-stepping loop.

code/modules.f90
Integers NGAUSS, NNODE, NEL, maxNSHL, NNODE_LOC in the defs_shell module were given explicit default initializations to 0. This prevents undefined-variable errors when solshell is .false. on non-master ranks.
code/shell_fem_find_elm.f90
Array slice notation tightened: FEM%B_NET_D(...) accesses now use explicit 1:3 column slices instead of the full array, avoiding a potential size-mismatch when B_NET_D has more than 3 columns.
code/solveFlow.f90
Decoupled the structural solve from the solshell guard: solveKLShell is now called whenever shel is true (on all participating ranks), while get_omega is guarded with if (solshell) — so angular velocity is only computed on the rank that owns the shell geometry, then broadcast to others.
