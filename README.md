# NURBScode_v3.0.4.VAWT_Acc_Fix
code/driver.f90
Turned the blade solver back on. It was disabled (solshell = .false., NBlade = 0). Now it runs on the main processor (NBlade = 1).
Fixed a communication bug. The code tells all processors what the current rotation angle is. Before, only the blade processor got that message — the others were left waiting forever. Now everyone gets it.
Fixed a missing end if. A code block wasn't closed properly, so the program logic was wrong inside the main time loop.
code/modules.f90
Set starting values to zero. Some integer variables had no starting value, which could cause random garbage behavior on processors not running the blade solver. Now they start at 0.
code/shell_fem_find_elm.f90
Fixed an array size mismatch. When grabbing coordinates from an array, the code now explicitly asks for columns 1–3 instead of the whole row, preventing a crash if the array has extra columns.
code/solveFlow.f90
Split one task into two guarded tasks. Solving the blade structure and computing the rotation speed were bundled together. Now the blade solve runs on all relevant processors, but only the blade-owning processor computes the rotation speed — then shares it with everyone else.
