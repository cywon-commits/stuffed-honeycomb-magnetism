## DO NOT treat the 36-spin high-Sz pilot as a finished low-Sz eigensolver

The current archive includes `kernel36pilot.cpp` and `pilot36.py` only as verified high-magnetization tests and realistic geometry. A prospective HPC programmer should use the detailed execution spec to implement a full reduced-basis solver, run memory and Hermiticity tests, then all-sector GS search. There is **no honest one-command SLURM script to launch an existing 36-spin GS calculation** in the provided code. A scheduler template is supplied for the 18-spin reproduction only.
