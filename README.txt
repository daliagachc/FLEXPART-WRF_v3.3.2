J. Brioude, Nov  13  2014
**************************************************************
If FLEXPART-WRF is used in your publication, please cite:
Brioude, J., D. Arnold, A. Stohl et al. (2013): The Lagrangian particle dispersion model FLEXPART-WRF version 3.1, Geosci. Model Dev., 6, 1889-1904
**************************************************************

To compile flexwrf, choose your compiler in makefile.mom (line 23), the path to the NetCDF library and then type
make -f makefile.mom mpi  for MPI+OPENMP hybrid run
make -f makefile.mom omp  for OPENMP parallel run
make -f makefile.mom serial for a serial run
********************************************************************
To run flexwrf, you can pass an argument to the executable that gives the name of the input file.
for instance
./flexwrf31_mpi /home/jbrioude/inputfile.txt
Otherwise, the file flexwrf.input in the current directory is read by default.

Examples of forward and backward runs are available in the examples directory.


*****************************************************************
Versions timeline
version 3.2: bug fix: particles are now dumped in binary format (the partposit* files) if netcdf output format is asked.
             bug fix on  memory leaks when netcdf output are used
             bug fix on  netcdf output when ioutforeachrelease is used in forward runs
             bug fix in readwind*.f90, read_ncwrfout.f90 gridcheck*.f90 that affect netcdf reading (Thanks to Alan Griffiths) 
version 3.1: bug fix on the sign of sshf in readwind.f90
	     modifications of advance.f90 to limit the vertical velocity from cbl scheme
	     bug fix in write_ncconc.f90
	     modifications of interpol*f90 routines to avoid crashes using tke_partition_hanna.f90 and tke_partition_my.f90
version 3.0    First public version 

version 2.4.1: New modifications on the wet deposition scheme from Petra Seibert

version 2.3.1: a NetCDF format output is implemented. 

version 2.2.7: CBL scheme is implemented. a new random generator is implemented.

version 2.0.6:
-map factors are used in advance.f90 when converting the calculated distance
into a WRF grid distance. 
-fix on the divergence based vertical wind

version 2.0.5:
the time over which the kernel is not used has been reduced from 10800 seconds
to 7200 seconds. Those numbers depend on the horizontal resolution, and a more
flexible solution might come up in a future version
version 2.0.4:
- bug fix for regular output grid
- IO problems in ASCII have been fixed
- add the option of running flexpart with an argument that gives the name of
  the inputfile instead of flexwrf.input
version 2.0.3:
- bug fix when flexpart is restarted.
-bug fix in coordtrafo.f90
- a new option that let the user decide if the time for the the time average
  fields from WRF has to be corrected or not.

version 2.0.2:
- bug fix in sendint2_mpi_old.f90
- all the *mpi*.f90 have been changed to handle more properly the memory.
- timemanager_mpi has changed accordingly. Some bug fix too
- bug fix in writeheader
- parallelization of calcpar and verttransform.f90, same for the nests.

version 2.0.1:
-1 option added in flexwrf.input to define the output grid with dxout and dyout
-fix in readinput.f90 to calculate maxpart more accurately

version 2.0: first OPENMP/MPI version

version 1.0: 
This is a fortran 90 version of FLEXPART.
Compared to PILT, the version from Jerome Fast available on the NILU flexpart website, several bugs and improvements have been made (not
necessarily commented) in the subroutines. 
non exhaustive list:
1) optimization of the kein-fritch convective scheme (expensive)
2) possibility to output the flexpart run in a regular lat/lon output grid.
flexwrf.input has 2 options to let the model know which coordinates are used
for the output domaine and the release boxes.
3) Differences in earth radius between WRF and WRF-chem is handled.
4) time averaged wind, instantaneous omega or a vertical velocity internally calculated in FLEXPART can be used now.
5) a bug fix in pbl_profile.f due to the variable kappa.

Turb option 2 and 3 from Jerome Fast's version lose mass in the model. Those
options are not recommended.

***********************************************************************
General comments on The hybrid version of flexpart wrf:
This version includes a parallelized hybrid version of FLEXPART that can be
used with:
- 1 node (1 computer) with multi threads using openmp in shared memory, 
- or several nodes (computers) in distributed memory (using mpi) and several threads in shared memory (using openmp).
if a mpi library is not available with your compiler, use makefile.nompi to compile flexwrf

The system variable OMP_NUM_THREADS has to be set before running the model to define the number of thread used. 
it can also be fixed in timemanager*f90. 
If not, flexwrf20_mpi will use 1 thread.

When submitting a job to several nodes, mpiexec or mpirun needs to know that 1 task has to be allocated per node to let openmp doing the work within each node in shared memory.
See submit.sh as an example. 

Compared to the single node version, this version includes modifications of:

- flexwrf.f90 that is renamed into flexwrf_mpi.f90
- timemanager.f90 that is renamed into timemanager_mpi.f90
- the interpol*f90 and hanna* has been modified.
- the routines *mpi*.f90 are used to send or receive data between nodes.

The most important modifications are in timemanager_mpi.f90, initialize.f90 and advance.f90.
search for JB in timemanager_mpi.f90 to have additional comments.
in advance.f90, I modified the way the random number is picked up (line 187). I use a simple count and the id of the thread instead of the random pick up that uses ran3.
If the series of random number is output for a give release box (uncomment lines 195 to 198), the distribution is quite good, and I don't see any bigger bias that the one in the single thread version.
of course, the distribution is less and less random when you increase the number of nodes or threads.


*********************************************************
performance:
this is the performance of the loop line 581 in timemanager_mpi.f90 that calculates the trajectories.
I use the version v74 as the reference (single thread, fortran 77).
There is a loss in performance between v74 and v90 because of the temporary variables th_* that has to be used as private variables in timemanager_mpi.f90
		v74
v90 1thread	0.96
v90 2threads	1.86
v90 4threads	3.57
v90 8threads	6.22

performance of the communication between nodes:
depends on the system. The super computer that I use can transfer about 1Gb in 1 second.
in timemanager_mpi.f90, the output lines 540 and 885 give the time needed by the system to communicate between nodes. using 100 millions particles and say 4 nodes, it takes about 1 second.

