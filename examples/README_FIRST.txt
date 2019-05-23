Those test cases can be run by copying a compiled binary and the flexwrf.input file of a test case in a directory. In addition, data from the directory src_flexwrf_v3.0/data are required for calculating dry deposition for the Cesium species.
The three test cases cover most of the format the input file should have if 2 output grid nests, ageclasses, several species with hourly emission variations are used.

an example of AVAILABLE file with the proper format is given for 2 different domains, assuming that each WRF output has 6 time frames with a time interval output of 30 minutes.
#################################
file flexwrf.input.forward1: 
- This case is a forward run (using LDIRECT=1) 
- The output gives the concentration and deposition from the total release of the 2 release location by using IOUTPUTFOREACHREL=0
- a passive airtracer and Cesium particles are released (NSPEC=2), without any hourly emission variation (EMITVAR=0)
- 2 output grid are defined. They are defined using latitude/longitude coordinates, the number of grid cells  and the grid spacing
- 2 ageclasses are asked by using LAGESPECTRA=1
- The NetCDF output format is selected, with 3 time frames per file.

Because the species #1 is a passive tracer, the dry and wet deposition fields
are empty. The wet and dry deposition fields exist for species #2

This test case can be run by typing
./flexwrf30_pgi_omp flexwrf.input.forward1
if, for instance, an OpenMP is compiled.
#################################
file flexwrf.input.forward2: 
- This case is a forward run (using LDIRECT=1)
- The output gives the concentration and deposition from each release location by using IOUTPUTFOREACHREL=1
- a passive airtracer and Cesium particles are released (NSPEC=2), with hourly and daily emission variation (EMITVAR=1)
- 1 output grid is defined. It is defined base on metric coordinates of the WRF output
- No ageclass is used because  LAGESPECTRA=0
- The NetCDF output format is selected, with 1 time frame per file. because IOUTPUTFOREACHREL=1, an additional dimension is used in the netcdf file compared to the case flexwrf.input.forward1

This test case can be run by typing
./flexwrf30_pgi_omp flexwrf.input.forward2
if, for instance, an OpenMP is compiled.
#################################
file flexwrf.input.backward1:
- This case is a backward run (using LDIRECT=-1)
- Only one species is emitted, the passive airtracer.
- The output gives the source-receptor relationship of the 3 release locations over Los Angeles Basin. Their vertical location is defined based on pressure
- 1 output grid is defined. 
- No ageclass is used because LAGESPECTRA=0
- The binary format is selected. The header file (binary format), latlon*.txt files (ascii format) and grid_time* (binary format) files are generated.

This test case can be run by typing
./flexwrf30_pgi_omp flexwrf.input.backward1
if, for instance, an OpenMP is compiled.

#################################
file flexwrf.input.backward2:
Identical to flexwrf.input.backward1 except that 2 WRF nested output are read. (The nested WRF output files is not available on http://www.flexpart.eu)
In such case, the OUTPUT grid can be defined in WRF coordinates but relative to the outer WRF domain.

