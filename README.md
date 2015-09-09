# runoff_tools
Some tools to make the generation of runoff files fo MOM[45] a breeze...

create_model_coast.f90   ! Produces a file which contains the points adjacent to land.

create_model_wet.f90    ! Produce a file containing all the wet points of a from a mosaic description of of a grid.

create_runoff_nn.f90     ! Create files with nearest neighbour of the source runoff which lie on the coast.

create_runoff_weights.f90 ! Spread the runoff to nearest coastal neighbors of the points found in create_nn.f90
                           Either get N nearest neighbours or all the neighbours within a radius
                           Weights can be evenly distributed or based on a gaussian weighting based on the distance
                           from either the source or its nearest neighbour. Requires the results from running
                           create_model_coast.
                           
create_runoff_weights_spread.f90 ! Spread the runoff to any wet neighbour of the target point. Requires output from
                           create_model_wet. This allows spread out to open ocean which may be necessary for high 
                           volume discharge like the Amazon River. The amount of spread can be modified with weighting.
                           
process_runoff.f90        ! Create the runoff files using masks and weights from above.

kdtree2.f90               ! module from the wonderful KDTREE 2 package written by Matthew Kennel.
                           http://arxiv.org/abs/physics/0408067
                           https://github.com/jmhodges/kdtree2
                           
                           
General notes.

Creating the connection file (runouff source to model coast).

A kdtree is used to identify the nearest neighbours to the source. This requires a mask of the source grid based
on the nonzero runoff. I calculate a time invariant mask (runoff_mask.nc) based on the whole time series. You could
recalculate the nearest neighbours for every time step as it's so quick but this allows for more flexibility.
There's also the problem of making sure the correct nearest neighbours are picked (e.g. wrong side of a peninsular,
think Mekong River).  We also want all the coastal points of the model.

I save the results (connected indices, positions and the model coast) in an intermediate file. This allows the user
to check the connections.

To look at the connections between points here's some example Ferret code.

=========================================

use runoff_connection_nn.nc

def axis/y/from_data/name=ygr j[j=1:3]

let MASK = IF J LE 2 THEN 0*SOURCE_I + 0*J[GY=YGR]+1

let/title="xpoints"  PNT_X = IF J EQ 1 THEN MOD(SOURCE_X+360,360)*MASK ELSE MOD(TARGET_X+360,360)*MASK

let title="ypoints" PNT_Y = IF J EQ 1 THEN SOURCE_Y*MASK ELSE TARGET_Y*MASK

! Plot connections around Indochina

polygon/j=1:3/i=1:3987/line/key/hlim=80:120/vlim=10:30/col=red pnt_x,pnt_y

go land_detail

! Overlay with source points

plot/ov/sym=1/vs/col=blue source_x,source_y

! Overlay with ALL coastal points

plot/ov/sym=11/vs/col=green coast_x,coast_y

! Plot circles for the actual target NNs.

 plot/sym=27/vs/col=purple/ov mod(target_x+360,360),target_y
 
 ======================================

Creating the weights.

This takes the file  runoff_connection_nn.nc  and again builds a kdtree for the model coastal/ocean points. This time
we find either the nearest N neighbours or the nearest neighbours within a certain radius of the nearest neighbours
identified in step 1. Weights are either even or can be distance weighted.

Transferring the runoff to the model grid.

Each timestep is read in and the total runoff for each point is calculated according to our weights. I write out using
Netcdf4 compression, deflation level 1 with chunking of 200 point tiles.

The size of the new files. 60 years of the Dai and Trenberth runoff output for the GFDL 0.1 degree tripolar grid.
This took about 2 minutes on a single processor.

[raf599@raijin6 gfdl_grid]$ ls -lh new_run.nc runoff_connection_nn.nc runoff_weights.nc

-rw-r----- 1 raf599 p93 424M Aug 26 17:15 new_run.nc

-rw-r----- 1 raf599 p93 3.7M Aug 26 15:57 runoff_connection_nn.nc

-rw-r----- 1 raf599 p93 2.5M Aug 26 15:57 runoff_weights.nc
