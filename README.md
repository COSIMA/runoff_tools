# runoff_tools
Some tools to make the generation of runoff files fo MOM[45] a breeze...

create_runoff_nn.f90     ! Creatle files with nearest neighbour of the source runoff which lie on the coast

create_runoff_weights.f90 ! Spread the runoff to nearest coastal neighbors of the points found in create_nn.f90
                           Either get N nearest neighbours or all the neighbours within a radius
                           Weights can be evenly distributed or based on a gaussian weighting based on the distance
                           from either the source or its nearest neighbour.

process_runoff.f90        ! Create the runoff files

kdtree2.f90               ! module from the wonderful KDTREE 2 package written by Matthew Kennel.
                           http://arxiv.org/abs/physics/0408067
                           https://github.com/jmhodges/kdtree2
                          
