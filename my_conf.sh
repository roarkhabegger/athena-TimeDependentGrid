# Easy to edit configure script for athena++ 
python configure.py\
			 --prob blast\
			 --flux roe\
			 --coord cartesian \
			 --nghost=2 \
			 -debug \
			 -hdf5 \
			 --include /usr/include/hdf5/openmpi/ \
			 --fluxcl hlle\
			 -cl \
			 -mpi \
			 -omp 
			 #--include /usr/include/hdf5/openmpi/ \
#			 --hdf5_path=/usr/lib/x86_64-linux-gnu/
