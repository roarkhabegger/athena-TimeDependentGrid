# Easy to edit configure script for athena++ 
python configure.py\
			 --prob blast\
			 --flux hlle\
			 --nghost=3 \
			 -debug \
			 -hdf5 \
			 -cl \
			 -mpi \
			 --include /usr/include/hdf5/openmpi/ \
			# --include /usr/include/hdf5/serial/ \
#			 --hdf5_path=/usr/lib/x86_64-linux-gnu/
