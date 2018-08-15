# Easy to edit configure script for athena++ 
python configure.py\
			 --prob shock_tube\
			 --flux hlle\
			 --fluxcl hlle\
			 --nghost=3 \
			 -debug \
			 -hdf5 \
			 -cl \
			 --include /usr/include/hdf5/serial/ \
			 #--include /usr/include/hdf5/openmpi/ \
#			 --hdf5_path=/usr/lib/x86_64-linux-gnu/
