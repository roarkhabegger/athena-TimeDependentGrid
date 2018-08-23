# Easy to edit configure script for athena++ 
python configure.py\
			 --prob blast\
			 --flux roe\
			 --nghost=3 \
			 -debug \
			 -hdf5 \
			 --include /usr/include/hdf5/serial/ \
			 --fluxcl hlle\
			 -cl \
			 #--include /usr/include/hdf5/openmpi/ \
#			 --hdf5_path=/usr/lib/x86_64-linux-gnu/
