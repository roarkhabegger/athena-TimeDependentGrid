# Easy to edit configure script for athena++ 
python configure.py\
			 --prob blast\
			 --flux roe\
			 --include /usr/include/hdf5/serial/ \
			 --nghost=2 \
			 -cl \
			 -hdf5 \
#			 --hdf5_path=/usr/lib/x86_64-linux-gnu/
