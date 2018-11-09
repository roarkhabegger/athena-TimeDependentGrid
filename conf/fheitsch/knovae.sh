python3 configure.py --prob knovae --coord spherical_polar -hdf5 --cxx icc -de
python3 configure.py --prob knovae --coord spherical_polar -hdf5 --cxx g++ -de
python3 configure.py --prob knovae --coord spherical_polar -hdf5 --cxx g++ -de --hdf5_path /usr/local/hdf5
python3 configure.py --prob knovae --coord spherical_polar -mpi --ns 2 -hdf5 --cxx icc --ccmd /nas/longleaf/apps-dogwood/hdf5/1.10.2/openmpi/bin/h5pcc --cflag="-DH5_HAVE_PARALLEL"
python3 configure.py --prob knovae --coord cartesian -mpi --ns 2 -hdf5 --cxx icc --ccmd /nas/longleaf/apps-dogwood/hdf5/1.10.2/openmpi/bin/h5pcc --cflag="-DH5_HAVE_PARALLEL"
python3 configure.py --prob knovae --coord spherical_polar --ns 2 --cxx icc -debug --cflag="-std=c++11"
