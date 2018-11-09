python3 configure.py --prob spiralarm --coord cartesian -hdf5 --cxx icc -mpi 
python3 configure.py --prob spiralarm --coord cartesian --eos isothermal -mpi -hdf5 --cxx icc --ccmd /nas/longleaf/apps-dogwood/hdf5/1.10.2/openmpi/bin/h5pcc --cflag="-DH5_HAVE_PARALLEL"
python3 configure.py --prob spiralarm --coord cartesian --eos isothermal --cxx icc -debug --cflag='-g -std=c++11'
python3 configure.py --prob spiralarm --coord cartesian --eos isothermal -mpi -hdf5 --cxx icc --ccmd /nas/longleaf/apps-dogwood/hdf5/1.10.2/openmpi/bin/h5pcc --cflag="-DH5_HAVE_PARALLEL" -fft --grav fft
