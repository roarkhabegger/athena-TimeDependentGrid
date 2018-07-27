#!/usr/bin/python
import numpy as np
import sys
import os
import subprocess as sbp
from functools import reduce
# matplotlib
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from mpl_toolkits.axes_grid1 import make_axes_locatable
# scipy
from scipy.interpolate import interp2d 

import argparse
from argparse import RawTextHelpFormatter
import time as t
import re

# Import from correct directory
import socket as s
comp = s.gethostname()
if comp == 'thetis': 
    sys.path.insert(0,'/afs/cas.unc.edu/users/j/d/'
                  'jdupuy26/Johns_work/'
                  'misc_scripts/plotting_tools')
elif comp == 'debpad':
    sys.path.insert(0,'/home/jdupuy26/Johns_work/'
                      'Grad_Work/Research/codes/'
                      'my_athena++/vis/python/')
else: 
    print('[init]: Computer %s not recognized!' % comp)
import athena_read as ar


#===================================================================
#
#  Code: app_plot.py
#
#  Purpose: Reads in athena++ hdf5 dumps and plots
#           quantities. Designed to be a general plotter
#           that plots all quantities, and any coordinate system,
#           including collisionless variables. Very similar to 
#           plot_sims.py 
#
#  Keywords: python app_plot.py -h   
#
#  Usage: python app_plot.py quant   
#
#  WARNING: THIS MUST BE RUN FROM SIMULATION DIRECTORY 
#          
#  Author: John Dupuy 
#          UNC Chapel Hill
#  Date:    07/20/18 
#  Updated: 07/24/18 
#====================================================================

#============ FUNCTIONS ==============================

#\func natural_sort()
# does natural (i.e. human) sorting of a list 
# of strings
# cf. https://stackoverflow.com/questions/4836710/
#            does-python-have-a-built-in-function-for-string-natural-sort
def natural_sort(l):
    convert      = lambda text: int(text) if text.isdigit() else text.lower()
    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ]
    return sorted(l, key = alphanum_key)  

#------------ FUNCTIONS FOR READING DATA -------------#
#\func get_files()
# returns list of all files that are to be read   
def get_files(mydir='./'):
    
    #  directory is special since it contains other stuff
    files  = sorted([fname for fname in
                     os.listdir(mydir+'/')
                     if  '.athdf' in fname and '.xdmf' not in fname])
    return files  

# \func get_log():
# takes the log of data 
def get_log(data):
    data[np.where(data <= 0.0)] = 1e-80
    data = np.log10(data)
    return data 

#\func get_athinput():
# reads athinput and returns base, params
def get_athinput(cwd=-1):
    if cwd == -1:
        cwd = os.getcwd()
    athin = cwd + '/'+[fnm for fnm in os.listdir(os.getcwd()+'/')
                                        if fnm.startswith('athinput')][0]
    data = ar.athinput(athin) 
    return data 

#\func get_appquant():
# Given my quant name, converts it into app_quant name
def get_appquant(quant):
    app_quant = [] 
    # Conserved quantities 
    if quant == 'd':
        app_quant = ['dens']
    elif quant == 'M1':
        app_quant = ['mom1']
    elif quant == 'M2':
        app_quant = ['mom2']
    elif quant == 'M3':
        app_quant = ['mom3']
    elif quant == 'E':
        app_quant = ['Etot']
    elif quant == 'ie':
        app_quant = ['Eint']
    # Primitive quantities 
    elif quant == 'v1':
        app_quant = ['mom1','dens']
    elif quant == 'v2':
        app_quant = ['mom2','dens']
    elif quant == 'v3': 
        app_quant = ['mom3','dens']
    elif quant == 'eint':
        app_quant = ['Etot','mom1','mom2','mom3','dens']
    elif quant == 'ediff':
        app_quant = ['Eint','Etot','mom1','mom2','mom3','dens'] 
    else:
        print('[get_appquant]: quant %s not understood, exiting... ' %(quant))
        quit() 

    derived = False

    if len(app_quant) > 1:
        derived = True

    return app_quant, derived  

#\func get_quant 
# gets quantity 
def get_quant(file,quant,lev,derived=False,myquant='None'):
    
    #read in hdf5 file  
    # note: specifying level=0 restricts all data to coarsest level
    data = ar.athdf(file,quantities=quant,level=lev)

    if derived:
        if myquant == 'v1' or myquant == 'v2' or quant == 'v3':
            qdat = data[quant[0]]/data[quant[1]]
        elif myquant == 'eint':
            qdat = data[quant[0]] - 0.5 * (data[quant[1]]**2. + 
                                           data[quant[2]]**2. +
                                           data[quant[3]]**3.)/data[quant[4]]
        elif myquant == 'ediff':
            eint = data[quant[1]] - 0.5 * (data[quant[2]]**2. + 
                                           data[quant[3]]**2. +
                                           data[quant[4]]**2.)/data[quant[5]] 
            ie   = data[quant[0]]
            qdat = np.abs(eint - ie)

            rms  = np.sqrt(np.sum( (eint-ie)**2.)/eint.size)

            print('[get_quant]: time = %1.3f ediff_RMS = %13.5e' % (data['Time'], rms))
    else:
        qdat = data[quant[0]] 


    return data['Time'], data['x1v'], data['x2v'], data['x3v'], qdat   

#\func get_alldata
# given files, this will get all the data 
def get_alldata(files,quant,myfrms,dims,lev,**kwargs):
    mydir    = './'
    
    # read athinput
    params = get_athinput()  
   
    # get no. of points 
    nx3, nx2, nx1 = dims
    
    # get extent of grid
    mn1 = params['mesh']['x1min']
    mx1 = params['mesh']['x1max']
    mn2 = params['mesh']['x2min']
    mx2 = params['mesh']['x2max']
    mn3 = params['mesh']['x3min']
    mx3 = params['mesh']['x3max']

    # get no. of files
    nf     = len(myfrms) 

    # Make array to hold quant for all time
    quant_arr = np.zeros((nf, nx3, nx2, nx1))
    # Make time arrays 
    tarr = np.zeros(nf)

    app_quant, dflag = get_appquant(quant) 

    # Define processor index
    i = 0 
    for iff in myfrms:
        tarr[i],\
        x1, x2, x3,\
        quant_arr[i] = get_quant(files[iff], app_quant, lev, dflag, quant) 

        i += 1
   
    return tarr, x1, x2, x3, quant_arr

# \func get_fc()
# returns face centered grid & coordinate system 
# currently only necessary for 2d 
def get_fc(file,lev):
    # read hdf5 file
    all_things = ['dens']
    data = ar.athdf(file, quantities=all_things, level=lev) 
    
    return ( data['Coordinates'], data['x1f'], data['x2f'], data['x3f'], 
             data['dens'].shape )

def get_labels(quant,dim,log=False):

    # Define dictionary for quants
    lab_quants = {  # Fluid variables  
                  'E':'Energy density', 'ie':'Internal energy density', 
                  'eint': 'Derived internal energy density',
                  'ediff': '| IE - eint |', 
                  'd':'$\\rho$',
                  'n':'Column density', 'p':'Surface pressure',
                  'pie':'Surface pressure (from U.IE)', 
                  'T':'T',
                  'M':'M$_{tot}$',
                  'v1':'v$_1$','v2':'v$_2$','v3':'v$_3$',
                  'M1':'M$_1$','M2':'M$_2$','M3':'M$_3$',
                  'phi':'$\Phi_G$',
                  'cs':'c$_s$', 
                  'v':'v$_{tot}$',
                  's1':'$\Sigma_c$','s1c':'$M_c (R < 0.5 \; {\\rm [kpc])}/M_c$',
                  'jl':'$\lambda_J$','vlos':'$v_{\\rm los}$',
                    # Collisionless variables 
                  'dcl':'$\\rho_{\\rm cl}$', 
                  'v1cl':'v$_{1,\\rm cl}$','v2cl':'v$_{2,\\rm cl}$',
                  'v3cl':'v$_{3,\\rm cl}$','M1cl':'M$_{1,\\rm cl}$',
                  'vcl':'v$_{\\rm cl}$','Mcl':'M$_{\\rm cl}$', 
                  'M2cl':'M$_{2,\\rm cl}$','M3cl':'M$_{3,\\rm cl}$',
                  'P11':'P$_{11}$','P22':'P$_{22}$','P33':'P$_{33}$',
                  'p11':'P$_{ie,11}$','p22':'P$_{ie,22}$','p33':'P$_{ie,33}$',
                  'P12':'P$_{12}$','P23':'P$_{23}$','P13':'P$_{13}$',
                  'E11':'E$_{11}$','E22':'E$_{22}$','E33':'P$_{33}$',
                  'E12':'E$_{12}$','E23':'E$_{23}$','E13':'P$_{13}$',
                  'detP':'det(P$_{ij}$', 'detE':'det(E$_{ij}$',
                  'normP':'$|| P_{ij} ||_F$', 'normE':'$|| E_{ij} ||_F$'}

    # Define dictionary for units 
    units = ' [comp units]' 

    if dim == 1:
        xlabel = 'x'+units
        ylabel = ''
    if dim == 2 or dim == 3:
        xlabel = 'x'+units
        ylabel = 'y'+units

    cbar_l = lab_quants[quant]
    cbar_l += units 

    

    if log:
        cbar_l = 'log$_{10}$('+cbar_l+')'

    return xlabel, ylabel, cbar_l 


#\func get_args()
# this function parses CMD line args
def get_args():
    # Read in system arguments
    parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter) 
    parser.add_argument("quant",type=str,
                        help="Plotting options:\n"
                              " Any quant available in get_quant\n")
    parser.add_argument("--anim", dest="anim",action='store_true',
                        default=False,
                        help="Switch to do animation\n")
    parser.add_argument("--iani", dest="iani",nargs=2,required=False,
                        default=[0,0],type=int,
                        help="Animate from frame iani[0] to iani[1]\n")
    parser.add_argument("--qminmax", dest="qminmax",nargs='+',required=False,
                        default=-1,type=float,
                        help="Min/max value for imshow")
    parser.add_argument("--ifrm", dest="ifrm",type=int,default=[0],
                        nargs='+',
                        help="Frame of simulation to plot:\n"
                             "  0: tsim = 0 \n"
                             "  1: tsim = dt_dump\n"
                             "  2: tsim = 2*dt_dump\n"
                             "  .               \n"
                             "  .               \n"
                             "  .               \n"
                             "  n: tsim = n*dt_dump\n"
                             " Note: by entering multiple integers, will make\n"
                             " panel plots for each ifrm")
    parser.add_argument("--mnmx",dest="mnmx", type=float,nargs=2,
                        required=False,default=[-np.pi,np.pi],
                        help="Plotting range in x and y\n"
                             "Note: assumes x, y range are equivalent")
    parser.add_argument("--save", dest="save",action='store_true',
                        default=False,
                        help="Switch to save anim or figure")
    parser.add_argument("--log",dest="log",action='store_true',
                        default=False,
                        help="Switch to take log of images, default is False")
    parser.add_argument("--units",dest="units",type=int,required=False,
                        default=0, help="units: 0-comp,1-cgs,2-SI")  
    parser.add_argument("--grid", dest="grid",action='store_true',
                        default=False, help="Switch to make plot to show grid")
    parser.add_argument("--fmt", dest="fmt", default='eps',
                        type=str, help='format for saving graphics, default: eps') 
    parser.add_argument("--vvec",dest="vvec",action='store_true',
                        default=False, required=False,
                        help="Overplot velocity vectors\n")
    parser.add_argument("--noplot",dest="noplot",action='store_true',
                        default=False, required=False,
                        help="Switch to return only stitched together array\n"
                             "To be used if this file is imported from another file\n")
    parser.add_argument("--sliced1d",dest="sliced1d",action='store_true',
                        default=False, required=False,
                        help="Switch to take 1D slice of 2D array along DIAGONAL\n")
    parser.add_argument("--slicel1d",dest="slicel1d",action='store_true',
                        default=False, required=False,
                        help="Switch to take 1D slice of 2D array along a line\n")
    parser.add_argument("--col",dest="col",action='store_true',
                        default=False, required=False,
                        help="Sum 3d simuations along the z-axis, so they can be viewed as" 
                             "2d plots\n") 
    parser.add_argument("--slice",dest="slc",type=int, required=False,
                        default='-1', help="Slice 3d array into 2d along INT axis\n")
    parser.add_argument("--nocyl",dest="nocyl",action='store_true',
                        default=False, required=False,
                        help="Switch for plotting cylindrical simulations as R, \phi") 
    parser.add_argument("--lev", dest="lev", type=int, required=False, default=None,
                        help="Either prolongate or restrict data to lev for SMR, default\n"
                             "is to prolongate to finest level") 
    return parser.parse_args() 
    
 

#------------- MAIN FUNCTION ------------------------#
def main(args):
    ctable = 'magma'
    plt.rcParams['image.cmap'] = ctable

    # parsing arguments            
    quant = args.quant
    anim  = args.anim
    iani  = args.iani
    ifrm  = args.ifrm
    save  = args.save
    log   = args.log
    mnx, mxx = args.mnmx
    iunit    = args.units 
    qminmax  = args.qminmax
    grid     = args.grid
    fmt      = args.fmt
    vvec     = args.vvec
    noplot   = args.noplot
    sliced1d = args.sliced1d
    slicel1d = args.slicel1d
    col      = args.col
    slc      = args.slc   
    nocyl    = args.nocyl
    lev      = args.lev

    # get files 
    files = get_files()
    # get the coordinate system, face-centered grid & dimension of data 
    params = get_athinput()
    coordsys, x1f, x2f, x3f, dims = get_fc(files[0],lev)

    # Get qminmax flag 
    qflag = True if np.size(qminmax) > 1 else False
    # Get panel flag
    pflag = True if np.size(ifrm) > 1 else False 
    # Get mnmx flag
    mnmxflag = False if mxx == np.pi else True  
    # Get slice flag
    slice2d  = False if slc == -1 else True 
    # Set the coordinate system flag
    cyl      = True if coordsys == 'cylindrical' and not nocyl else False 

    if np.size(ifrm) == 1: ifrm = ifrm[0]
    
    # Change default iani values
    if iani[1] == 0:
        iani[1] = len(files) 

    # determine myframes
    if anim:
        myfrms = range(iani[0],iani[1])
    elif pflag:
        myfrms = ifrm
    else: 
        myfrms = [ifrm] 
    # get data
    tarr, x1, x2, x3, imgs = get_alldata(files,quant,myfrms,dims,lev)
    
    if log:
        imgs = get_log(imgs) 


    # Get dimensions of data
    nf, nx3, nx2, nx1 = imgs.shape 
    # Determine dimensional plotting flag
    flag3d, flag2d, flag1d = False, False, False  
    if nx3 > 1:
        flag3d = True
        dim    = 3
    elif (nx3 == 1) and (nx2 > 1):
        flag2d = True
        dim    = 2
    else:
        flag1d = True 
        dim    = 1

    # Change flags for slicing 
    if sliced1d or slicel1d:
        if flag2d:
            flag2d = False
            flag1d = True 
        if flag3d:
            flag3d = False
            flag1d = True 
    if col:
        flag3d = False
        flag2d = True
    if slice2d:
        flag3d = False
        flag2d = True 
    
    # Get face-centered grid if cylindrical  
    if flag2d and cyl: 
        x1f, x2f   = np.meshgrid(x1f, x2f, indexing='xy')
        x1cf, x2cf = x1f*np.cos(x2f), x1f*np.sin(x2f) 


    # Determine labels 
    xlab, ylab, clab = get_labels(quant,dim,log)

    # Now plot the data 
    if flag3d:
        from mpl_toolkits.mplot3d import Axes3D 
        fig = plt.figure()
        ax1 = fig.add_subplot(111, projection='3d')  
    else:
        fig = plt.figure(figsize=(7.5,5.5),facecolor='white') 
        ax1 = fig.add_subplot(111) 

    if flag1d: 
        if sliced1d:
            # Note this assumes nx1 = nx2 = nx3 
            if dim == 2:
                imgs = imgs[:,0,:,:]
                # Take slice 
                vals = []
                for j in range(len(imgs)):
                    vals.append(np.array([imgs[j,i,i] for i in range(len(x1))]))
                imgs = vals 
                # Get distance along diagonal 
                rad   = np.sqrt(x1**2. + x2**2.) 
                rad[x1 < 0] = -rad[x1 < 0] 
                x1    = rad.copy() 
            if dim == 3:
                # Take slice
                vals = []
                for j in range(len(imgs)):
                    vals.append(np.array([imgs[j,i,i,i] for i in range(len(x1))]))
                imgs = vals 
                # Get distance along diagonal
                rad = np.sqrt(x1**2. + x2**2. + x3**2.)
                rad[x1 < 0] = -rad[x1 < 0]
                x1 = rad.copy() 

        elif slicel1d:
            if dim == 2:
                imgs = imgs[:,0,:,:]
                # Take slice 
                vals = []
                for j in range(len(imgs)):
                    vals.append(np.array([imgs[j,nx2/2,i] for i in range(len(x1))]))
                imgs = vals 
            if dim == 3:
                # Take slice
                vals = []
                for j in range(len(imgs)):
                    vals.append(np.array([imgs[j,nx3/2,nx2/2,i] for i in range(len(x1))]))
                imgs = vals 
            
        # Get rid of unnecessary dimensions 
        else:
            imgs = imgs[:,0,0,:] 

        # Handle animation
        if anim:
            if qflag:
                qmin, qmax = qminmax[0], qminmax[1]
            else:
                qmin, qmax = np.min(imgs), np.max(imgs) 
            
            # Set labels 
            ax1.set_title('t = %1.2f' % (tarr[0]) ) 
            ax1.set_xlabel(xlab)
            ax1.set_ylabel(clab) 
            # Plot first frame 
            ax1.plot(x1, imgs[0], '.')
            
            def animate(ifrm):
                # Clear figure
                ax1.cla()
                # Set title & labels 
                ax1.set_title('t = %1.2f' % (tarr[ifrm]) )
                ax1.set_xlabel(xlab)
                ax1.set_ylabel(clab) 
                # Set xmin, xmax 
                if mnmxflag:
                    ax1.set_xlim(mnx, mxx)
                #ax1.set_ylim(qmin,qmax)
                # Now plot
                ax1.plot(x1, imgs[ifrm],'.')
                return
            ani = animation.FuncAnimation(fig, animate, range(len(myfrms)),
                                          repeat=False) 

        # Handle plotting a single frame 
        else:
            # plot 
            ax1.plot(x1, imgs[0],'.')  
            ax1.set_xlabel(xlab)
            ax1.set_ylabel(clab)
            ax1.set_title('t = %1.2f' % (tarr[0]) ) 
    
    elif flag2d:
        # Get extent of grid 
        if cyl:
            mnx1, mxx1, mnx2, mxx2 = ( np.min(x1cf), np.max(x1cf),
                                       np.min(x2cf), np.max(x2cf) )
        else: 
            mnx1, mxx1, mnx2, mxx2 = ( np.min(x1f), np.max(x1f),
                                       np.min(x2f), np.max(x2f) ) 

        # Determine colorbar 
        div = make_axes_locatable(ax1)
        cax = div.append_axes('right', '5%', '5%') 


        # Get rid of unnecessary dimensions 
        if col:
            imgs = np.sum(imgs, axis=1) 
        elif slice2d:
            imgs = imgs[:,slc,:,:] 
        else:
            imgs = imgs[:,0,:,:] 

        if qflag:
            qmin, qmax = qminmax[0], qminmax[1]
        else:
            qmin, qmax = np.min(imgs), np.max(imgs) 
    
        # Handle animation
        if anim: 
            ax1.set_title('t = %1.2f' %(tarr[0]) )
            if cyl:
                im   = ax1.pcolorfast(x1cf,x2cf, imgs[0], vmin=qmin,vmax=qmax)
                ax1.set_aspect('equal')
            else: 
                im   = ax1.imshow(imgs[0], extent=(mnx1, mxx1, mnx2, mxx2),
                                           vmin=qmin,vmax=qmax, origin='lower',
                                           interpolation='None')
            im.set_rasterized(True) 
            cbar = fig.colorbar(im,label=clab,cax=cax) 

            def animate(ifrm):
                # Clear figure 
                ax1.cla()
                cax.cla()
                
                # Set title 
                ax1.set_title('t = %1.2f' %(tarr[ifrm]) )
                # Plot 
                if cyl:
                    im   = ax1.pcolorfast(x1cf,x2cf, imgs[ifrm], vmin=qmin,vmax=qmax)
                else:
                    im   = ax1.imshow(imgs[ifrm], extent=(mnx1, mxx1, mnx2, mxx2),
                                                  vmin=qmin,vmax=qmax, origin='lower',
                                                  interpolation='None')
                im.set_rasterized(True) 
                # Set labels for x,y
                ax1.set_xlabel(xlab)
                ax1.set_ylabel(ylab)
                # Set xmin, xmax
                if mnmxflag:
                    ax1.set_xlim(mnx,mxx)
                    ax1.set_ylim(mnx,mxx)
                else:
                    ax1.set_xlim(mnx1,mxx1)
                    ax1.set_ylim(mnx2,mxx2) 

                # Set aspect ratio
                ax1.set_aspect('equal')
                # Set colorbar
                cbar = fig.colorbar(im,label=clab, cax=cax)
                return   

            ani = animation.FuncAnimation(fig, animate, range(len(myfrms)),
                                          repeat=False)
        # Handle a single frame 
        else:
            print(np.mean(imgs[0]))  
            ax1.set_xlabel(xlab)
            ax1.set_ylabel(ylab)

            ax1.set_title('t = %1.2f' %(tarr[0]) )
            if cyl:
                im   = ax1.pcolorfast(x1cf,x2cf, imgs[0], vmin=qmin,vmax=qmax)
                ax1.set_aspect('equal')
            else:
                im   = ax1.imshow(imgs[0], extent=(mnx1, mxx1, mnx2, mxx2),
                                           vmin=qmin,vmax=qmax, origin='lower',
                                           interpolation='None')
            im.set_rasterized(True) 
            cbar = fig.colorbar(im,label=clab,cax=cax) 

    # Handle 3d plotting (using mayavi) 
    elif flag3d: 
        # Get extent of grid 
        mnx1, mxx1, mnx2, mxx2, mnx3, mxx3 = ( np.min(x1), np.max(x1),
                                               np.min(x2), np.max(x2),
                                               np.min(x3), np.max(x3) ) 
        dx1, dx2, dx3 = x1[1]-x1[0], x2[1]-x2[0], x3[1]-x3[0]

        mnx1 -= 0.5*dx1; mxx1 += 0.5*dx1
        mnx2 -= 0.5*dx2; mxx2 += 0.5*dx2
        mnx3 -= 0.5*dx3; mxx3 += 0.5*dx3

        x1, x2, x3 = np.meshgrid(x1,x2,x3) 
        
        ax1.scatter(x1,x2,x3, c=imgs[0].ravel(),cmap=plt.hot()) 
    else: 
        print("[main]: Unsure what to plot  :( ")  
            

    if save:
        #mydir  = '/srv/analysis/jdupuy26/figures/'
        mydir = os.getcwd()+'/'
        # Create file name (have to make sure it is unique for each sim to avoid overwrites)  
        #myname = os.path.basename(os.path.dirname(os.path.realpath('bgsbu.log')))
        #myname = os.getcwd().split('longevity_study/',1)[1].replace('/','_') 
        myname = '' 
        if anim:
            print("[main]: Saving animation...")
            ani.save(mydir+myname+'_'+base+'_'+quant+'.gif',fps=15.
                     ,writer='imagemagick')

        else:
            print("[main]: Saving frame...")
            plt.savefig(mydir+myname+'_'+base+'_'+quant+str(ifrm)+'.'+fmt, format=fmt,bbox_inches='tight')
    else:
        plt.show() 

        
if __name__ == '__main__':
   # If this file is called from cmd line 
   args = get_args() 
   main(args) 

