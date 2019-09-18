#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 23 22:08:48 2018

@author: nicgaspar
"""

## Linear diffusion and channels on a raster

## Import what is needed
from landlab import RasterModelGrid # required for creating grids
from landlab import CLOSED_BOUNDARY, FIXED_VALUE_BOUNDARY # reuqired for bounday condition
from landlab.components import LinearDiffuser, FlowRouter,DrainageDensity # Linear Diffuser and Flowrouter are two different component used for Linear approach for tectonic revolution, and surface water routing processes
from landlab.components import FastscapeEroder   #Compute fluvial erosion using stream power theory (“fastscape” algorithm)
from landlab.plot import imshow_grid # for plotting
from matplotlib import pyplot as plt
import numpy as np
from numpy import newaxis
import time

#counting time
#t=time.time()

for seed_n in range(71, 72):

    n_node_col=100
    n_node_row=100
    grid_space=20
    
    ## Make a grid that is 100 by 100 with dx=dy=100. m
    rmg1 = RasterModelGrid((n_node_col,n_node_col),grid_space)
    rmg2 = RasterModelGrid((n_node_col,n_node_col),grid_space)
    #slps=RasterModelGrid((n_node_col,n_node_col),grid_space)
    ## Add elevation field to the grid.
    # Note that you need to create the fields that a components takes as inputs before instantiating a component - though you can put values into the arrays later if you need to; by CC from tutorial 
    # create the fiels amd fill and field with ones; define the element first, i.e., 'node', 'link','cell', the second id the name to give the field
    z1 = rmg1.add_ones('node','topographic__elevation')
    z2 = rmg2.add_ones('node','topographic__elevation')
    ## Add noise to initial landscape
    # np.random.rand() gives a uniform distribution over [0,1] with the give shape
    #z1[rmg1.core_nodes] += np.random.rand(len(rmg1.core_nodes))
    #z2[rmg2.core_nodes] += np.random.rand(len(rmg2.core_nodes))
    
    np.random.seed(seed_n)
    z1[rmg1.core_nodes] += np.random.rand(len(rmg1.core_nodes))
    z2[rmg2.core_nodes] += np.random.rand(len(rmg2.core_nodes))
    
    ## Set boundary conditions
    ## Close all perimeter nodes first
    rmg1.set_status_at_node_on_edges(right=CLOSED_BOUNDARY, top=CLOSED_BOUNDARY, \
                                  left=CLOSED_BOUNDARY, bottom=CLOSED_BOUNDARY)
    rmg2.set_status_at_node_on_edges(right=CLOSED_BOUNDARY, top=CLOSED_BOUNDARY, \
                                  left=CLOSED_BOUNDARY, bottom=CLOSED_BOUNDARY)
    ## Now open the corner one
    rmg1.status_at_node[0] = FIXED_VALUE_BOUNDARY
    rmg2.status_at_node[0] = FIXED_VALUE_BOUNDARY
    
    ###############################################################################
    ################### Waterhsed 1  ##############################################
    #### Instantiate process components
    ld1 = LinearDiffuser(rmg1, linear_diffusivity=0.3)                             ### parameters defined
    fr1 = FlowRouter(rmg1, method='D8')
    #slp1 = FlowRouter(slps,method='D8')
    fse1 = FastscapeEroder(rmg1, K_sp = 1e-5, m_sp=0.5, n_sp=1.)                   ### parameters defined                  
    
    ## Set some variables
    rock_up_rate = 1e-4 # m/yr # uplift rate                                       ### parameter defined
    dt = 1000 # yr
    rock_up_len = dt*rock_up_rate # m
    
    ## Time loop where evolution happens
    for i in range(50000):
        z1[rmg1.core_nodes] += rock_up_len #uplift only the core nodes
        ld1.run_one_step(dt) #linear diffusion happens.
        fr1.run_one_step() #flow routing happens, time step not needed, as this is not time sensitive
        fse1.run_one_step(dt) #fluvial incision happens
        ## optional print statement
    #    print('i', i)
        
    ## Plotting the topography
    plt.figure(1)
    imshow_grid(rmg1, 'topographic__elevation', grid_units = ['m','m'],
                     var_name='Elevation (m)')
    plt.savefig('diffusive_%i.png' %seed_n)
    plt.close()
    #need to run for about 4000 time steps, or 4,000,000 years to reach SS
    
    
    ################### Waterhsed 2  ##############################################
    ###### Instantiate process components
    ld2 = LinearDiffuser(rmg2, linear_diffusivity=0.1)
    fr2 = FlowRouter(rmg2, method='D8')
    fse2 = FastscapeEroder(rmg2, K_sp = 1e-5, m_sp=0.5, n_sp=1.)
    ## Time loop where evolution happens
    for i in range(50000):
        z2[rmg2.core_nodes] += rock_up_len #uplift only the core nodes
        ld2.run_one_step(dt) #linear diffusion happens.
        fr2.run_one_step() #flow routing happens, time step not needed, as this is not time sensitive
        fse2.run_one_step(dt) #fluvial incision happens
    #    slps['node']['topographic__elevation']=rmg1['node']['topographic__elevation']
    #    slp1.run_one_step()
        
        ## optional print statement
    #    print('i', i)
        
    ## Plotting the topography
    plt.figure(2)
    imshow_grid(rmg2, 'topographic__elevation', grid_units = ['m','m'],
                     var_name='Elevation (m)')
    plt.savefig('fluvial_%i.png' %seed_n)
    plt.close()
    ##############################################################################
    ## extract slope component by x and y direction, and write into pfb files
    
    file_path_out = '/home/chaochen/payette/chaochen/CZO/Chao/Demonstration_example/landlab_R2/parflow_input' # the file location defined for DEM.sx.pfb 
    plot_saving_path='/home/chaochen/payette/chaochen/CZO/Chao/Demonstration_example/landlab_R2/preprocessing'
    soil_indi_path_out = '/home/chaochen/payette/chaochen/CZO/Chao/Demonstration_example/landlab_R2/parflow_input/geo_indi'
    #file_path_out = '/home/cchen/payette/chaochen/CZO/Chao/Demonstration_example/landlab/parflow_input'
    #plot_saving_path='/home/cchen/payette/chaochen/CZO/Chao/Demonstration_example/landlab/preprocessing'
    
    
    ################################ write into text file, DEM
    ##organize the data first
    DEM_WS1=rmg1.at_node['topographic__elevation'][rmg1.core_nodes]
    DEM_WS2=rmg2.at_node['topographic__elevation'][rmg2.core_nodes]
    
    nrows_valid = n_node_col-2
    ncols_valid = n_node_row-2
    #nrows_valid = n_node_col
    #ncols_valid = n_node_row
    DEM_WS1.shape=(nrows_valid,ncols_valid)
    np.flipud(DEM_WS1) ## flip it from bottom to top
    #DEM_WS1_valid=DEM_WS1[1:nrows_valid-1,1:ncols_valid-1]
    DEM_WS1_valid=DEM_WS1
    
    DEM_WS2.shape=(nrows_valid,ncols_valid)
    np.flipud(DEM_WS2) ## flip it from bottom to top
    #DEM_WS2_valid=DEM_WS2[1:nrows_valid-1,1:ncols_valid-1]
    DEM_WS2_valid=DEM_WS2
    
    #write to file
    np.matrix(DEM_WS1_valid)
    [row,col]=DEM_WS1_valid.shape
    filename_WS1=('DEMfil_WS1_%i.dem.sa' %(seed_n))
    filename_matrix_WS1=('dem_CC1.txt')

    #    DEMfile.write('ncols %i\n' %col)
    #    DEMfile.write('nrows %i\n' %row)
    #    DEMfile.write('xllcorner 0.00\n')
    #    DEMfile.write('yllcorner 0.00\n')
    #    DEMfile.write('cellsize %i\n' %grid_space)
#    with file(filename_matrix_WS1, 'w') as DEMfile:
#        np.savetxt(DEMfile,DEM_WS1_valid,fmt='%.2f')        
    with file(filename_WS1, 'w') as DEMfile:
        DEMfile.write('%i %i 1\n' % (col,row))
        for line in DEM_WS1_valid:
            np.savetxt(DEMfile, line, fmt='%.2f')
    
    np.matrix(DEM_WS2_valid)
    [row,col]=DEM_WS2_valid.shape
    
    filename_WS2=('DEMfil_WS2_%i.dem.sa' %(seed_n))
    filename_matrix_WS2=('dem_CC2.txt')

    #    DEMfile.write('ncols %i\n' %col)
    #    DEMfile.write('nrows %i\n' %row)
    #    DEMfile.write('xllcorner 0.00\n')
    #    DEMfile.write('yllcorner 0.00\n')
    #    DEMfile.write('cellsize %i\n' %grid_space)
#    with file(filename_matrix_WS2, 'w') as DEMfile:
#        np.savetxt(DEMfile,DEM_WS2_valid,fmt='%.2f')
    with file(filename_WS2, 'w') as DEMfile:
        DEMfile.write('%i %i 1\n' % (col,row))
        for line in DEM_WS2_valid:
            np.savetxt(DEMfile, line, fmt='%.2f')
    
    
    
    
    # this way of calculating slope in two direction may not be valid, as the elevation just calculated has not been 
            # pit-filled. That being said there may have some puddles in the DEM
    from CalandConv import CalandConv_slope as CF_slope
    #rmg1_3D=rmg1[:,:,newaxis]
    #rmg2_3D=rmg2[:,:,newaxis]
    
    CF_slope(rmg1,n_node_col,n_node_row,grid_space,file_path_out,'DEMfil_WS1',plot_saving_path) # 'DEMfil' is the defined for the name of DEMfil.sx.pfb,DEMfil.sy.pfb
    CF_slope(rmg2,n_node_col,n_node_row,grid_space,file_path_out,'DEMfil_WS2',plot_saving_path) # 'DEMfil' is the defined for the name of DEMfil.sx.pfb,DEMfil.sy.pfb
    
    from CalandConv import CalandConv_DrainageDensity as CF_dd
    dd_thrshld=10000
    DrainageDensity_WS1 = CF_dd(rmg1,dd_thrshld,'DEMfil_WS1')
    DrainageDensity_WS2 = CF_dd(rmg2,dd_thrshld,'DEMfil_WS2')
    
    from CalandConv import CalandConv_SlopeVsArea as CF_sa
    CF_sa(rmg1,plot_saving_path,'DEMfil_WS1')
    CF_sa(rmg2,plot_saving_path,'DEMfil_WS2')
    
    from CalandConv import CalandConv_curvature as CF_cur
    fileName_curvature1=('Curvature_WS1_%i' %seed_n)
    fileName_curvature2=('Curvature_WS2_%i' %seed_n)
    Curvature_WS1 = CF_cur(rmg1,n_node_row,n_node_col,plot_saving_path,fileName_curvature1)  # the output is a vector
    Curvature_WS2 = CF_cur(rmg2,n_node_row,n_node_col,plot_saving_path,fileName_curvature2)  # the output is a vector
    
    ###############################################################################
    # calculat the Total Mobilr Regolith (TMR) using Nich patton's work
    # TMR=a*curvature+b
    # the coefficients for the linear corresltion are subjected to change
    a = 16.86
    b = 0.57
    
    # make the array into matrix and cut only the core area
    from CalandConv import CalandConv_TMR as CF_TMR
    TMR_WS1 = CF_TMR(rmg1,n_node_row,n_node_col,Curvature_WS1, a,b,plot_saving_path,'TMR_WS1')
    TMR_WS2 = CF_TMR(rmg2,n_node_row,n_node_col,Curvature_WS2, a,b,plot_saving_path,'TMR_WS2')
    
    ################################################################################
    # calculate soil indicator file
    from Soil_parameterInParfow_tclwritting import soil_indi_cal
    Soil_indi_WS1=('Soil_indi_WS1_%i' %seed_n)
    Soil_indi_WS2=('Soil_indi_WS2_%i' %seed_n)
    soil_indi_cal(TMR_WS1,soil_indi_path_out,Soil_indi_WS1,n_node_row,n_node_col,grid_space)
    soil_indi_cal(TMR_WS2,soil_indi_path_out,Soil_indi_WS2,n_node_row,n_node_col,grid_space)


## load filled dem from the output of GIS
#
#filleddem1=np.loadtxt('filleddem_cc1.txt', skiprows=6)
#with file(filename_WS1, 'w') as DEMfile:
#    DEMfile.write('%i %i 1\n' % (col,row))
#    for line in filleddem1:
#        np.savetxt(DEMfile, line, fmt='%.2f')