from landlab import RasterModelGrid # required for creating grids
from landlab import CLOSED_BOUNDARY, FIXED_VALUE_BOUNDARY # reuqired for bounday condition
from landlab.components import LinearDiffuser, FlowRouter,DrainageDensity # Linear Diffuser and Flowrouter are two different component used for Linear approach for tectonic revolution, and surface water routing processes
from landlab.components import FastscapeEroder
from landlab.plot import imshow_grid # for plotting
from matplotlib import pyplot as plt
import numpy as np
import math



############################################################################
################# Calculate TMR upon the generated DEM#######################
############### slop calculation by components ###########
def CalandConv_slope(rmg,n_node_col,n_node_row,grid_space,file_path_out,slopefilename,sfpath):
    
    rmg.at_node.keys()  # see the existing fields at nodes
    rmg['node']['topographic__elevation']
    #slop_array=rmg['node']['topographic__steepest_slope'] # not used in this case, coz we need slope x and slope y
    #core_slop=slop_array[rmg.core_nodes]
    
    # calculate the slope at x and y, using the function .calc_slope_at_node
    (array_of_magnitude1,(array_of_slope_x_radians1, array_of_slope_y_radians1))=rmg.calc_slope_at_node(elevs='topographic__elevation', method='patch_mean', ignore_closed_nodes=True, return_components=True)
    slope_x1=array_of_magnitude1*array_of_slope_x_radians1
    slope_y1=array_of_magnitude1*array_of_slope_y_radians1
    #landlab list # check the available component
    slope_x_active=slope_x1[rmg.core_nodes]
    slope_y_active=slope_y1[rmg.core_nodes]
    
    rmg.number_of_core_nodes
    
    
    ######################################################################
    ################ write pfb file, and plot for visual checking #################
    #####################################################################
    ################# write parflow required files- soil ################
    ################ define the value 
    sx = slope_x_active  #slope x
    sy = slope_y_active #slope y
    # define the dimension of the file
    nrows = n_node_col-2
    ncols = n_node_row-2
    sx.shape=(nrows,ncols)
    np.flipud(sx) ## flip it from bottom to top
    #sx=sx[1:nrows-1,1:ncols-1]
    sy.shape=(nrows,ncols)
    np.flipud(sy)
    #sy=sy[1:nrows-1,1:ncols-1]
    # nlayers= 1 #defined in each computation as NZ before calling the function
    cellsize = grid_space
    
    ################################# sx, write pfb file and plot it into png for visual checking
    slopename_full= '%s.sx.pfb' %(slopefilename)
    fn = '%s/%s' % (file_path_out,slopename_full)

    import pfbwriting as pfbwtn
    #pfbwtn.pfb_writer_onelayer(0, 0, 0, ncols-2, nrows-2, 1, cellsize, cellsize, 0, sx, 0, 0, 0, ncols-2, nrows-2,  1, 0, 0, 0, 1, fn)
    pfbwtn.pfb_writer_onelayer(0, 0, 0, ncols, nrows, 1, cellsize, cellsize, 0, sx, 0, 0, 0, ncols, nrows,  1, 0, 0, 0, 1, fn)
    print 'DEMfil.sx.pfb: created' 
    
    # plot
    plt.figure(figsize=(10, 7))
    plt.pcolor(sx,cmap='jet')
    fname ='%s.sx' %(slopefilename)
    #sfpath = '/home/chaochen/payette/chaochen/CZO/Chao/Demonstration_example/landlab/preprocessing'
    plt.colorbar()
    plt.savefig('%s/%s.png'% (sfpath,fname))
    plt.close()
    
    #################################### sy, write pfb file and plot it into png for visual checking
    slopename_full= '%s.sy.pfb' %(slopefilename)   
    fn = '%s/%s' % (file_path_out,slopename_full)
    #pfbwtn.pfb_writer_onelayer(0, 0, 0, ncols-2, nrows-2, 1, cellsize, cellsize, 0, sy, 0, 0, 0, ncols-2, nrows-2, 1, 0, 0, 0, 1, fn)
    pfbwtn.pfb_writer_onelayer(0, 0, 0, ncols, nrows, 1, cellsize, cellsize, 0, sy, 0, 0, 0, ncols, nrows, 1, 0, 0, 0, 1, fn)
    print 'DEMfil.sy.pfb: created' 
    
    plt.figure(figsize=(10, 7))
    plt.pcolor(sy<0.000001,cmap='jet')
    fname ='%s.sy' %(slopefilename)
    #sfpath = '/home/chaochen/payette/chaochen/CZO/Chao/Demonstration_example/landlab/preprocessing'
    plt.colorbar()
    plt.savefig('%s/%s.png'% (sfpath,fname))
    plt.close()


###############################################################################
###############  ###########
## DrainageDensity: Calculate drainage density from topography
def CalandConv_DrainageDensity(rmg,DD_thrshld,filename):
    rmg.at_node['topographic__elevation'][rmg.core_nodes]
    channels = np.array(rmg.at_node['drainage_area'] > DD_thrshld, dtype=np.uint8)
    dd = DrainageDensity(rmg,channel__mask=channels)
    mean_drainage_density = dd.calc_drainage_density()
    print 'the mean drainage density of %s is: %f' % (filename,mean_drainage_density)
    
    return mean_drainage_density
#    np.isclose(mean_drainage_density,0.0297364)

###############################################################################
# Claculatet he slope versus area correlation
def CalandConv_SlopeVsArea(rmg,plotsave_path,filename):
    from matplotlib.pyplot import loglog,figure, show, plot, xlabel, ylabel, title 
    steepest_slope=rmg.at_node['topographic__steepest_slope']
    drainage_area=rmg.at_node['drainage_area']
      
    x = steepest_slope[rmg.core_nodes]
    print 'the largest slope is: %f' % (x.max())
    print 'the smallest slope is: %f' %(x.min())
    y = drainage_area[rmg.core_nodes]
    #y = np.log10(y)
    
    
    #plt.figure(figsize=(10, 7))
    #plt.scatter(y,x)
    figure('final slope-area plot')
    loglog(y,x,'.')
    xlabel('Drainage (km**2)')
    ylabel('Local slope')
    title('Slope-Area plot for whole landcape')
    
    fname='%s_SlopeVsArea.png' %(filename)
    plt.savefig('%s/%s'% (plotsave_path,fname))
    plt.close()
###############################################################################
# calculate the curvature,which is the flux divergence when D =1
def CalandConv_curvature(rmg,n_node_row,n_node_col,plotsave_path,filename):
    z= rmg.at_node['topographic__elevation'] # extract the elevation
    slop_link=rmg.calc_grad_at_link(z)       # calcalute the gradience at link, which is the input of the flux divegent function
    curvature=rmg.calc_flux_div_at_node(slop_link)     # calculate the flux divergent, when no -D*slope_link defined, which means D =1, and this is curvature
    curvature_core=curvature[rmg.core_nodes]
    # make it into a matrix for plotting
    curvature_core.shape=(n_node_row-2,n_node_col-2)
    np.flipud(curvature_core) 
    curvature_cut = curvature_core[1:n_node_row-3,1:n_node_col-3]
    
    
    
    #plot
    plt.figure(figsize=(10, 7))
    plt.pcolor(curvature_cut>0.0025,cmap='jet')
    plt.colorbar()
    plt.savefig('%s/%s.png'% (plotsave_path,filename))
    plt.close()
    
    #output the streams
    stream_network=np.zeros((96,96),dtype=int)
    stream_location=curvature_cut>0.0025
    stream_network[stream_location]=1
#    np.matrix(stream_network)
 #   [row,col]=DEM_WS1_valid.shape
    filename_ST=('streamflow_network%s.txt' %filename)
    with file(filename_ST, 'w') as ST:
        np.savetxt(ST,stream_network,fmt='%i')
        #ST.write('%i %i 1\n' % (col,row))
        #for line in DEM_WS1_valid:
            #np.savetxt(DEMfile, line, fmt='%.2f')
            
    return curvature 
###############################################################################
# calculate the TMR based on linear correaltion with two parameters
def CalandConv_TMR(rmg,n_node_row,n_node_col,curvature, a,b,plotsave_path,filename):
    TMR = a*curvature+b
    TMR_core=TMR[rmg.core_nodes]
    TMR_core.shape=(n_node_row-2,n_node_col-2)
    np.flipud(TMR_core)   
    TMR_cut = TMR_core[1:n_node_row-3,1:n_node_col-3]
    
    #plot
    plt.figure(figsize=(10, 7))
    plt.pcolor(TMR_cut,cmap='jet')
    plt.colorbar()
    plt.savefig('%s/%s.png'% (plotsave_path,filename))
    #plt.close()    
    return TMR_cut
    
