##############################
# this pyhton file is tended to generate indicatorfile, which is a 3-D ssubsurface media index file
# this 3D index, especially the vertical discretization, is depending the veritical soil discretization
# however, this vertical soil discretization is defined in the parflow tcl script
# so before or after excute this file to genrate the soil indicator filw, alwasys check the vetical discretization
# to see if it is the same as in the Parflow tcl script
import numpy as np
import pandas as pd

# using the developed TMR for soil
# define the vertical delineation for the subsurface, 1 for soil, 2 for saprolit, 3 for bedrock
# define soil goes as deep as 10 layers
def soil_indi_cal(TMR,file_path_out,WS_name,n_node_row,n_node_col,grid_space):
#TMR=TMR_WS1
#TMR=TMR_WS2

    n_layer_soil = 10
    
    # define saprolite goes as deep as 15 layers, for the convenience of comuputation convergence
    n_layer_saprolite = 5
    depth_saprolite = 30
    # define bedrock as the 16 layer
    n_layer_bedrock = 1
    depth_bedrock = 40
    ###
    # define the dzScale.value, which should be exactly the same as in the Parflow script
    # make sure this is the same as defined in the tcl script
    # totoally 16 layers as defined here
    dzScalr_Valur_list_set = [2,
                          2,
                          1,
                          1,
                          1,
                          0.5,
                          0.08,
                          0.06,
                          0.06,
                          0.06,
                          0.06,
                          0.06,
                          0.05,
                          0.04,
                          0.02,
                          0.01]  
    dzScalr_Valur_list=dzScalr_Valur_list_set[::-1]
    
    thickness_soil_list = []# creat an empty list for the weight
    Indicator_matrix = np.zeros((n_node_row-4,n_node_col-4,(n_layer_soil+n_layer_saprolite+n_layer_bedrock)),dtype=np.int)
    temp_value=0
    ############## detemine the thickness of each soil layer
    # distribute the thickness according to weight
    for i in range(n_layer_soil):   #10 layer here
        soil_index_sum = sum(dzScalr_Valur_list[0:n_layer_soil]) #0:9 here
        temp_value += dzScalr_Valur_list[i]/soil_index_sum*TMR.max()
        thickness_soil_list.append(temp_value)
    
    thickness_saprolite_list = []
    temp_value=0
    ############### detemine the thickness of each saprolite layer
    #for i in range (n_layer_soil,n_layer_soil+n_layer_saprolite):  #10:15
    #    saprolite_index_sum = sum(dzScalr_Valur_list[n_layer_soil:n_layer_soil+n_layer_saprolite]) #10:14
    #    temp_value += dzScalr_Valur_list[i]/saprolite_index_sum*(depth_saprolite-TMR.max()) # here refers to the thickness below the deepest soil layer
    #    thickness_saprolite_list.append(temp_value) 
    #
    ############### detemine the thickness of bedrock layer
    #thickness_bedrock = depth_bedrock-depth_saprolite
    
    
    ############ construct a indicator matrix
    
    str_soil=6
    end_soil=16
    for i in range(n_node_col-4):
        for j in range(n_node_col-4):
            Indicator_matrix[i,j,1:5] = 2
            Indicator_matrix[i,j,0] = 3
            for layer_n in range (str_soil,end_soil):
                if TMR[i,j] >= (thickness_soil_list[layer_n-6]):
                    Indicator_matrix[i,j,(end_soil-layer_n+5)] = 1 # layer 1,
                else:
                    Indicator_matrix[i,j,(end_soil-layer_n+5)] = 2
    
    ###############################################################################
    ########################### write pfb file ####################################
    # define the output path
    #file_path_out = '/home/chaochen/payette/chaochen/CZO/Chao/Demonstration_example/landlab_R2/parflow_input/geo_indi'
#    fn = '%s/%s.indi.pfb' % (file_path_out,WS_name)
#    from pfbwriting import pfb_writer_multilayer as pfb_wrt_multi
#    pfb_wrt_multi(0, 0, 0, n_node_row-4, n_node_col-4, 16, grid_space, grid_space, 5, Indicator_matrix, 0, 0, 0, n_node_row-4, n_node_col-4, 16, 0, 0, 0, 1,fn)
#    print 'DEMfil_WS.indi.pfb: created' 
    
    
    # make it into 100 by 100, by extedning the edge values
    x_left  =Indicator_matrix[0,:,:]
    x_right =Indicator_matrix[95,:,:]

    Indicator_matrix_copy=Indicator_matrix
    temp1=np.insert(Indicator_matrix_copy,0,x_left,axis=0)
    temp2=np.insert(temp1,96,x_right,axis=0)

    y_lower =temp2[:,0,:]
    y_upper =temp2[:,95,:]
    temp3=np.insert(temp2,0,y_lower,axis=1)
    temp4=np.insert(temp3,96,y_upper,axis=1)
    
   
    
    
    fn = '%s/%s.indi.pfb' % (file_path_out,WS_name)
    from pfbwriting import pfb_writer_multilayer as pfb_wrt_multi
    pfb_wrt_multi(0, 0, 0, n_node_row-2, n_node_col-2, 16, grid_space, grid_space, 5, temp4, 0, 0, 0, n_node_row-2, n_node_col-2, 16, 0, 0, 0, 1,fn)
    print 'DEMfil_WS.indi.pfb: created'


