import Gadget as G
import numpy  as np

def reading_snap(_path,_snap,_parttype):
    ids_  = G.read_block(_path+_snap,"ID  ",parttype=_parttype)
    pos_  = G.read_block(_path+_snap,"POS ",parttype=_parttype)
    mass_ = G.read_block(_path+_snap,"MASS",parttype=_parttype)
    rho_  = G.read_block(_path+_snap,"RHO ",parttype=_parttype)
    vel_  = G.read_block(_path+_snap,"VEL ",parttype=_parttype)
    bfld_ = G.read_block(_path+_snap,"BFLD",parttype=_parttype)
    vrms_ = G.read_block(_path+_snap,"VRMS",parttype=_parttype)
    hsml_ = G.read_block(_path+_snap,"HSML",parttype=_parttype)
    inte_ = G.read_block(_path+_snap,"U   ",parttype=_parttype)
    
    x_     =[row[0] for row in pos_]
    y_     =[row[1] for row in pos_]
    z_     =[row[2] for row in pos_]
    xbfld_ =[row[0] for row in bfld_]
    ybfld_ =[row[1] for row in bfld_]
    zbfld_ =[row[2] for row in bfld_]
    vx_    =[row[0] for row in vel_]
    vy_    =[row[1] for row in vel_]
    vz_    =[row[2] for row in vel_]
 
    return ids_,x_,y_,z_,mass_,rho_,vx_,vy_,vz_,xbfld_,ybfld_,zbfld_,vrms_,hsml_,inte_

def sorting_by_ids(_id,_x,_y,_z,_mass,_rho,_vx,_vy,_vz,_xbfld,_ybfld,_zbfld,_vrms,_hsml,_inte):
    id_,x_,y_,z_,mass_,rho_,vx_,vy_,vz_,xbfld_,ybfld_,zbfld_,vrms_,hsml_,inte_ = zip(*sorted(zip(_id,_x,_y,_z,_mass,_rho,_vx,_vy,_vz,_xbfld,_ybfld,_zbfld,_vrms,_hsml,_inte)))
    return np.array(id_),np.array(x_),np.array(y_),np.array(z_),np.array(mass_),np.array(rho_),np.array(vx_),np.array(vy_),np.array(vz_),np.array(xbfld_),np.array(ybfld_),np.array(zbfld_),np.array(vrms_),np.array(hsml_),np.array(inte_)
    #path=path_save+'IMSHOW_VRMS/z_ABS3/IM_VRMSxz_'+sindex[sn]+'.jpg'
    #path=path_save+'IMSHOW_VRMS/z_ABS3/IM_VRMSxz_'+sindex[sn]+'.jpg'






