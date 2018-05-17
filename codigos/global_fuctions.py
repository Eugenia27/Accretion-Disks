import numpy as np


def sobgroup(_var,_limiter_column,_limits):
    var_ = []
    lc   = _limiter_column
    lim  = _limits
    for j in range(len(_var)):
        var_.append([])
    for i in range(len(_var)):
        var_[i].append(_var[i][(_var[lc]<lim[1])&(_var[lc]>lim[0])])

    return var_

def subgroup_all(_sgvar):
    sgx_     =_sgvar[0][0]
    sgy_     =_sgvar[1][0]
    sgz_     =_sgvar[2][0]
    sgr_     =_sgvar[3][0]
    sgrho_   =_sgvar[4][0]
    sgvx_    =_sgvar[5][0]
    sgvy_    =_sgvar[6][0]
    sgvz_    =_sgvar[7][0]
    sgbx_    =_sgvar[8][0]
    sgby_    =_sgvar[9][0]
    sgbz_    =_sgvar[10][0]
    sgvrms_  =_sgvar[11][0]
    sgeint_  =_sgvar[12][0] 
    return    sgx_,sgy_,sgz_,sgr_,sgrho_,sgvx_,sgvy_,sgvz_,sgbx_,sgby_,sgbz_,sgvrms_,sgeint_ 


def subgroup_all_wm(_sgvar):
    sgx_     =_sgvar[0][0]
    sgy_     =_sgvar[1][0]
    sgz_     =_sgvar[2][0]
    sgr_     =_sgvar[3][0]
    sgrho_   =_sgvar[4][0]
    sgmass_  =_sgvar[5][0]
    sgvx_    =_sgvar[6][0]
    sgvy_    =_sgvar[7][0]
    sgvz_    =_sgvar[8][0]
    sgbx_    =_sgvar[9][0]
    sgby_    =_sgvar[10][0]
    sgbz_    =_sgvar[11][0]
    sgvrms_  =_sgvar[12][0]
    sgeint_  =_sgvar[13][0] 
    return    sgx_,sgy_,sgz_,sgr_,sgrho_,sgmass_,sgvx_,sgvy_,sgvz_,sgbx_,sgby_,sgbz_,sgvrms_,sgeint_ 

def subgroup_all_extended(_sgvar):
    sgx_     =_sgvar[0][0]
    sgy_     =_sgvar[1][0]
    sgz_     =_sgvar[2][0]
    sgr_     =_sgvar[3][0]
    sgrho_   =_sgvar[4][0]
    sgmass_  =_sgvar[5][0]
    sgvx_    =_sgvar[6][0]
    sgvy_    =_sgvar[7][0]
    sgvz_    =_sgvar[8][0]
    sgbx_    =_sgvar[9][0]
    sgby_    =_sgvar[10][0]
    sgbz_    =_sgvar[11][0]
    sgvrms_  =_sgvar[12][0]
    sgeint_  =_sgvar[13][0]
    sgides_  =_sgvar[14][0] 
    return    sgx_,sgy_,sgz_,sgr_,sgrho_,sgmass_,sgvx_,sgvy_,sgvz_,sgbx_,sgby_,sgbz_,sgvrms_,sgeint_,sgides_ 

def subgroup_all_crun(_sgvar):
    sgx_     =_sgvar[0][0]
    sgy_     =_sgvar[1][0]
    sgz_     =_sgvar[2][0]
    sgr_     =_sgvar[3][0]
    sgrho_   =_sgvar[4][0]
    sgmass_  =_sgvar[5][0]
    sgvx_    =_sgvar[6][0]
    sgvy_    =_sgvar[7][0]
    sgvz_    =_sgvar[8][0]
    sgvrms_  =_sgvar[9][0]
    sgeint_  =_sgvar[10][0]
    return    sgx_,sgy_,sgz_,sgr_,sgrho_,sgmass_,sgvx_,sgvy_,sgvz_,sgvrms_,sgeint_ 

