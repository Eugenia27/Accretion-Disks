import matplotlib as mpl 
import matplotlib.pyplot as plt
from matplotlib.mlab import griddata
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy       as np
from mpl_toolkits.axes_grid.inset_locator import (inset_axes, InsetPosition, mark_inset)




def plots_IMSHOW(_snap,_x1,_x2,_x3,_property,_limx1,_limx2,_limx3,_limproperty,_ngrid,_ncont_levels,_prop_type,_paths,_filename):
    print 'estoy en 1'
    limx1= _limx1
    limx2= _limx2
    limx3= _limx3
    limp =_limproperty
    
    x1p = _x1[      (abs(_x1)>limx1[0]) & (abs(_x1)<limx1[1])]
    x2p = _x2[      (abs(_x1)>limx1[0]) & (abs(_x1)<limx1[1])]
    x3p = _x3[      (abs(_x1)>limx1[0]) & (abs(_x1)<limx1[1])]
    pp  = _property[(abs(_x1)>limx1[0]) & (abs(_x1)<limx1[1])]

    x1  = x1p[(x2p<limx2[1]) & (x2p>limx2[0]) & (x3p<limx3[1]) & (x3p>limx3[0])]
    x2  = x2p[(x2p<limx2[1]) & (x2p>limx2[0]) & (x3p<limx3[1]) & (x3p>limx3[0])]
    x3  = x3p[(x2p<limx2[1]) & (x2p>limx2[0]) & (x3p<limx3[1]) & (x3p>limx3[0])]
    p   =  pp[(x2p<limx2[1]) & (x2p>limx2[0]) & (x3p<limx3[1]) & (x3p>limx3[0])]

    ng = _ngrid
    ncl= _ncont_levels
    
    #plot
    fig = plt.figure(figsize=(20,10))
    ax = fig.add_subplot(111)

    [X, Y]    = np.meshgrid(np.linspace(min(x1),max(x1),ng), np.linspace(min(x2),max(x2),ng))
    Z         = griddata(x1,x2,p,X,Y,interp='linear')
    map_prop  = ax.imshow(Z,vmin=limp[0],vmax=limp[1])
    cont_prop = ax.contour(Z,levels=np.linspace(limp[0],limp[1],ncl), cmap='gist_rainbow_r' )

    locsx,labelsx = plt.xticks()
    labelsx = [str(i*(2*limx1[1]/float(ng))-limx1[1]) for i in locsx]
    tick_locsx = locsx
    tick_lblsx = labelsx
    plt.xticks(tick_locsx[1:], tick_lblsx[1:])

    locsy,labelsy = plt.yticks()
    labelsy = [str(i*(2*limx2[1]/float(ng))-limx2[1]) for i in locsy]
    tick_locsy = locsy
    tick_lblsy = labelsy
    plt.yticks(tick_locsy[1:], tick_lblsy[1:])

    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)

    pt=_prop_type
    
    if pt=="rho":
        textbar = r"$log$ $\rho/(M_{\odot}/kpc^{3})$"
    elif pt=="vrms":    
        textbar = "$log$ $VRMS/(km/s)$"
    elif pt=="u":
        textbar = "$log$  $U/(km/s)^{2}$"

    cbar = fig.colorbar(map_prop,cax=cax, label=textbar,ticks=np.linspace(limp[0],limp[1],ncl))
    cbar.add_lines(cont_prop)
    
    time=_snap[9::]
    ax.minorticks_on()
    plt.text(0.87, 0.92,"t="+time, ha='center', va='center', transform=ax.transAxes,color="#ffa74e")
    
    path=_paths+pt+str(time)+'.jpg'
    plt.savefig(path, dpi=100, bbox_inches='tight')
    plt.show()


def plots_IMSHOW(_snap,_x1,_x2,_property,_limx1,_limx2,_limproperty,_cond1,_cond2,_limcond1,_limcond2,_ngrid,_ncont_levels,_prop_type,_paths,_filename):
    print 'estoy en 2'
    limx1= _limx1
    limx2= _limx2
    limp =_limproperty
    limc1=_limcond1
    limc2=_limcond2

    x1p = _x1
    x2p = _x2
    pp  = _property
    
    x1  = x1p[(_cond1<limc1[1]) & (_cond1>limc1[0]) & (_cond2<limc2[1]) & (_cond2>limc2[0])]
    x2  = x2p[(_cond1<limc1[1]) & (_cond1>limc1[0]) & (_cond2<limc2[1]) & (_cond2>limc2[0])]
    p   =  pp[(_cond1<limc1[1]) & (_cond1>limc1[0]) & (_cond2<limc2[1]) & (_cond2>limc2[0])]

    ng = _ngrid
    ncl= _ncont_levels
    
    #plot
    fig = plt.figure(figsize=(20,10))
    ax = fig.add_subplot(111)

    [X, Y]    = np.meshgrid(np.linspace(min(x1),max(x1),ng), np.linspace(min(x2),max(x2),ng))
    Z         = griddata(x1,x2,p,X,Y,interp='linear')
    map_prop  = ax.imshow(Z,vmin=limp[0],vmax=limp[1])
    cont_prop = ax.contour(Z,levels=np.linspace(limp[0],limp[1],ncl), cmap='gist_rainbow_r' )

    locsx,labelsx = plt.xticks()
    labelsx = [str(i*(2*limx1[1]/float(ng))-limx1[1]) for i in locsx]
    tick_locsx = locsx
    tick_lblsx = labelsx
    plt.xticks(tick_locsx[1:], tick_lblsx[1:])

    locsy,labelsy = plt.yticks()
    labelsy = [str(i*(2*limx2[1]/float(ng))-limx2[1]) for i in locsy]
    tick_locsy = locsy
    tick_lblsy = labelsy
    plt.yticks(tick_locsy[1:], tick_lblsy[1:])

    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)

    pt=_prop_type
    
    if pt=="rho":
        textbar = r"$log$ $\rho/(M_{\odot}/kpc^{3})$"
    elif pt=="vrms":    
        textbar = "$log$ $VRMS/(km/s)$"
    elif pt=="u":
        textbar = "$log$  $U/(km/s)^{2}$"

    cbar = fig.colorbar(map_prop,cax=cax, label=textbar,ticks=np.linspace(limp[0],limp[1],ncl))
    cbar.add_lines(cont_prop)
    
    time=_snap[9::]
    ax.minorticks_on()
    plt.text(0.87, 0.92,"t="+time, ha='center', va='center', transform=ax.transAxes,color="#ffa74e")
    
    path=_paths+pt+str(time)+'.jpg'
    plt.savefig(path, dpi=100, bbox_inches='tight')
    plt.show()



def plots_IMSHOW(_snap,_x1,_x2,_x3,_property,_c1,_c2,_conditions,_limx1,_limx2,_limx3,_limproperty,_ngrid,_ncont_levels,_prop_type,_paths,_filename):

    print 'estoy en 3'
    limx1= _limx1
    limx2= _limx2
    limx3= _limx3
    limp =_limproperty
    cond = _conditions
    
    x1p = _x1[      (abs(_x1)>limx1[0]) & (abs(_x1)<limx1[1])]
    x2p = _x2[      (abs(_x1)>limx1[0]) & (abs(_x1)<limx1[1])]
    x3p = _x3[      (abs(_x1)>limx1[0]) & (abs(_x1)<limx1[1])]
    pp  = _property[(abs(_x1)>limx1[0]) & (abs(_x1)<limx1[1])]
    c1p = _c1[      (abs(_x1)>limx1[0]) & (abs(_x1)<limx1[1])]
    c2p = _c2[      (abs(_x1)>limx1[0]) & (abs(_x1)<limx1[1])]

    x1  = x1p[(x2p<limx2[1]) & (x2p>limx2[0]) & (x3p<limx3[1]) & (x3p>limx3[0]) & (c1p>cond[0]) & (c2p>cond[1])]
    x2  = x2p[(x2p<limx2[1]) & (x2p>limx2[0]) & (x3p<limx3[1]) & (x3p>limx3[0]) & (c1p>cond[0]) & (c2p>cond[1])]
    x3  = x3p[(x2p<limx2[1]) & (x2p>limx2[0]) & (x3p<limx3[1]) & (x3p>limx3[0]) & (c1p>cond[0]) & (c2p>cond[1])]
    p   =  pp[(x2p<limx2[1]) & (x2p>limx2[0]) & (x3p<limx3[1]) & (x3p>limx3[0]) & (c1p>cond[0]) & (c2p>cond[1])]

    ng = _ngrid
    ncl= _ncont_levels
    
    #plot
    fig = plt.figure(figsize=(20,10))
    ax = fig.add_subplot(111)

    [X, Y]    = np.meshgrid(np.linspace(min(x1),max(x1),ng), np.linspace(min(x2),max(x2),ng))
    Z         = griddata(x1,x2,p,X,Y,interp='linear')
    map_prop  = ax.imshow(Z,vmin=limp[0],vmax=limp[1])
    cont_prop = ax.contour(Z,levels=np.linspace(limp[0],limp[1],ncl), cmap='gist_rainbow_r' )

    locsx,labelsx = plt.xticks()
    labelsx = [str(i*(2*limx1[1]/float(ng))-limx1[1]) for i in locsx]
    tick_locsx = locsx
    tick_lblsx = labelsx
    plt.xticks(tick_locsx[1:], tick_lblsx[1:])

    locsy,labelsy = plt.yticks()
    labelsy = [str(i*(2*limx2[1]/float(ng))-limx2[1]) for i in locsy]
    tick_locsy = locsy
    tick_lblsy = labelsy
    plt.yticks(tick_locsy[1:], tick_lblsy[1:])

    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)

    pt=_prop_type
    
    if pt=="rho":
        textbar = r"$log$ $\rho/(M_{\odot}/kpc^{3})$"
    elif pt=="vrms":    
        textbar = "$log$ $VRMS/(km/s)$"
    elif pt=="u":
        textbar = "$log$  $U/(km/s)^{2}$"

    cbar = fig.colorbar(map_prop,cax=cax, label=textbar,ticks=np.linspace(limp[0],limp[1],ncl))
    cbar.add_lines(cont_prop)
    
    time=_snap[9::]
    ax.minorticks_on()
    plt.text(0.87, 0.92,"t="+time, ha='center', va='center', transform=ax.transAxes,color="#ffa74e")
    
    path=_paths+pt+str(time)+'.jpg'
    plt.savefig(path, dpi=100, bbox_inches='tight')
    plt.show()


def plots_IMSHOW_cc(_snap,_x1,_x2,_pr,_limx1,_limx2,_limpr,_ngrid,_ncont_levels,_text,_paths,_filename):

    limx1   = _limx1
    limx2   = _limx2
    limp    = _limpr
    textbar = _text

    x1  = _x1
    x2  = _x2
    p   = _pr

    ng = _ngrid
    ncl= _ncont_levels

    #plot
    fig = plt.figure(figsize=(10,10))
    ax = fig.add_subplot(111)

    [X, Y]    = np.meshgrid(np.linspace(np.min(x1),np.max(x1),ng), np.linspace(np.min(x2),np.max(x2),ng))
    Z         = griddata(x1,x2,p,X,Y,interp='linear')
    map_prop  = ax.imshow(Z,vmin=limp[0],vmax=limp[1],cmap='jet',origin='lower')
    cont_prop = ax.contour(Z,levels=np.linspace(limp[0],limp[1],ncl),cmap='Greys_r',linewidths=1.2)#cmap='gist_rainbow_r',lw=0.2 )
    plt.clabel(cont_prop, inline=1, fontsize=20, colors='k')
    
    locsx,labelsx = plt.xticks()
    labelsx = [str(i*(2*limx1[1]/float(ng))-limx1[1]) for i in locsx]
    tick_locsx = locsx
    tick_lblsx = labelsx
    plt.xticks(tick_locsx[1:], tick_lblsx[1:])

    locsy,labelsy = plt.yticks()
    labelsy = [str(i*(2*limx2[1]/float(ng))-limx2[1]) for i in locsy]
    tick_locsy = locsy
    tick_lblsy = labelsy
    plt.yticks(tick_locsy[1:], tick_lblsy[1:])

    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)       

    cbar = fig.colorbar(map_prop,cax=cax, label=textbar,ticks=np.linspace(limp[0],limp[1],ncl))
    cbar.add_lines(cont_prop)    
    time=_snap[9::]
    ax.minorticks_on()
    plt.text(0.87, 0.92,"t="+time, ha='center', va='center', transform=ax.transAxes,color="#ffa74e")    
    path=_paths+'cc_'+'_'+_filename+'_'+str(time)+'.jpg'
    plt.savefig(path, dpi=100, bbox_inches='tight')
    #plt.show()



def visible_ticks(_xticks,_yticks):
    xticks = _xticks
    for i in range(len(xticks)):
        if i%2==0:
            xticks[i].label1.set_visible(False)            
    yticks = _yticks
    for i in range(len(yticks)):
        if i%2!=0:
            yticks[i].label1.set_visible(False)



def plot_profile(ax,binn,val,disp,limx,limy,alabels,llabels,snap,s1,s2,save):
    ax.fill_between(binn,val-disp,val+disp,edgecolor='grey',hatch='//////',facecolor='none',alpha=0.3,label=llabels[0])
    ax.plot(binn,val,c='navy',lw=1,marker='o',ms=3,mfc='white',label=llabels[1],mec='navy')
    ax.set_xlabel(alabels[0])
    ax.set_ylabel(alabels[1])
    ax.set_xlim(limx[0],limx[1])
    ax.set_ylim(limy[0],limy[1])
    xticks,yticks = ax.xaxis.get_major_ticks(),ax.yaxis.get_major_ticks()       
    xticks[0].label1.set_visible(False) 
    yticks[0].label1.set_visible(False) 
    xticks[-1].label1.set_visible(False) 
    yticks[-1].label1.set_visible(False) 
    ax.minorticks_on()    
    time=snap[9::]
    plt.text(0.87, 0.05,"t="+time, ha='center', va='center', transform=ax.transAxes,color="black")
    ax.legend(loc='upper right', frameon=False,ncol=1, scatterpoints=1,borderpad=0.1, labelspacing=0.2, handletextpad=0.1)
    ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    if save==1:
        plt.savefig(s1+'.pdf', dpi=100, bbox_inches='tight')
        plt.savefig(s2+'.png', dpi=100, bbox_inches='tight')
    
    plt.show()


def plot_profile_with_subplot(ax,binn,val,disp,limx,limy,limxsp,limysp,alabels,llabels,snap,s1,s2,save):
    ax.fill_between(binn,val-disp,val+disp,edgecolor='grey',hatch='//////',facecolor='none',alpha=0.3,label=llabels[0])
    ax.plot(binn,val,c='navy',lw=1,marker='o',ms=3,mfc='white',label=llabels[1],mec='navy')
    ax.set_xlabel(alabels[0])
    ax.set_ylabel(alabels[1])
    ax.set_xlim(limx[0],limx[1])
    ax.set_ylim(limy[0],limy[1])
    xticks,yticks = ax.xaxis.get_major_ticks(),ax.yaxis.get_major_ticks()       
    xticks[0].label1.set_visible(False) 
    yticks[0].label1.set_visible(False) 
    xticks[-1].label1.set_visible(False) 
    yticks[-1].label1.set_visible(False) 
    ax.minorticks_on()    
    time=snap[9::]
    ax.text(0.87, 0.05,"t="+time, ha='center', va='center', transform=ax.transAxes,color="black")
    ax.legend(loc='upper left', frameon=False,ncol=1, scatterpoints=1,borderpad=0.1, labelspacing=0.2, handletextpad=0.1)
    ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    
    
    ax_sp = plt.axes([0,0,1,1])
    ip = InsetPosition(ax, [0.45,0.65,0.5,0.30])
    ax_sp.set_axes_locator(ip)
    ax_sp.plot(binn,val,c='navy',lw=1,marker='o',ms=3,mfc='white',mec='navy')
    ax_sp.set_xlim(limxsp[0], limxsp[1]) 
    ax_sp.set_ylim(limysp[0], limysp[1])
    xticks_sp,yticks_sp = ax_sp.xaxis.get_major_ticks(),ax_sp.yaxis.get_major_ticks() 
    visible_ticks(xticks_sp,yticks_sp) 
    ax_sp.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    ax_sp.tick_params(labelsize=10)
    ax_sp.minorticks_on()    
    
    
    if save==1:
        plt.savefig(s1+'.pdf', dpi=100, bbox_inches='tight')
        plt.savefig(s2+'.png', dpi=100, bbox_inches='tight')
    
    plt.show()








