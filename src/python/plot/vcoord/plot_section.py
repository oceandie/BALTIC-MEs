#!/usr/bin/env python

#            |--------------------------------------------------------------------------|
#            | This module contains all the functions needed to plot data for or from   |
#            | the NEMO model. In particular, it contains the following functions:      |
#            |                                                                          |
#            |          *) mpl_sec_get_section                                          |
#            |          *) mpl_sec_int_s2z                                              |
#            |          *) mpl_sec_bathy                                                |
#            |          *) mpl_sec_settings                                             |
#            |          *) mpl_sec_draw_scalar                                          |
#            |          *) mpl_sec_draw_levels                                          |
#            |          *) mpl_sec_draw_grid_pnts                                       |
#            |          *) mpl_sec_figure                                               |
#            |          *) mpl_lev                                                      |
#            |          *) mpl_sec                                                      |
#            |          *) mpl_plot_sec_loop                                            |
#            |                                                                          |
#            | Author: D. Bruciaferri                                                   |
#            | Date and place: 13th Feb 2016, Plymouth University, United Kindom        |
#            |--------------------------------------------------------------------------|

import numpy as np
import utils as utl

from scipy.interpolate import interp1d

import matplotlib as matpl
import matplotlib.pyplot as plt
import matplotlib.tri as tri
import matplotlib.path as path
import matplotlib.patches as patches
import matplotlib.ticker as mtick
import matplotlib.colors as colors

import cartopy.crs as ccrs
import cartopy.feature as feature

transform = ccrs.PlateCarree()
proj = ccrs.Mercator()

#===================================================================================================

def mpl_sec_get_section(lon3, lat3, rbat2, zenv2, sec_i, sec_j, coord_type, \
                        var3=[], tmsk3=[], tdep3=[], wdep3=[], tpnt3=[], wlev3=[], var3_aux=[], msk_mes=None):

    '''
    This function returns x-axis coordinates data for a given section
    specified by sec_i and sec_j.

    If wlev3, rbat2 and/or zenv2 are given, the same section for these fields 
    is returned as well. 
    '''
   

    nk = lon3.shape[0]
    nj = lon3.shape[1]
    ni = lon3.shape[2]

    if isinstance(sec_i, list) and isinstance(sec_j, list):

       if not len(sec_i) == len(sec_j):
          print("")
          print("mpl_sec_get_section() ERROR:")
          print("section indexes lists do not")
          print("have the same length.")
          print("Please check")
          print("")
          return

       secJ, secI  = utl.get_poly_line_ij(sec_i, sec_j)
       sec_name    =  "lon_" + str(sec_i[0]) + "-" + str(sec_i[-1]) + \
                     "_lat_" + str(sec_j[0]) + "-" + str(sec_j[-1])
    else:
    
       if sec_i == -1:

          secI        = np.arange(ni)
          secJ        = np.multiply(np.ones(ni, dtype=np.int),sec_j)
          sec_name    = "j_" + str(sec_j)

       if sec_j == -1:

          secI        = np.multiply(np.ones(nj, dtype=np.int),sec_i)
          secJ        = np.arange(nj)
          sec_name    = "i_" + str(sec_i)

    sectLen = len(secI)
    
    xcoord  = np.zeros(shape=(nk,sectLen))

    if rbat2 != []: 
       ycoordB = np.zeros(shape=(sectLen))
    else:
       ycoordB = []

    if zenv2 != []:
       num_env = len(zenv2) 
       ycoordH = np.zeros(shape=(num_env,sectLen))
    else:
       ycoordH = []

    if var3  != []: 
       var2 = np.zeros(shape=(nk,sectLen))
    else:
       var2 = []

    if tmsk3 != []: 
       tmsk2 = np.zeros(shape=(nk,sectLen))
    else:
       tmsk2 = []

    if tdep3 != []: 
       tdep2 = np.zeros(shape=(nk,sectLen))
    else:
       tdep2 = []

    if wdep3 != []: 
       wdep2 = np.zeros(shape=(nk,sectLen))
    else:
       wdep2 = []

    if tpnt3 != []:
       tpnt2 = np.zeros(shape=(nk,sectLen))
    else:
       tpnt2 = []

    if wlev3 != []:
       wlev2 = np.zeros(shape=(nk,sectLen))
    else:
       wlev2 = []

    if var3_aux  != []:
       var2_aux = np.zeros(shape=(nk,sectLen))
    else:
       var2_aux = []

    if msk_mes is not None:
       msk_mes2 = np.zeros(shape=(nk,sectLen))
    else:
       msk_mes2 = []

    lon2 = lon3[0,:,:]
    lat2 = lat3[0,:,:]

    tot_dist = 0.0
    last_lon = lon2[secJ[0],secI[0]]
    last_lat = lat2[secJ[0],secI[0]]

    lon_sec = []
    lat_sec = []
    lon_sec.append(last_lon)
    lat_sec.append(last_lat)
    # CYCLING OVER THE COORDINATES ------------------------------------------
    for s in np.arange(sectLen):
        i=secI[s]
        j=secJ[s]

        if i < ni:
           i1 = i + 1
        else:
           i1 = i - 1
        if j < nj:
           j1 = j + 1
        else:
           j1 = j - 1
        if (lon2[j ,i ] == 0 and lat2[j ,i ] == 0) and \
           (lon2[j1,i1] == 0 and lat2[j1,i1] == 0):
           break

        this_lon = lon2[j,i]
        this_lat = lat2[j,i]
        lon_sec.append(this_lon)
        lat_sec.append(this_lat)
        if var2     != []: var2[:,s]     = var3[:,j,i]
        if tmsk2    != []: tmsk2[:,s]    = tmsk3[:,j,i]
        if tdep2    != []: tdep2[:,s]    = tdep3[:,j,i]
        if wdep2    != []: wdep2[:,s]    = wdep3[:,j,i]
        if tpnt2    != []: tpnt2[:,s]    = tpnt3[:,j,i]
        if wlev2    != []: wlev2[:,s]    = wlev3[:,j,i]
        if var2_aux != []: var2_aux[:,s] = var3_aux[:,j,i]
        if msk_mes2 != []: msk_mes2[:,s] = msk_mes[j,i]

        if coord_type == "dist":
           dist_from_last = utl.hvrsn_dst(last_lon, last_lat, this_lon, this_lat) / 1000.; # in km
           tot_dist = tot_dist + dist_from_last;
           xcoord[:,s] = tot_dist
           last_lon = this_lon
           last_lat = this_lat

        elif coord_type == "index":
           if isinstance(sec_i, list) and isinstance(sec_j, list):
              print("mpl_sec_get_section() ERROR:")
              print("You have a user defined section.")
              print('Please chose coord_type = "dist"')
              return
           if sec_i == -1: xcoord[:,s] = i
           if sec_j == -1: xcoord[:,s] = j

        if rbat2 != []: ycoordB[s] = rbat2[j,i]
        if zenv2 != []: 
           for nenv in range(num_env):
               ycoordH[nenv,s] = zenv2[nenv][j,i]
    #------------------------------------------------------------------------
    if coord_type == "dist": xlabel = "Distance [Km]"

    if coord_type == "index":
       if sec_i == -1: xlabel = "i index"    
       if sec_j == -1: xlabel = "j index"       

    return xcoord, ycoordB, ycoordH, xlabel, sec_name, var2, tmsk2, tdep2, wdep2, tpnt2, wlev2, var2_aux, lon_sec, lat_sec, msk_mes2   

#===================================================================================================

def mpl_sec_int_s2z(cn_RVlevel, var2, tdep2, x2):

    '''
    This function interpolates values defined on
    sigma levels into z levels.
    '''

    if cn_RVlevel == 'auto':

       reg_dep1 = np.mean(tdep2, axis=1)

    else:

       reg_dep1 = cn_RVlevel

    reg_dep2 = np.repeat(reg_dep1[:, np.newaxis], x2.shape[1],       axis=1)
    reg_x2   = np.repeat(x2[0,:][np.newaxis,:],   reg_dep2.shape[0], axis=0)

    reg_var2 = np.zeros(shape=(reg_x2.shape[0],reg_x2.shape[1]))

    #----------------------------------------------------------------------------------------
    # The interpolation function interpolate.interp1d has been compared with the function
    # com.vert_1Dslice and the provide the exact same reults, but interpolate.interp1d is
    # much faster.
    # ---------------------------------------------------------------------------------------
    # com.vert_1Dslice:
    #
    #          for k in np.arange(1,reg_x2.shape[0]):
    #              reg_var2[k,:] = com.vert_1Dslice(Tvar2new, Tdep2new, reg_dep2[k,0])
    #          reg_var2[0,:] = Tvar2new[0,:]
    #
    # interpolate.interp1d:
    #
    for i in np.arange(reg_x2.shape[1]):
        f_var          = interp1d(Tdep2new[:,i], Tvar2new[:,i], kind='linear')
        reg_var2[:,i]  = f_var(reg_dep2[:,i])
    #
    # ---------------------------------------------------------------------------------------

    return reg_x2, reg_dep2, reg_var2 

#===================================================================================================

def mpl_sec_bathy(rbat, rbat_fill, zenv, x2, tdep2, wdep2, tmsk2, mbat_fill, \
                  vlevel="", PlotType="", rbat_max=0, msk_loc=[]):

    '''
    This function returns patches of
    real bathymetry (rbat) and
    model bathymetry (mbat) to be used
    in sections plots.
    '''

    rpatch = []
    mpatch = []
    RB     = []
    MB     = []
    ZE     = []

    # 1. REAL BATHYMETRY VERTEXES

    if rbat != []:
       rbat_x = []
       rbat_z = []
       for i in np.arange(len(rbat)):
           rbat_x.append(x2[0,i])    
           if np.isnan(rbat[i]): 
              rbat_z.append(0.)
           else:
              rbat_z.append(rbat[i])
       rbat_x.append(x2[-1,-1])
       rbat_x.append(x2[-1,0])
       if tdep2 == [] or wdep2 == [] or tmsk2 == []:
          rbat_z.append(rbat_max)
          rbat_z.append(rbat_max) 
       if tdep2 != [] and wdep2 != [] and tmsk2 != []:
          if vlevel == "S_re": # or vlevel == "MES":
             rbat_z.append(rbat_max)
             rbat_z.append(rbat_max)
          else:
             rbat_z.append(tdep2[-1,-1] + 0.5*(tdep2[-1,-1]-tdep2[-2,-1]))
             rbat_z.append(tdep2[-1, 0] + 0.5*(tdep2[-1, 0]-tdep2[-2, 0]))
    # 2. REAL BATHYMETRY LINE AND PATCH
    #RB = plt.plot(rbat_x[:-2], rbat_z[:-2], color='black', label='Bathymetry')
    ##RB = plt.plot(rbat_x[:], rbat_z[:], color='black', label='Bathymetry')
    if rbat_fill == "true":
       vertexes = list(zip(rbat_x, rbat_z))
       nverts = len(vertexes)
       codes     = np.ones(nverts, int) * path.Path.LINETO
       codes[0]  = path.Path.MOVETO
       polyg     = path.Path(vertexes, codes)
       rpatch    = patches.PathPatch(polyg, facecolor='black', edgecolor='black', alpha=1.,zorder=10)

    # 3. MODEL BATHYMETRY VERTEXES

    if tdep2 != [] and wdep2 != [] and tmsk2 != [] and vlevel != "no":

       msk_oce = np.zeros(shape=tmsk2.shape)
       msk_oce[tmsk2 == 0] = 1
       msk_oce[tmsk2 == 1] = np.nan
       if msk_loc !=[]:
          msk_oce[msk_loc == 0] = np.nan
          msk_oce[-1,:] = 1

       # PCOLOR OR LEVELS PLOT CASES
       #if PlotType != "contourf" or vlevel == "Z_ps":

       if vlevel == "S_re": #or vlevel == "MES":
          #mbat_x, mbat_z = utl.create_model_bathy_sec(vlevel, msk_oce, x2, tdep2, wdep2, rbat_max)
          mbat_x = rbat_x
          mbat_z = rbat_z
       elif vlevel == "MES" or vlevel == "SZT":
          mbat_x, mbat_z = utl.create_model_bathy_sec(vlevel, msk_oce, x2, tdep2, wdep2, rbat_max)
       else:
          mbat_x, mbat_z = utl.create_model_bathy_sec(vlevel, msk_oce, x2, tdep2, wdep2)
       # 4. MODEL BATHYMETRY LINE AND PATCH 
       MB = plt.plot(mbat_x[:-2], mbat_z[:-2], color='grey', label='Model Bathymetry')
       plt.setp(MB, 'linewidth', 0.2)
       if mbat_fill == "true":
             vertexes = list(zip(mbat_x, mbat_z))
             nverts = len(vertexes)
             codes     = np.ones(nverts, int) * path.Path.LINETO
             codes[0]  = path.Path.MOVETO
             polyg     = path.Path(vertexes, codes)
             mpatch    = patches.PathPatch(polyg, facecolor='grey', edgecolor='grey', alpha=1.,zorder=1)

    # 5. Z-ENV VERTEXES

    if zenv != []:
       num_env = zenv.shape[0]
       for nenv in range(num_env):
           zenv_x = []
           zenv_z = []
           for i in np.arange(len(zenv[nenv,:])):
               zenv_x.append(x2[0,i])
               zenv_z.append(zenv[nenv,i])
           # 6. Z-ENV LINE
           ZE = plt.plot(zenv_x, zenv_z, color='magenta', label='Z-Envelope')
           plt.setp(ZE, 'linewidth', 3., zorder=20)    

    return rpatch, mpatch, RB, MB, ZE

#===================================================================================================

def mpl_sec_settings(var2, x2, dep2, varlim, xlim, ylim):

    if var2 != [] and varlim != []:

       if isinstance(varlim, (str, bytes)):
          vlim_max = np.nanmax(var2)
          vlim_min = np.nanmin(var2)
       else:
          vlim_min = varlim[0]
          vlim_max = varlim[1]
    else:

       vlim_max = []
       vlim_min = []
#-------------------------------------------------
    if isinstance(xlim, (str, bytes)):
       xlim_min = np.nanmin(x2)
       xlim_max = np.nanmax(x2)
    else:
       xlim_min = xlim[0]
       xlim_max = xlim[1]
#--------------------------------------------------

    if isinstance(ylim, (str, bytes)):

       ylim_min = np.nanmin(dep2)
       ylim_max = np.nanmax(dep2)

    else:

       ylim_min = ylim[0]
       ylim_max = ylim[1]


    return vlim_min, vlim_max, xlim_min, xlim_max, ylim_min, ylim_max

#===================================================================================================

def mpl_sec_draw_scalar(PlotType, var2, x2, dep2, vlim_min, vlim_max, \
                        colmap, cn_line, cn_level, cn_label, cn_color, cn_fill, var2_aux=[]):

    '''
    This function draws section scalar plots.
    '''
 
    var2_masked = np.ma.masked_where(np.isnan(var2), var2)
    if var2_aux != []: 
       var2_aux_masked = np.ma.masked_where(np.isnan(var2_aux), var2_aux)

    if isinstance(colmap, (str, bytes)):
       cmap = plt.cm.get_cmap(colmap)
    else:
       cmap = colmap

    #cmap.set_under('white')
    #cmap.set_bad("white")    

    # ------------------------------------------------------------------------------------------
    if PlotType == "pcolor":

       pc = plt.pcolormesh(x2, dep2, var2_masked, cmap=cmap, vmin=vlim_min, vmax=vlim_max)

    # ------------------------------------------------------------------------------------------
    if PlotType == "contourf":    

       if cn_line == 'true':

          if cn_color == 'true':
             pc = plt.contour(x2, dep2, var2_masked, cn_level, cmap = cmap,  \
                            vmin = vlim_min, vmax = vlim_max, extend = 'both')
          else:
             pc = plt.contour(x2, dep2, var2_masked, cn_level, colors = 'k', \
                            vmin = vlim_min, vmax = vlim_max, extend = 'both')

          if cn_label == 'true':
             plt.clabel(pc, inline=1, fontsize=25, fmt='%1.1f', colors='white')

       if cn_fill == 'true':

          norm = colors.BoundaryNorm(cn_level, cmap.N, clip=True)
          pc = plt.contourf(x2, dep2, var2_masked, cn_level, cmap = cmap, norm=norm,\
                            vmin = vlim_min, vmax = vlim_max, extend = 'both')

          #pc = plt.contourf(x2, dep2, var2_masked, cn_level, cmap = cmap, \
          #                  vmin = vlim_min, vmax = vlim_max, extend = 'neither')

       if var2_aux != []:
          rho_lev = [27.25,27.5,27.6,27.7,27.8,27.9,28.0,28.1,28.2]
          pc_aux = plt.contour(x2, dep2, var2_aux_masked, rho_lev, colors = 'k', 
                               linewidths=2., extend = 'both')
          pc_ovf = plt.contour(x2, dep2, var2_aux_masked, [27.8], colors = 'magenta', 
                               linewidths=5.0, extend = 'both')
          plt.clabel(pc_aux, inline=1, fontsize=25, fmt='%1.1f', colors='black')

    return cmap, pc

#===================================================================================================

def mpl_sec_draw_levels(x2, wlev2):

    nk = wlev2.shape[0]

    for k in np.arange(nk):

        W_l  = plt.plot(x2[k,:], wlev2[k,:], linestyle='--', color='k')
        plt.setp(W_l, 'linewidth', 0.5)

    return W_l

#===================================================================================================

def mpl_sec_draw_grid_pnts(x2, dep2, var2=[], vlim_min=[], vlim_max=[], cmap=[]):

    '''
    This function draws grid points, 
    eventually also scalar values 
    associated with them.
    '''

    if (var2 != [] and vlim_min != [] and vlim_max != [] and cmap != []):

       plt.scatter(np.ravel(x2), np.ravel(dep2), s=10, c=np.ravel(var2),\
                      marker='o',cmap=cmap, vmin=vlim_min, vmax=vlim_max)

    else:

       plt.scatter(np.ravel(x2), np.ravel(dep2), s=10, c="black", marker='+')

#===================================================================================================

def mpl_sec_figure(fig, pc, var_strng, title, xlabel, ylabel, cn_level=None):

    '''
    This function sets property of
    section plots fig.
    '''

    ax = fig.gca() # Get the current Axes instance on the current
                   # figure matching the given keyword args, or create one.

    ax.set_title(title, y=1.02, fontsize='25')

    ax.invert_yaxis()

    ax.set_xlabel(xlabel, fontsize='30')
    ax.set_ylabel(ylabel, fontsize='30')

    ax.tick_params(axis='x', labelsize=20)
    ax.tick_params(axis='y', labelsize=20)

    ax.xaxis.labelpad = 20
    ax.yaxis.labelpad = 20

    if pc != []:

       if cn_level is not None:
          cb = plt.colorbar(pc, ticks=cn_level)
       else:
          cb = plt.colorbar(pc)
       cb.ax.tick_params(labelsize=25)
       cb.set_label(var_strng, labelpad=20, size=25)

    else:

       cb = []

    return ax, cb

#===================================================================================================

def mpl_lev(exp_name, fig_type, sec_i, sec_j, tlon3, tlat3, tdep3, wdep3, tmsk3, m, \
            coord_type, vlevel, check='true', rbat2=[], zenv2=[], rbat2_fill="false", \
            mbat_ln="false", mbat_fill="false", xlim='no', ylim='no', xgrid='false', 
            msk_mes=None, first_zlv=None):

            if msk_mes is not None: print('LOCALISING')

            funcID   = "sec_"

            i3D_flag = True

            if tdep3 == [] or wdep3 == [] or tmsk3 == []:
               i3D_flag = False


            if i3D_flag: 
               print("TRANSECT OF NEMO-OGCm VERTICAL GEOMETRY")
               var_strng = 'vertical_levels'
            else:
               print("TRANSECT OF NEMO-OGCm BATHYMETRY")
               var_strng = 'bathymetry' 

            #==============================================
            # CHECKING CONSISTENCY OF ARRAYS DIMENSIONS

            # t-coordinates
            if i3D_flag:
               if not (tdep3.shape == tlat3.shape and \
                       tdep3.shape == tlon3.shape and \
                       tdep3.shape == wdep3.shape and \
                       tdep3.ndim == 3):

                       print("mpl_lev() ERROR:")
                       print("")
                       print("Please check t-coordinates")
                       print("matrixes dimensions")
                       print("i3D_flag = True")
                       return
            else:
               if not (tlon3.shape == tlat3.shape and \
                       tlon3.ndim == 3):

                       print("mpl_lev() ERROR:")
                       print("")
                       print("Please check t-coordinates")
                       print("matrixes dimensions")
                       print("i3D_flag = False")
                       return

            nk = tlon3.shape[0]
            nj = tlon3.shape[1]
            ni = tlon3.shape[2]

            print("MODEL DIMENSIONS: ", nk, " x ", nj, " x ", ni)

            #==============================================
            # GETTING VERTICAL LEVELS

            tpnt3 = []
            wlev3 = []
            if i3D_flag:
               tpnt3, wlev3 = utl.nemo_vgeom(vlevel, tdep3, wdep3)

            #==============================================
            # SELECTING THE REQUESTED SECTION

            if i3D_flag:
               var3 = np.copy(tmsk3)
               x2, rbat2_sec, zenv2_sec, xlabel, sec_name, var2, tmsk2, tdep2, wdep2, tpnt2, wlev2, var2_aux, lon_sec, lat_sec, msk_mes2 = \
               mpl_sec_get_section(tlon3, tlat3, rbat2, zenv2, sec_i, sec_j, coord_type,\
                                   var3 , tmsk3, tdep3, wdep3, tpnt3, wlev3, msk_mes=msk_mes)
            else:
               x2, rbat2_sec, zenv2_sec, xlabel, sec_name, var2, tmsk2, tdep2, wdep2, tpnt2, wlev2, var2_aux, lon_sec, lat_sec, msk_mes2 = \
               mpl_sec_get_section(tlon3, tlat3, rbat2, zenv2, sec_i, sec_j, coord_type,\
                                   [], [], [], [], tpnt3, wlev3, msk_mes=msk_mes)

            #==============================================
            # OPENING FIGURE

            fig = plt.figure(figsize = (34.,17.5), dpi=100)
      
            #==============================================
            # SETTING BATHYMETRY PATCH

            if i3D_flag:
               # We plot the model bathymetry creating a patch where the grid cells have T points
               # in the middle. Z_ps is a special case and it is maneged in a specific way.
               Tvar2new, Tx2new, Tdep2new = utl.nemo_ver_Tgrid4pcolor(var2, x2, tdep2, wdep2)

               if vlevel == "Z_ps":
                  for k in range(Tdep2new.shape[0]):
                      Tdep2new[k,:] = Tdep2new[k,int(Tdep2new.shape[1]/2.)]
                  rpatch, mpatch, RB, MB, ZE = mpl_sec_bathy(rbat2_sec, rbat2_fill, zenv2_sec, \
                                                             x2       , tdep2     , wdep2    , \
                                                             tmsk2    , mbat_fill , vlevel   , "")
               elif vlevel == "MES" or vlevel == "SZT":

                  if vlevel == "SZT" and first_zlv is not None:
                     msk_mes2[first_zlv::,:] = 0

                  rbat2_max = np.nanmax(rbat2+2000.)
                  rpatch, mpatch, RB, MB, ZE = mpl_sec_bathy(rbat2_sec, rbat2_fill, zenv2_sec, \
                                                             x2       , tdep2     , Tdep2new , \
                                                             tmsk2    , mbat_fill , vlevel   , \
                                                             ""       , rbat2_max, msk_mes2)
                  if msk_mes is not None:
                     msk_zps2 = np.ones(msk_mes2.shape) - msk_mes2
                     Tdep2loc = np.copy(Tdep2new)
                     for k in range(Tdep2loc.shape[0]):
                         Tdep2loc[k,:] = Tdep2loc[k,int(Tdep2loc.shape[1]/2.)]
                     rbat2_max = np.nanmax(rbat2+1000.)
                     rpatch_loc, mpatch_loc, RB_loc, MB_loc, ZE_loc = mpl_sec_bathy(rbat2_sec, rbat2_fill, zenv2_sec, \
                                                                                    x2       , tdep2     , wdep2    , \
                                                                                    tmsk2    , mbat_fill , "Z_ps"   , \
                                                                                    ""       , rbat2_max, msk_zps2)
               else:
                  rpatch, mpatch, RB, MB, ZE = mpl_sec_bathy(rbat2_sec, rbat2_fill, zenv2_sec, \
                                                             x2       , tdep2     , Tdep2new , \
                                                             tmsk2    , mbat_fill , vlevel   , "")
            else:
               rbat2_max = np.nanmax(rbat2)
               rpatch, mpatch, RB, MB, ZE = mpl_sec_bathy(rbat2_sec, rbat2_fill, zenv2_sec, \
                                                          x2       , tdep2     , tdep2    , \
                                                          tmsk2    , mbat_fill , ""       , \
                                                          ""       , rbat2_max)

            #==============================================
            # SETTING PLOT LIMITS

            if i3D_flag:
               vlim_min, vlim_max, xlim_min, xlim_max, ylim_min, ylim_max = \
               mpl_sec_settings([], x2, tdep2, [], xlim, ylim)
            else:
               vlim_min, vlim_max, xlim_min, xlim_max, ylim_min, ylim_max = \
               mpl_sec_settings([], x2, rbat2, [], xlim, ylim)

            #==============================================
            # DRAWING W VERTICAL LEVELS IF REQUESTED

            if i3D_flag:
               if msk_mes is not None:
                  wlev2_loc = np.copy(wlev2)
                  wlev2_loc[msk_zps2 == 0] = np.nan
                  for k in range(nk):
                      wlev_loc = wlev2_loc[k,:]
                      wlev_loc = wlev_loc[np.isfinite(wlev_loc)]
                      if wlev_loc.shape[0] > 0:
                         wlev2_loc[k,:] = np.unique(wlev_loc)[-1]
                  wlev2[msk_zps2 == 1] = wlev2_loc[msk_zps2 == 1]
               W_l = mpl_sec_draw_levels(x2, wlev2)

            #==============================================
            # DRAWING T GRID POINTS IF REQUESTED

            if i3D_flag and check == "true": 
               mpl_sec_draw_grid_pnts(x2, tpnt2)

            #==============================================
            # DRAWING X-GRID IF REQUESTED
     
            if i3D_flag and xgrid == "true": 
                                  
               for i in range(Tx2new.shape[1]):
                   X = plt.plot(Tx2new[:,i], Tdep2new[:,i], linestyle='--', color='k')
                   plt.setp(X, 'linewidth', 0.5)
               #X = plt.plot(Tx2new[:,220], Tdep2new[:,220], linestyle='-', color='blue')
               #plt.setp(X, 'linewidth', 2.5)
               #X = plt.plot(Tx2new[:,355], Tdep2new[:,i], linestyle='-', color='r')
            #==============================================
            # APPLYING PATCH IF REQUESTED

            title  = exp_name + "\n" + "Transect along " + sec_name

            ylabel = "Depth [m]"
            ax, cb = mpl_sec_figure(fig, [], "", title, xlabel, ylabel)


            if mpatch != []: ax.add_patch(mpatch)
            if rpatch != []: ax.add_patch(rpatch)
            if msk_mes is not None: ax.add_patch(mpatch_loc)

            #==============================================
            # LIMITS OF MAIN SECTION

            plt.axis([xlim_min, xlim_max, ylim_max, ylim_min])
            pos1 = ax.get_position()
            
            #==============================================
            # INSET MAP

            # if proj != []:
            #    lon_sec = np.asarray(lon_sec)
            #    lat_sec = np.asarray(lat_sec)
            #    lon0 = np.nanmin(lon_sec)
            #    lon1 = np.nanmax(lon_sec)
            #    lat0 = np.nanmin(lat_sec)
            #    lat1 = np.nanmax(lat_sec)
            #    if np.absolute(lon0-lon1) < 1. or np.absolute(lat0-lat1) < 1.:
            #       stp = 5.
            #    else:
            #       stp = 2.
            #    lon0 += -stp
            #    lon1 +=  stp
            #    lat0 += -stp
            #    lat1 +=  stp
            #    map_lims = [lon0, lon1, lat0, lat1]
            #    a = plt.axes([pos1.x0+0.005, pos1.y0+0.01, .2, .2], projection=proj)
            #    a.coastlines()
            #    a.add_feature(feature.LAND, color='gray',edgecolor='gray',zorder=1)
            #    MAP = plt.plot(lon_sec, lat_sec, c="red", transform=transform)
            #    plt.setp(MAP, 'linewidth', 2.5)
            #    a.set_extent(map_lims)

            #==============================================
            # SAVING FIGURE
            name_fig = funcID + "_" + var_strng + "_" + sec_name + "_" +\
                       "_maxdepth_" + str(ylim_max) + fig_type

            fig.savefig(name_fig,bbox_inches="tight", pad_inches=0)

            plt.close()

            return

#===================================================================================================

def mpl_sec(exp_name, fig_type, var_strng, unit_strng, date, time, timestep, PlotType, sec_i, sec_j,
            tlon3, tlat3, tdep3, wdep3, tmsk3, var3, m, coord_type, vlevel, rbat2=[], zenv2=[],
            rbat2_fill="false", mbat_ln="false", mbat_fill="false", xlim='no', ylim='no', varlim='no',
            colmap='jet', check='false', check_val='false', xgrid='false', cn_level=20, cn_line='true',
            cn_label='true', cn_color='false', cn_fill='true', cn_RegVertLev='false', cn_RVlevel='auto', 
            var_aux=None, msk_mes=None):

    funcID   = "sec_"

    print("TRANSECT of NEMO-OGCm OUTPUT: ") #+ var_strng
    if msk_mes is not None: print('LOCALISING')

    #==============================================
    # CHECKING CONSISTENCY OF ARRAYS DIMENSIONS

    # t-coordinates
    if not (tdep3.shape == tlat3.shape and \
            tdep3.shape == tlon3.shape and \
            tdep3.shape == wdep3.shape and \
            tdep3.ndim == 3):
       print("mpl_sec() ERROR:")
       print("")
       print("Please check t-coordinates")
       print("matrixes dimensions")
       return

    if not (var3.shape == tdep3.shape and \
            tmsk3.shape == tdep3.shape):
       print("mpl_sec() ERROR:")
       print("")
       print("Please check scalar variable")
       print("matrix dimensions")
       return

    nk = tlon3.shape[0]
    nj = tlon3.shape[1]
    ni = tlon3.shape[2]

    print("3D VARIABLES SHAPE: ", nk, " x ", nj, " x ", ni)

    #==============================================
    # MANAGING VERTICAL LEVELS

    tpnt3 = []
    wlev3 = []
    if vlevel != "no":
       tpnt3, wlev3 = utl.nemo_vgeom(vlevel, tdep3, wdep3)
    elif msk_mes is not None:
       tpnt3, wlev3 = utl.nemo_vgeom('MES', tdep3, wdep3)

    #==============================================
    # SELECTING THE REQUESTED SECTION    

    if var_aux is not None: 
       var3_aux = var_aux
    else:
       var3_aux = []

    x2, rbat2_sec, zenv2_sec, xlabel, sec_name, var2, tmsk2, tdep2, wdep2, tpnt2, wlev2, var2_aux, lon_sec, lat_sec, msk_mes2 = \
    mpl_sec_get_section(tlon3, tlat3, rbat2, zenv2, sec_i, sec_j, coord_type,\
                        var3, tmsk3, tdep3, wdep3, tpnt3, wlev3, var3_aux, msk_mes)

    #==============================================
    # CREATING VERTICAL MATRIXES CENTERED ON T POINTS

    if PlotType == "pcolor":

       Tvar2new, Tx2new, Tdep2new = utl.nemo_ver_Tgrid4pcolor(var2, x2, tdep2, wdep2)
       Tx2   = np.copy(Tx2new)
       Tdep2 = np.copy(Tdep2new)
       Tvar2 = np.copy(Tvar2new)
       Tmsk2 = np.copy(tmsk2)

    elif PlotType == "contourf":

       # Extending matrixes for values at surface and bottom
       Tx2new   = np.zeros(shape=(x2.shape[0]+2,    x2.shape[1]))
       Tdep2new = np.zeros(shape=(tdep2.shape[0]+2, tdep2.shape[1]))
       Tvar2new = np.zeros(shape=(var2.shape[0]+2,  var2.shape[1]))
       Tmsk2new = np.zeros(shape=(tmsk2.shape[0]+2, tmsk2.shape[1]))

       Tx2new[1:-1,:]   = x2
       Tx2new[0,:]      = x2[0,:]
       Tx2new[-1,:]     = x2[-1,:]

       Tdep2new[1:-1,:] = tdep2
       Tdep2new[0,:]    = 0
       Tdep2new[-1,:]   = tdep2[-1,:]

       Tvar2new[1:-1,:] = var2
       Tvar2new[0,:]    = var2[0,:]
       Tvar2new[-1,:]   = var2[-1,:]

       if var2_aux != []:
          Tvar2_auxnew = np.zeros(shape=(var2_aux.shape[0]+2,  var2_aux.shape[1]))
          Tvar2_auxnew[1:-1,:] = var2_aux
          Tvar2_auxnew[0,:]    = var2_aux[0,:]
          Tvar2_auxnew[-1,:]   = var2_aux[-1,:]

       Tmsk2new[1:-1,:] = tmsk2
       Tmsk2new[0,:]    = tmsk2[0,:]
       Tmsk2new[-1,:]   = tmsk2[-1,:]

       #-------------------------------------------------------------
       # Interpolation on regular vertical levels if required
       if cn_RegVertLev == 'true':
          reg_x2, reg_dep2, reg_var2 = mpl_sec_int_s2z(cn_RVlevel, \
                                          Tvar2new, Tdep2new, Tx2new)
          Tx2   = np.copy(reg_x2)
          Tdep2 = np.copy(reg_dep2)
          Tvar2 = np.copy(reg_var2)
       #-------------------------------------------------------------

       else:

          Tx2   = np.copy(Tx2new)
          Tdep2 = np.copy(Tdep2new)
          Tvar2 = np.copy(Tvar2new)
          if var2_aux != []:
             Tvar2_aux = np.copy(Tvar2_auxnew)

       Tmsk2 = np.copy(Tmsk2new)

    #==============================================
    # MANAGING DEPTH for ZPS

    if vlevel == "Z_ps":
       for k in range(Tdep2.shape[0]):
           Tdep2[k,:] = Tdep2[k,int(Tdep2.shape[1]/2.)]
             
    #==============================================
    # OPENING FIGURE

    fig = plt.figure(figsize = (34.,17.5), dpi=100)

    #==============================================
    # SETTING BATHY PATCHES

    if rbat2 != [] or mbat_ln == "true":

       if msk_mes is not None:
          VLEV = 'MES'
       else:
          VLEV = vlevel

       if msk_mes is not None:
          msk_zps2 = np.ones(msk_mes2.shape) - msk_mes2
          Tdep2loc = np.copy(tdep2)
          for k in range(Tdep2loc.shape[0]):
              Tdep2loc[k,:] = Tdep2loc[k,int(Tdep2loc.shape[1]/2.)]
          rbat2_max = np.nanmax(rbat2+500.)
          rpatch_loc, mpatch_loc, RB_loc, MB_loc, ZE_loc = mpl_sec_bathy(rbat2_sec, rbat2_fill, zenv2_sec, \
                                                                         x2       , tdep2     , wdep2    , \
                                                                         tmsk2    , mbat_fill , "Z_ps"   , \
                                                                         ""       , rbat2_max, msk_zps2)
          del rpatch_loc, RB_loc, MB_loc, ZE_loc

       if PlotType == "pcolor":

             if VLEV == "Z_ps":
                rpatch, mpatch, RB, MB, ZE = mpl_sec_bathy(rbat2_sec, rbat2_fill, zenv2_sec, \
                                                           x2, tdep2, wdep2, tmsk2, mbat_fill, \
                                                           VLEV, PlotType)
             elif VLEV == "S_re" or VLEV == "MES":
                #rbat2_max = np.nanmax(rbat2)
                rbat2_max = np.nanmax(tdep3) + 0.5*(np.nanmax(tdep3[-1,:,:])-np.nanmax(tdep3[-2,:,:]))
                rpatch, mpatch, RB, MB, ZE = mpl_sec_bathy(rbat2_sec, rbat2_fill, zenv2_sec, \
                                                        x2, tdep2, Tdep2new, tmsk2, mbat_fill, \
                                                        VLEV, PlotType, rbat2_max, msk_mes2)
             else:
                rpatch, mpatch, RB, MB, ZE = mpl_sec_bathy(rbat2_sec, rbat2_fill, zenv2_sec, \
                                                           x2, tdep2, Tdep2new, tmsk2, mbat_fill, \
                                                           VLEV, PlotType)
        
       elif PlotType == "contourf":
             rbat2_max = np.nanmax(rbat2+500.)
             rpatch, mpatch, RB, MB, ZE = mpl_sec_bathy(rbat2_sec, rbat2_fill, zenv2_sec, \
                                                        x2, tdep2, wdep2, tmsk2, mbat_fill, \
                                                        VLEV, PlotType, rbat2_max, msk_mes2)

    #==============================================
    # SETTING PLOT LIMITS

    vlim_min, vlim_max, xlim_min, xlim_max, ylim_min, ylim_max = \
    mpl_sec_settings(Tvar2, Tx2, Tdep2, varlim, xlim, ylim)
  
    #==============================================
    # DRAWING SCALAR FIELD IF REQUESTED

    patch_msk = False

    if PlotType == "contourf":
       if mbat_fill == "true": #and vlevel != "Z_ps":

          #cmap, pc = mpl_sec_draw_scalar(PlotType, Tvar2, Tx2, Tdep2, vlim_min, vlim_max, \
          #                               colmap, cn_line, cn_level, cn_label, cn_color, cn_fill)
          #patch_msk = True
       #else:
          Tvar2_sol = np.copy(Tvar2)
          Tvar2_auxsol = np.copy(Tvar2_aux)
          for i in np.arange(Tvar2_sol.shape[1]):
              for k in np.arange(1,Tvar2_sol.shape[0]):
                  if Tmsk2[k,i] == 0 and Tmsk2[k-1,i] == 1:
                     Tvar2_sol[k::,i] = Tvar2_sol[k-1,i]
                     if Tvar2_auxsol != []: 
                        Tvar2_auxsol[k::,i] = Tvar2_auxsol[k-1,i]
                     break
          cmap, pc = mpl_sec_draw_scalar(PlotType, Tvar2_sol, Tx2, Tdep2, vlim_min, vlim_max, \
                                         colmap, cn_line, cn_level, cn_label, cn_color, cn_fill, Tvar2_auxsol)
               
    elif PlotType == "pcolor":
          cmap, pc = mpl_sec_draw_scalar(PlotType, Tvar2, Tx2, Tdep2, vlim_min, vlim_max, \
                                         colmap, cn_line, cn_level, cn_label, cn_color, cn_fill)
    else:
       cmap = []
       pc   = []
    
    #==============================================
    # DRAWING VERTICAL LEVELS IF REQUESTED

    #if vlevel != "no":
    #   if msk_mes is not None:
    #      wlev2_loc = np.copy(wlev2)
    #      wlev2_loc[msk_zps2 == 0] = np.nan
    #      for k in range(nk):
    #          wlev_loc = wlev2_loc[k,:]
    #          wlev_loc = wlev_loc[np.isfinite(wlev_loc)]
    #          if wlev_loc.shape[0] > 0:
    #             wlev2_loc[k,:] = np.unique(wlev_loc)[-1]
    #      wlev2[msk_zps2 == 1] = wlev2_loc[msk_zps2 == 1]
    #   W_l = mpl_sec_draw_levels(x2, wlev2)

    #==============================================
    # DRAWING GRID POINTS IF REQUESTED

    if check == "true":
       if check_val == 'true': 
          mpl_sec_draw_grid_pnts(x2, tpnt2, \
                                 var2=var2, vlim_min=vlim_min, vlim_max=vlim_max, cmap=cmap)
       else:
          mpl_sec_draw_grid_pnts(x2, tpnt2)

    #==============================================
    # DRAWING X-GRID IF REQUESTED
 
    if xgrid == "true":
       if PlotType == "contourf": 
          x = Tx2 + 0.5*(Tx2[0,1] - Tx2[0,0])
       else:
          x = Tx2
       for i in np.arange(x.shape[1]):
           X = plt.plot(x[:,i], Tdep2[:,i], linestyle='--', color='k')
           plt.setp(X, 'linewidth', 0.5)

    #==============================================
    # APPLYING PATCH IF REQUESTED

    title  = exp_name + ", " + date +  "  " + time + "\n" +\
             "Transect along " + sec_name 

    minval  = "%02e" % np.nanmin(Tvar2)
    maxval  = "%02e" % np.nanmax(Tvar2)

    meanval = "%02e" % np.mean(np.ma.masked_array(Tvar2, np.isnan(Tvar2)))
    #meanval = "%02e" % mv.filled(np.nan)

    title = title + "\n" +\
            "Range = [" + minval + ", " + maxval + "], " + \
            "Average = " + meanval

    ylabel = "Depth [m]" 
    clabel = var_strng + " " + unit_strng

    if PlotType == "contourf": 
       ax, cb = mpl_sec_figure(fig, pc, clabel, title, xlabel, ylabel, cn_level)
    else:     
       ax, cb = mpl_sec_figure(fig, pc, clabel, title, xlabel, ylabel)

    if mpatch != []: ax.add_patch(mpatch)
    if rpatch != []: ax.add_patch(rpatch)
    if patch_msk: plt.gca().patch.set_color('silver') 
    if msk_mes is not None: ax.add_patch(mpatch_loc)

    #==============================================
    # SETTING LIMITS

    plt.axis([xlim_min, xlim_max, ylim_max, ylim_min])
    pos1 = ax.get_position()

    #==============================================
    # INSET MAP

    if m != []:
       a = plt.axes([pos1.x0+0.005, pos1.y0+0.01, .2, .2])
       m.drawcoastlines()
       m.fillcontinents()
       x2, y2 = m(np.asarray(lon_sec), np.asarray(lat_sec))
       MAP = plt.plot(x2, y2, c="red")
       plt.setp(MAP, 'linewidth', 2.5)

    #==============================================
    # SAVING FIGURE

    name_fig = funcID + var_strng + "_" + sec_name + "_" +\
               PlotType + "_maxdepth_" + str(ylim_max) + "_" + timestep + fig_type

    fig.savefig(name_fig,bbox_inches="tight", pad_inches=0)

    plt.close()

    return

#===================================================================================================

def mpl_sec_loop(exp_name, fig_type, var_strng, unit_strng, date, timeres, timestep, PlotType, \
                 sec_i, sec_j, tlon3, tlat3, tdep3, wdep3, tmsk3, var4, proj, coord_type, vlevel, \
                 rbat2=[], zenv2=[], rbat2_fill="false", mbat_ln="false", mbat_fill="false", \
                 xlim='no', ylim='no', varlim='no', check='false', check_val='false', \
                 xgrid='false', colmap='jet', cn_level=20, cn_line='true', cn_label='true', \
                 cn_color='false', cn_fill='true', cn_RegVertLev='false', cn_RVlevel='auto', 
                 var_aux=None, msk_mes=None, first_zlv=None):

    '''
    This function loops over timesteps and section list
    to create section plots.
    '''

    if var4 == []:
       print("")
       print("=========================================")
       print("NO VARIABLES SELECTED")
       print("")
       plot_geom = True
    else:
       plot_geom = False
       print("")
       print("=========================================")
       print("VARIABLES: ", var_strng)
       print("")
       print(PlotType)
       print("")
       # MANAGING TIME-STEPS -----------------------------------
       if timestep == []:
          timestep = "init"
       else:
          if isinstance(timestep[0], (str, bytes)):
             if timestep[0] == "all":
                if var4 != []:
                   timestep = np.arange(var4.shape[0])
                else:
                   print("mpl_plot_sec_loop() ERROR:")
                   print("No variable selected, please check!!")
                   return
       # MANAGING LIMITS ---------------------------------------
       if isinstance(varlim, (str, bytes)):
          if varlim == "maxmin":
             varmax = np.nanmax(var4)
             varmin = np.nanmin(var4)
             print("var4 max:", varmax)
             print("var4 min:", varmin)
             varlim = [varmin, varmax]

       else:
          print("var4 max:", varlim[1])
          print("var4 min:", varlim[0])

    #===================================================================================
    # LOOP over SECTIONS

    for s in np.arange(len(sec_i)):
        I = sec_i[s]
        J = sec_j[s]
        secI = []
        secJ = []
        #print("SEC I INDEX: ", I)
        #print("SEC J INDEX: ", J)
        if isinstance(I, list) and isinstance(J, list):
           for p in np.arange(len(I)):
               if not isinstance(I[p], (int)):
                  j, i = utl.get_ij_from_lon_lat(I[p], J[p], tlon3[0,:,:], tlat3[0,:,:])
                  secI.append(i)
                  secJ.append(j)
               else:
                  secI.append(I[p])
                  secJ.append(J[p])
           #print(secI)
           #print(secJ)
        else:
           if (I == -1 and J == -1) or \
              (I != -1 and J != -1):
              print("mpl_plot_sec_loop() ERROR:")
              print("No fixed index for section chosen: please check!!")
              return

           if I == -1:
              secI = I
              if not isinstance(J, (int)):
                 j, i = utl.get_ij_from_lon_lat(tlon3[0,0,0], J, tlon3[0,:,:], tlat3[0,:,:])
                 #print(j)
                 secJ = j          
              else:
                 secJ = J
           if J == -1:
              secJ = J            
              if not isinstance(I, (int)):
                 j, i = utl.get_ij_from_lon_lat(I, tlat3[0,0,0], tlon3[0,:,:], tlat3[0,:,:])
                 #print(j, i)
                 secI = i
              else:
                 secI = I       

        if plot_geom:

           mpl_lev(exp_name, fig_type, secI, secJ, tlon3, tlat3, tdep3, wdep3, tmsk3, proj, \
                   coord_type, vlevel, check, rbat2, zenv2, rbat2_fill, mbat_ln, mbat_fill, \
                   xlim, ylim, xgrid, msk_mes, first_zlv)

        else:
           # LOOP over TIMESTEPS
           if isinstance(timestep, (str, bytes)):
              time = "00:00"
              # in this case var4 is 3D
              mpl_sec(exp_name, fig_type, var_strng, unit_strng, date, time, timestep, PlotType, \
                      secI, secJ, tlon3, tlat3, tdep3, wdep3, tmsk3, var4, proj, coord_type, \
                      vlevel, rbat2, zenv2, rbat2_fill, mbat_ln, mbat_fill, \
                      xlim, ylim, varlim, colmap, check, check_val, xgrid, \
                      cn_level, cn_line, cn_label, cn_color, cn_fill, \
                      cn_RegVertLev, cn_RVlevel, var_aux)
           else:
              for t in timestep:
                  print("time step: ", t)
                  if var4 != []:
                     var3              = var4[t,:,:,:]
                     var3[tmsk3 == 0.] = np.nan
                  else:
                     var3 = []
                  if var_aux is not None:
                     var3_aux = var_aux[t,:,:,:]
                     var3_aux[tmsk3 == 0.] = np.nan
                  else:
                     var3_aux = var_aux

                  if   timeres == "1h": 
                       time = "%02d" % t + ":30"
                  elif timeres == "1d": 
                       time = "12:30" 
                  elif timeres == "6h": 
                       time = "%02d" % 3*(2*t+1)
                  elif timeres == "4h": # for jmmp only one file output
                       time = "" #"%02d" % t*4  

                  mpl_sec(exp_name, fig_type, var_strng, unit_strng, date, time, str(t+1), PlotType,
                          secI, secJ, tlon3, tlat3, tdep3, wdep3, tmsk3, var3, proj, coord_type,
                          vlevel, rbat2, zenv2, rbat2_fill, mbat_ln, mbat_fill,
                          xlim, ylim, varlim, colmap, check, check_val, xgrid, cn_level,
                          cn_line, cn_label, cn_color, cn_fill, cn_RegVertLev, cn_RVlevel, 
                          var3_aux, msk_mes, first_zlv)

    return


