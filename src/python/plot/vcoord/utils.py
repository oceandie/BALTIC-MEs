#!/usr/bin/env python

#     |------------------------------------------------------------|
#     | Author: Diego Bruciaferri                                  |
#     | Date and place: 07-09-2021, Met Office, UK                 |
#     |------------------------------------------------------------|


import os
from os.path import join, isfile, basename, splitext
import glob
import numpy as np
import xarray as xr
from scipy.interpolate import interp1d

# =====================================================================================================

def compute_masks(ds_domain, merge=False):
    """
    Compute masks from domain_cfg Dataset.
    If merge=True, merge with the input dataset.
    Parameters
    ----------
    ds_domain: xr.Dataset
        domain_cfg datatset
    add: bool
        if True, merge with ds_domain
    Returns
    -------
    ds_mask: xr.Dataset
        dataset with masks
    """

    # Extract variables
    k = ds_domain["z_c"] + 1
    top_level = ds_domain["top_level"]
    bottom_level = ds_domain["bottom_level"]

    # Page 27 NEMO book.
    # I think there's a typo though.
    # It should be:
    #                  | 0 if k < top_level(i, j)
    # tmask(i, j, k) = | 1 if top_level(i, j) ≤ k ≤ bottom_level(i, j)
    #                  | 0 if k > bottom_level(i, j)
    tmask = xr.where(np.logical_or(k < top_level, k > bottom_level), 0, np.nan)
    tmask = xr.where(np.logical_and(bottom_level >= k, top_level <= k), 1, tmask)
    tmask = tmask.rename("tmask")

    tmask = tmask.transpose("z_c","y_c","x_c")

    # Need to shift and replace last row/colum with tmask
    # umask(i, j, k) = tmask(i, j, k) ∗ tmask(i + 1, j, k)
    umask = tmask.rolling(x_c=2).prod().shift(x_c=-1)
    umask = umask.where(umask.notnull(), tmask)
    umask = umask.rename("umask")

    # vmask(i, j, k) = tmask(i, j, k) ∗ tmask(i, j + 1, k)
    vmask = tmask.rolling(y_c=2).prod().shift(y_c=-1)
    vmask = vmask.where(vmask.notnull(), tmask)
    vmask = vmask.rename("vmask")

    # Return
    masks = xr.merge([tmask, umask, vmask])
    if merge:
        return xr.merge([ds_domain, masks])
    else:
        return masks

# =====================================================================================================

def hvrsn_dst(lon1, lat1, lon2, lat2):
    '''
    This function calculates the great-circle distance in meters between 
    point1 (lon1,lat1) and point2 (lon2,lat2) using the Haversine formula 
    on a spherical earth of radius 6378.137km. 

    The great-circle distance is the shortest distance over the earth's surface.
    ( see http://www.movable-type.co.uk/scripts/latlong.html)

    --------------------------------------------------------------------------------

    If lon2 and lat2 are 2D matrixes, then dist will be a 2D matrix of distances 
    between all the points in the 2D field and point(lon1,lat1).

    If lon1, lat1, lon2 and lat2 are vectors of size N dist wil be a vector of
    size N of distances between each pair of points (lon1(i),lat1(i)) and 
    (lon2(i),lat2(i)), with 0 => i > N .
    '''

    deg2rad = np.pi / 180.
#    ER = 6378.137 * 1000. # Earth Radius in meters
    ER = 6372.8 * 1000. # Earth Radius in meters

    dlon = np.multiply(deg2rad, (lon2 - lon1))
    dlat = np.multiply(deg2rad, (lat2 - lat1))

    lat1 = np.multiply(deg2rad, lat1)
    lat2 = np.multiply(deg2rad, lat2)

    # Computing the square of half the chord length between the points:
    a = np.power(np.sin(np.divide(dlat, 2.)),2) + \
        np.multiply(np.multiply(np.cos(lat1),np.cos(lat2)),np.power(np.sin(np.divide(dlon, 2.)),2))

    # Computing the angular distance in radians between the points
    angle = np.multiply(2., np.arctan2(np.sqrt(a), np.sqrt(1. -a)))

    # Computing the distance 
    dist = np.multiply(ER, angle)

    return dist

# =====================================================================================================

def get_ij_from_lon_lat(lon, lat, glamt, gphit):
    '''
    get_ij_from_lon_lat find closest model grid point i/j from given lat/lon

    Syntax:
      [i, j = get_ij_from_lon_lat(lon, lat, glamt, gphit)
   
    Description:
      returns the i,j model grid position which is closest to the given
      lat/lon point. The model grid is given by 2D lon and lat matrix glamt
      and gphit.
    '''

    dist = hvrsn_dst(lon, lat, glamt, gphit)

    min_dist = np.amin(dist)

    find_min = np.where(dist == min_dist)
    sort_j = np.argsort(find_min[1])

    j_indx = find_min[0][sort_j]
    i_indx = find_min[1][sort_j]

    return j_indx[0], i_indx[0]

# =====================================================================================================

def bresenham_line(x0, x1, y0, y1):
    '''
    point0 = (y0, x0), point1 = (y1, x1)

    It determines the points of an n-dimensional raster that should be 
    selected in order to form a close approximation to a straight line 
    between two points. Taken from the generalised algotihm on

    http://en.wikipedia.org/wiki/Bresenham%27s_line_algorithm
    '''
    steep = abs(y1 - y0) > abs(x1 - x0)
   
    if steep:
       # swap(x0, y0)
       t  = y0
       y0 = x0
       x0 = t
       # swap(x1, y1)    
       t  = y1
       y1 = x1
       x1 = t

    if x0 > x1:
       # swap(x0, x1)
       t  = x1
       x1 = x0
       x0 = t
       # swap(y0, y1)
       t  = y1
       y1 = y0
       y0 = t

    deltax = np.fix(x1 - x0)
    deltay = np.fix(abs(y1 - y0))
    error  = 0.0

    deltaerr = deltay / deltax
    y = y0

    if y0 < y1:
       ystep = 1
    else:
       ystep = -1

    c=0
    pi = np.zeros(shape=[x1-x0+1])
    pj = np.zeros(shape=[x1-x0+1])
    for x in np.arange(x0,x1+1) :
        if steep:
           pi[c]=y
           pj[c]=x
        else:
           pi[c]=x
           pj[c]=y
        error = error + deltaerr
        if error >= 0.5:
           y = y + ystep
           error = error - 1.0
        c += 1

    return pj, pi

# =====================================================================================================

def get_poly_line_ij(points_i, points_j):
    '''
    get_poly_line_ij draw rasterised line between vector-points
    
    Description:
    get_poly_line_ij takes a list of points (specified by 
    pairs of indexes i,j) and draws connecting lines between them 
    using the Bresenham line-drawing algorithm.
    
    Syntax:
    line_i, line_j = get_poly_line_ij(points_i, points_i)
    
    Input:
    points_i, points_j: vectors of equal length of pairs of i, j
                        coordinates that define the line or polyline. The
                        points will be connected in the order they're given
                        in these vectors. 
    Output:
    line_i, line_j: vectors of the same length as the points-vectors
                    giving the i,j coordinates of the points on the
                    rasterised lines. 
    '''
    line_i=[]
    line_j=[]

    line_n=0

    if len(points_i) == 1:
       line_i = points_i
       line_j = points_j
    else:
       for fi in np.arange(len(points_i)-1):
           # start point of line
           i1 = points_i[fi]
           j1 = points_j[fi]
           # end point of line
           i2 = points_i[fi+1]
           j2 = points_j[fi+1]
           # 'draw' line from i1,j1 to i2,j2
           pj, pi = bresenham_line(i1,i2,j1,j2)
           if pi[0] != i1 or pj[0] != j1:
              # beginning of line doesn't match end point, 
              # so we flip both vectors
              pi = np.flipud(pi)
              pj = np.flipud(pj)

           plen = len(pi)

           for PI in np.arange(plen):
               line_n = PI
               if len(line_i) == 0 or line_i[line_n-1] != pi[PI] or line_j[line_n-1] != pj[PI]:
                  line_i.append(int(pi[PI]))
                  line_j.append(int(pj[PI]))


    return line_j, line_i

# =====================================================================================================

def create_model_bathy_sec(vlevel, msk2, xcoord2, tdep2, wdep2, max_dep=[]):
    ''' 
    This function returns the vertexs which 
    define the model bathymetry. These vertexes 
    can be used to create a patch of the model 
    bathymetry.

    msk2: 2D matrix, slice of tmask.
          N.B this matrix has to be modified
              so that
                       land  = 1.
                       ocean = np.nan
    xcoord2: 2D matrix, x coordinates of the slice
    tdep2:   2D matrix, slice of gdept_0
    wdep2:   2D matrix, slice of gdepw_0
    vlevel:  type of model vertical geometry
    '''
    nz = msk2.shape[0]
    nx = msk2.shape[1]

    z_indx = []

    for i in np.arange(nx):
        for k in np.arange(1,nz):
            if k == 1 and msk2[k-1,i] == 1:
               #print tdep2[k-1,i], wdep2[k-1,i]
               #print i, k-1
               z_indx.append(k-1)
               break
            else:
               if msk2[k,i] == 1 and np.isnan(msk2[k-1,i]):
                  #print tdep2[k,i], wdep2[k,i] 
                  #print i, k
                  z_indx.append(k)
                  break
    vert_x = []
    vert_z = []

    if vlevel == "Z_fs" or vlevel == "MES" or vlevel == "SZT":
    # For MES this is not very accurate: 
    # if levels are very wide, the plotting 
    # is very inaccurate.
       for i in range(nx):
           if i == 0:
              left_z  = wdep2[z_indx[i],i]
              right_z = wdep2[z_indx[i],i+1]
              left_x  = xcoord2[z_indx[i],i]
              right_x = xcoord2[z_indx[i],i] + 0.5*(xcoord2[z_indx[i],i+1] - xcoord2[z_indx[i],i])
           elif i > 0 and i < nx-1:
              left_z  = wdep2[z_indx[i],i]
              right_z = wdep2[z_indx[i],i+1]
              left_x  = xcoord2[z_indx[i],i] - 0.5*(xcoord2[z_indx[i],i] - xcoord2[z_indx[i],i-1])
              right_x = xcoord2[z_indx[i],i] + 0.5*(xcoord2[z_indx[i],i+1] - xcoord2[z_indx[i],i])
           else:
              left_z  = wdep2[z_indx[i],i]
              right_z = wdep2[z_indx[i],i]
              left_x   = xcoord2[z_indx[i],i] - 0.5*(xcoord2[z_indx[i],i] - xcoord2[z_indx[i],i-1])
              right_x  = xcoord2[z_indx[i],i]
           vert_z.append(left_z)
           vert_z.append(right_z)
           vert_x.append(left_x)
           vert_x.append(right_x)
    else:
       for i in range(nx):
           vert_z.append(wdep2[z_indx[i],i])
           vert_z.append(wdep2[z_indx[i],i])
           if i == 0:
              left_x  = xcoord2[z_indx[i],i]
              right_x = xcoord2[z_indx[i],i] + 0.5*(xcoord2[z_indx[i],i+1] - xcoord2[z_indx[i],i])
           elif i > 0 and i < nx-1:
              left_x  = xcoord2[z_indx[i],i] - 0.5*(xcoord2[z_indx[i],i] - xcoord2[z_indx[i],i-1])
              right_x = xcoord2[z_indx[i],i] + 0.5*(xcoord2[z_indx[i],i+1] - xcoord2[z_indx[i],i])
           else:
              left_x   = xcoord2[z_indx[i],i] - 0.5*(xcoord2[z_indx[i],i] - xcoord2[z_indx[i],i-1])
              right_x  = xcoord2[z_indx[i],i]
           vert_x.append(left_x)
           vert_x.append(right_x)

    # Appending the borders of the plot to close the patch
    vert_x.append(xcoord2[-1,-1])
    vert_x.append(xcoord2[-1,0])

    if vlevel == "MES":
       if max_dep == []:
          print("create_model_bathy_sec() ERROR:")
          print("With regular sigma levels you must provide")
          print("the maximum depth of the basin")
          return
       else:
          vert_z.append(max_dep)
          vert_z.append(max_dep)
    else:
       vert_z.append(tdep2[-1,-1] + 0.5*(tdep2[-1,-1]-tdep2[-2,-1]))
       vert_z.append(tdep2[-1,0] + 0.5*(tdep2[-1,0]-tdep2[-2,0]))

    return vert_x, vert_z
 
# =====================================================================================

def nemo_vgeom(vlevel, tdep3, wdep3):

    '''
    This function returns depth matrixes needed
    to plot correctly vertical levels and T points
    of NEMO model in sections plots. 
    '''

    ni = tdep3.shape[2]   # number of grid points along i-direction
    nj = tdep3.shape[1]   # number of grid points along j-direction
    nk = tdep3.shape[0]

    # Z PARTIAL STEPS
    if vlevel == "Z_ps":

       tpnt3 = np.copy(tdep3)
       wlev3 = np.zeros(shape=(nk,nj,ni))

       for k in np.arange(nk):

           wlev3[k,:,:] = np.unique(wdep3[k,:,:])[-1]

    # Z FULL STEPS or S / MES COORDINATES
    else:

       tpnt3 = np.copy(tdep3)
       wlev3 = np.copy(wdep3)

    return tpnt3, wlev3

# =====================================================================================

def nemo_ver_Tgrid4pcolor(var, xcoord, tgdep, wgdep):

    if var.ndim != 2:
       print("ERROR: check the dimensions of var variable!")
       return
    else:
       var2 = var

    if xcoord.ndim != 2:
       print("ERROR: check the dimensions of xcoord variable!")
       return

    if tgdep.ndim != 2:
       print("ERROR: check the dimensions of tgdep variable!")
       return

    if wgdep.ndim != 2:
       print("ERROR: check the dimensions of wgdep variable!")
       return

    # ==================
    #     xcoord 
    # ==================

    indx       = np.arange(1,xcoord.shape[1]+1) - 0.5
    indx_new   = np.arange(xcoord.shape[1]+1)
    xcoord_new = np.zeros(shape=(xcoord.shape[0]+1,len(indx_new)))

    for k in range(xcoord_new.shape[0]-1):
        f_x                = interp1d(indx, xcoord[k,:])
        xcoord_new[k,1:-1] = f_x(indx_new[1:-1])
        xcoord_new[k,0]    = xcoord[k,0] - (xcoord_new[k,1] - xcoord[k,0])
        xcoord_new[k,-1]   = xcoord[k,-1] + (xcoord[k,-1] - xcoord_new[k,-2])
    xcoord_new[-1,:]       = xcoord_new[-2,:]

    # ==================
    #      depth
    # ==================

    depth      = np.zeros(shape=(wgdep.shape[0]+1,len(indx_new)))

    wgdep1         = np.zeros(shape=(wgdep.shape[0]+1,len(indx)))
    wgdep1[:-1,:]  = wgdep
    wgdep1[-1,:]   = tgdep[-1,:] + (tgdep[-1,:] - wgdep[-1,:])

    xcoord1        = np.zeros(shape=(xcoord.shape[0]+1,len(indx)))
    xcoord1[:-1,:] = xcoord
    xcoord1[-1,:]  = xcoord[-1,:]

    # CREATING arrays for interpolation
    wgdep1_vec      = np.zeros(shape=(wgdep1.shape[0]*(wgdep1.shape[1]-1), 2 ))
    wgdep1_vec[:,0] = np.ravel(wgdep1[:,:-1], order='F')
    wgdep1_vec[:,1] = np.ravel(wgdep1[:,1:], order='F')

    xcoord1_vec      = np.zeros(shape=(wgdep1.shape[0]*(wgdep1.shape[1]-1), 2 ))
    xcoord1_vec[:,0] = np.ravel(xcoord1[:,:-1], order='F')
    xcoord1_vec[:,1] = np.ravel(xcoord1[:,1:], order='F')

    xcoord_4int     = xcoord_new[:,1:-1]
    xcoord_4int_vec = np.ravel(xcoord_4int, order='F')

    # LINEAR INTERPOLATION of wgdep on xcoord_new
    wgdep_int = wgdep1_vec[:,0] + \
                np.multiply( wgdep1_vec[:,1]-wgdep1_vec[:,0],\
                             np.divide(xcoord_4int_vec-xcoord1_vec[:,0],\
                                       xcoord1_vec[:,1]-xcoord1_vec[:,0]) 
                           )

    depth[:,1:-1] = np.reshape(wgdep_int, (wgdep1.shape[0],wgdep1.shape[1]-1), order='F')

    # EXTRAPOLATING depth values for left and right boundaries
    depth[:,0]    = depth[:,1]
    depth[:,-1]    = depth[:,-2]

    # ==================
    #     var
    # ==================

    mskval = 1.e-20

    var_new = np.ones(shape=(wgdep.shape[0]+1,len(indx_new))) * mskval
    var_new[:-1,:-1] = var
    var_new[var_new == mskval] = np.nan


    return var_new, xcoord_new, depth

