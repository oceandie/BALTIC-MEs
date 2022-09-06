
import numpy as np
import netCDF4 as nc4
import sys
import os
import ast

def applymask(invar,bmdi):
    # apply mdi values to an array
    vsize=invar.shape
    invar=invar.flatten()
    if type(bmdi) == type([]):
        # we have a list of mdi values
        # set all to the last mdi value
        for item in bmdi:
            mdi=invar==item
            invar[mdi]=bmdi[-1]
    else:
        mdi=invar==bmdi
        invar[mdi]=bmdi
    return np.reshape(invar,vsize)


#####################################################################
## main rebuild script
#####################################################################

INFILE=sys.argv[1]
if os.path.dirname(INFILE) == '':
    # if no directory path is provided then assume $PWD
    INDIR=os.getcwd()
    INSTEM=INDIR+'/'+INFILE
else:
    INSTEM=INFILE

# check if the infile will exist (assuming the first file is INSTEM_0000.nc)
if not os.path.isfile(INSTEM+'_0000.nc'):
    print('Supplied filename does not exist:')
    print(INSTEM)
    sys.exit(1) 

try:
    domask2=str(sys.argv[2])
except:
    domask2='False'

#maskdict={'True':True,'False':False}
#domask=maskdict[domask2]
domask=ast.literal_eval(domask2)

# pass OUTDIR if necessary
POUTDIR=sys.argv[3]
if os.path.dirname(POUTDIR) == '':
    # if no directory path is provided then assume $PWD
    ODIR=os.getcwd()
    OUTDIR=ODIR+'/'+POUTDIR
else:
    OUTDIR=POUTDIR

# create output directory if it doesn't exist
if not os.path.isdir(OUTDIR):
    os.mkdir(OUTDIR)

# get OUTFILE name from INSTEM
OUTFILE=OUTDIR+'/'+os.path.basename(INSTEM)


print(INSTEM)
print(domask)
print(OUTDIR)


if domask:
 MASKSTR='_mask'
else:
 MASKSTR=''


# load first file
print('Reading '+INSTEM+'_0000.nc')
f0=nc4.Dataset(INSTEM+'_0000.nc','r')
xx,yy=f0.DOMAIN_size_global
tdoms=f0.DOMAIN_number_total


# loop over variables
vnames    = []
vdata     = []
tt        = 0
zz        = 0
vmdi      = []
exceptmdi = -99
vert_flag = []

for v in f0.variables:

        #print v, f0.variables[v].ndim, f0.variables[v].shape
        
        vnames.append(v)

        vflag = False
        #if not v.endswith('lat') and not v.endswith('lon'): vflag = True
        cond1 = not v.endswith('lat') and not v.endswith('lon')
        cond2 = not v.endswith('x')   and not v.endswith('y')
        if cond1 and cond2: 
           vflag = True
        vert_flag.append(vflag)

        # determine length of time and depth dimension
        if f0.variables[v].ndim == 2 and vflag: 
           if tt == 0: tt=f0.variables[v].shape[0]
           if zz == 0: zz=f0.variables[v].shape[1]     
        if f0.variables[v].ndim == 3:
           if tt == 0: tt=f0.variables[v].shape[0]
        if f0.variables[v].ndim == 4:
           if tt == 0: tt=f0.variables[v].shape[0]
           #if zz == 0: 
           zz=f0.variables[v].shape[1]
           
        # create final arrays
        if f0.variables[v].ndim == 0:
            vdata.append(np.zeros(1,'float'))

        if f0.variables[v].ndim == 1:
            vdata.append(np.zeros(f0.variables[v].shape[0],'float'))

        if f0.variables[v].ndim == 2:
            if vflag:
               vdata.append(np.zeros((tt,zz),'float')) 
            else:
               vdata.append(np.zeros((yy,xx),'float'))

        if f0.variables[v].ndim == 3:
            # assume time,y,x
            vdata.append(np.zeros((tt,yy,xx),'float'))

        if f0.variables[v].ndim == 4:
            
            # assume time,z,y,x
            vdata.append(np.zeros((tt,zz,yy,xx),'float'))
 
        if domask:
            # determine and check missing data values for each variable
            try:
                thismdi=f0.variables[v]._FillValue
            except:
                print('')
                print('No FillValue attribute for '+v)
                thismdi=exceptmdi
            if f0.variables[v].ndim == 2:
                mymdi=f0.variables[v][0,0]
                try:
                    mymdi2=f0.variables[v][1,-2]
                except:
                    mymdi2=mymdi
                thisv=f0.variables[v][:,:]
            if f0.variables[v].ndim == 3:
                mymdi=f0.variables[v][0,0,0]
                try:
                    mymdi2=f0.variables[v][0,1,-2]
                except:
                    mymdi2=mymdi
                thisv=f0.variables[v][0,:,:]
            if f0.variables[v].ndim == 4:
                mymdi=f0.variables[v][0,-1,0,0]
                try:
                    mymdi2=f0.variables[v][0,-1,1,-2]
                except:
                    mymdi2=mymdi
                thisv=f0.variables[v][0,-1,:,:]
            # is mymdi value masked already suggesting a masked array
            if not np.ma.is_masked(mymdi):
                # are the 2 mdi values the same
                if thismdi != mymdi and thismdi != exceptmdi:
                    # are there any of the mdi values stated in the file
                    if thisv.flatten().tolist().count(thismdi) == 0:
                        # no values so check my guesses at fill value
                        # if mymdi and mymdi2 are the same then just store check one of them otherwise check both
                        if mymdi != mymdi2:
                            percentchk=10
                            # check a special case first where more than percentchk of the data is both mdi values
                            if thisv.flatten().tolist().count(mymdi)/float(thisv.size)*100. > percentchk and thisv.flatten().tolist().count(mymdi2)/float(thisv.size)*100. > percentchk:
                                # assume both are mdi values
                                # so save both
                                vmdi.append([mymdi,mymdi2])
                            elif thisv.flatten().tolist().count(mymdi)/float(thisv.size)*100. > percentchk and thisv.flatten().tolist().count(mymdi2)/float(thisv.size)*100. < percentchk:
                                # there is more mymdi values
                                vmdi.append(mymdi)
                            elif thisv.flatten().tolist().count(mymdi2)/float(thisv.size)*100. > percentchk and thisv.flatten().tolist().count(mymdi)/float(thisv.size)*100. < percentchk:
                                # there is more mymdi2 values
                                vmdi.append(mymdi2)
                            else:
                                # as a fallback just use the value from file as both less than percentchk of the data
                                vmdi.append(thismdi)
                        else:
                            # both mdi values are the same so just check one
                            if thisv.flatten().tolist().count(mymdi) > 1:
                                vmdi.append(mymdi)
                            else:
                                vmdi.append(thismdi)
                    else:
                        # use mdi value from file
                        vmdi.append(thismdi)
                else:
                    # mdi values are the same so skip the check and just store
                    vmdi.append(thismdi)
            else:
                # mdi values are the same so skip the check and just store
                vmdi.append(thismdi)

numvars=len(vnames)

if domask:
    print('')
    print('Using these mdi values: ')
    print(vmdi)
    print('')


# load each file and save the data
for f in range(tdoms):
    if f != 0:
        fnum="%04d" % f
        print('Reading '+INSTEM+'_'+fnum+'.nc')
        f0=nc4.Dataset(INSTEM+'_'+fnum+'.nc','r')

    x0=f0.DOMAIN_position_first[0]-1
    y0=f0.DOMAIN_position_first[1]-1
    x1=f0.DOMAIN_position_last[0]
    y1=f0.DOMAIN_position_last[1]
 
    if y1 > yy:
        y1=yy

    if x1 > xx:
        x1=xx

    xx1=x1-x0
    yy1=y1-y0

    for v in range(numvars):
        #print 'Remaining variables to extract = '+str(numvars-v)
        thisv=f0.variables[vnames[v]]

        # store variable data
        if thisv.ndim == 0: vdata[v][0]=thisv[0]

        if thisv.ndim == 1: vdata[v][:]=thisv[:]

        if thisv.ndim == 2:
           # apply mask if setting is active
           if domask:
               if vmdi[v] != exceptmdi:
                  thisv=applymask(thisv[:,:],vmdi[v])
           if vert_flag[v]:
              vdata[v][:,:]=thisv[:,:]
           else:    
              vdata[v][y0:y1,x0:x1]=thisv[0:yy1,0:xx1]

        if thisv.ndim == 3:
           # apply mask if setting is active
           if domask:
              if vmdi[v] != exceptmdi:
                 thisv=applymask(thisv[:,:,:],vmdi[v])
           vdata[v][:,y0:y1,x0:x1]=thisv[:,0:yy1,0:xx1]

        if thisv.ndim == 4:
           # apply mask if setting is active
           if domask:
              if vmdi[v] != exceptmdi:
                 thisv=applymask(thisv[:,:,:,:],vmdi[v])
           vdata[v][:,:,y0:y1,x0:x1]=thisv[:,:,0:yy1,0:xx1]

    if f != tdoms-1:
        f0.close()


# output result
print('Creating output file: '+OUTFILE+MASKSTR+'.nc')
outf=nc4.Dataset(OUTFILE+MASKSTR+'.nc','w',format="NETCDF4")

# create dimensions
if sys.version_info[0] < 3:
   iterator = f0.dimensions.iteritems()
else:
   iterator = f0.dimensions.items()

for name, dimension in iterator:
    if name == 'x':
       print('Creating dimension: '+name+' size = '+str(f0.DOMAIN_size_global[0]))
       outf.createDimension(name, (f0.DOMAIN_size_global[0] if not dimension.isunlimited() else None))
    elif name == 'y':
       print('Creating dimension: '+name+' size = '+str(f0.DOMAIN_size_global[1]))
       outf.createDimension(name, (f0.DOMAIN_size_global[1] if not dimension.isunlimited() else None))
    else:
       print('Creating dimension: '+name+' size = '+str(len(dimension)))
       outf.createDimension(name, (len(dimension) if not dimension.isunlimited() else None))

# create variables
vcount=-1
if sys.version_info[0] < 3:
   iterator = f0.variables.iteritems()
else:
   iterator = f0.variables.items()

for name, variable in iterator:
    print('Outputting variable: '+name)
    if name in vnames:
        vcount+=1
        # create nc variable
        if name == 'nav_lon' or name == 'nav_lat':
            x = outf.createVariable(name, variable.datatype, variable.dimensions)   
        else:
            if domask:
                if type(vmdi[vcount]) == type([]):
                    x = outf.createVariable(name, variable.datatype, variable.dimensions, fill_value=vmdi[vcount][-1])
                else:
                    x = outf.createVariable(name, variable.datatype, variable.dimensions, fill_value=vmdi[vcount])
            else:
                try:
                    x = outf.createVariable(name, variable.datatype, variable.dimensions, fill_value=variable._FillValue)
                except:
                    x = outf.createVariable(name, variable.datatype, variable.dimensions)
        # output data
        #print x
        #print vdata[vcount].shape
        x[:]=vdata[vcount][:]
        # copy variable attributes
        x.setncatts({k: variable.getncattr(k) for k in variable.ncattrs() if k != '_FillValue'})
    else:
        if f0.variables[name].dimensions:
            x = outf.createVariable(name, variable.datatype, variable.dimensions)
            x[:]=f0.variables[name][:]
            try:
                x.setncatts({k: variable.getncattr(k) for k in variable.ncattrs()})
            except:
                x.setncatts({k: variable.getncattr(k) for k in variable.ncattrs() if k != '_FillValue'})
        else:
            x = outf.createVariable(name, variable.datatype)
            x = f0.variables[name][:]

# copy global attributes
outf.setncatts({k: f0.getncattr(k) for k in f0.ncattrs() if not k.startswith('DOMAIN') })

outf.close()
f0.close()
