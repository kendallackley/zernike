Zernike Pipeline
=========

A Zernike Pipeline for the [GOTO](https://goto-observatory.org/) project.

See [the Wiki](https://github.com/GOTO-OBS/zernike/wiki) for a list of usage cases.


Python usage example
---------------------

````python
from zernike import zernike
from astropy.io import fits
import matplotlib.pyplot as plt
path = '/Path/To/Project/Dir/2019-08-14/'
inpath = path
catpath= path+'cat/'
outpath = path+'out/'
psfpath = path+'out/'
transpath = path+'out/'
parampath = '/Path/to/git/package/zernike/zernike/param/'
config = parampath + 'param_file.ini'

inimage = 'r0175104_UT1-median.fits'
refimage = 'r0175104_UT1-median.fits'

z = zernike.Shapelet(inimage=inimage,refimage=refimage,config=config,
                     inpath=inpath,parampath=parampath,psfpath=psfpath,
                     catpath=catpath,transpath=transpath,outpath=outpath,
                     nthreads=1,verbose=False, align=0,wcs_flag=0,subtract=1)

z.main()

### If align=0 and wcs_flag=0, then you must have the aligned template in the inpath.
### If align=1 and wcs_flag=0, then spalipy will run 
### If align=1 and wcs_flag=1, then iraf.wregister will run

### If subtract=1 then HOTPANTS will run, else the subtracted file must be in the outpath. 

### I suggest to have the script subtract for you

````
