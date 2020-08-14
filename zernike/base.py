import numpy as np
from astropy.table import Table
import subprocess

__all__ = ["zernike"]


class ZernObj():
    def __init__(self):
        return

    def _setatt__ (self,z_ord,z_mean,z_mean_std,z_med,z_med_std,z_gauss,
                    z_gauss_std,z_gauss_chi2):
        self.z_ord = [z_ord]
        self.z_mean = [z_mean]
        self.z_mean_std = [z_mean_std]
        self.z_med = [z_med]
        self.z_med_std = [z_med_std]
        self.z_gauss = [z_gauss]
        self.z_gauss_std = [z_gauss_std]
        self.z_gauss_chi2 = [z_gauss_chi2]

class TransObj():
    def __init__(self):
        return

    def set_att(self,x,y,snr,zdist,mag_inst,mag_inst_err,ra,dec,flux,flux_err,
                    mag_app,mag_app_err,seeing):
        self.x = x
        self.y = y
        self.snr = snr
        self.zdist = zdist
        self.mag_inst = mag_inst
        self.mag_inst_err = mag_inst_err
        self.ra = ra
        self.dec = dec
        self.flux = flux
        self.flux_err = flux_err
        self.mag_app = mag_app
        self.mag_app_err = mag_app_err
        self.seeing = seeing


class FileIO():

    def read_params_file(self, params_file, **kwargs):

        params_dict={}
        with open(params_file,'r') as params_f:
            for line in params_f:
                s=line.split()
                if s:
                    if s[0][0] != '#':
                        params_dict[s[0]]=s[1]
        return params_dict

    def read_file(self, filename, **kwargs):
        from astropy.io import ascii

        """
        Read in an ascii file
        """
        try:
            x_key = kwargs['x_key']
            y_key = kwargs['y_key']
        except:
            x_key = 'XWIN_IMAGE'
            y_key = 'YWIN_IMAGE'
        cat = ascii.read(filename)
        cat[x_key] -= 1
        cat[y_key] -= 1
        return cat

    def read_zern_file(self,filename, **kwargs):
        from astropy.io import ascii

        """
        Read in an ascii file
        """
        names = ['Z_ORD','MEAN','MEAN_STD','MED','MED_STD','GAUSS','GAUSS_STD',
                'CHI2']
        return ascii.read(filename,names=names)

    def write_file(self, filename, **kwargs):
        from astropy.io import ascii

        """
        Write an ascii file
        """
        return ascii.SExtractor.write(filename)

    def sextractor_script(self, infile, catfile, sex_param, align=False, **kwargs):
        """
        """
        import subprocess
        import shutil

        sargs = []
        for kw in kwargs:
            sargs.append(" -{} {} ".format(kw,kwargs[kw]))

        cmd = shutil.which('sextractor')
        if cmd is not None:
            pass
        elif cmd is None:
            cmd = shutil.which('sex')
            if cmd is not None:
                pass
            else:
                warnings.warn("No SExtractor installed")

        if align:
            sargs.append(" -{} {} ".format('CATALOG_TYPE','FITS_LDAC'))
            sargs.append(" -{} {} ".format('c',sex_param))
            sargs.append(" -{} {} ".format('CATALOG_NAME',infile.replace('.fits','.ldac.fits')))

        else:
            sargs.append(" -{} {} ".format('CATALOG_TYPE','ASCII_HEAD'))
            sargs.append(" -{} {} ".format('c',sex_param))
            sargs.append(" -{} {} ".format('CATALOG_NAME',catfile))



        args = [cmd+' '+infile+' '.join(sargs)]
        print(args)
        subprocess.call(args,shell=True)
        return

    # def sextractor_script(self,infile,catfile,sex_param,align=False,**kwargs):
    #     """
    #     """
    #     import astromatic_wrapper as aw
    #     import shutil
    #
    #     sargs = []
    #     for kw in kwargs:
    #         sargs.append(" -{} {} ".format(kw,kwargs[kw]))
    #
    #     files = {'image': infile}
    #
    #     cmd = shutil.which('sextractor')
    #     if cmd is not None:
    #         kwargs['cmd'] = cmd
    #     elif cmd is None
    #         cmd = shutil.which('sex')
    #         if cmd is not None:
    #             kwargs['cmd'] = cmd
    #         else:
    #             warnings.warn("No SExtractor installed")
    #
    #     if align:
    #         kwargs['config']['CATALOG_TYPE'] = 'FITS_LDAC'
    #         kwargs['config_file'] = sex_swarp_param
    #         sextractor = aw.api.Astromatic(**kwargs)
    #         sextractor.run_frames(files['image'],frames=[1])
    #
    #     else:
    #         kwargs['config']['CATALOG_TYPE'] = 'ASCII_HEAD'
    #         kwargs['config_file'] = sex_param
    #
    #
    #
    #
    #     # args = ['sex '+infile+' -c '+sex_param+' -CATALOG_NAME '+catfile+
    #     #         ' '.join(sargs)]
    #     # subprocess.call(args,shell=True)
    #     return

    def remap_swarp(self,inimage,refimage,output,**kwargs):
        swarp_kwargs = []
        for kw in kwargs:
            swarp_args.append(" -{} {} ".format(kw,kwargs[kw]))

        swarp_args.append(" -{} {} ".format('IMAGEOUT_NAME',output))
        swarp_args.append(" -{} {} ".format('COMBINE','N'))
        swarp_args.append(" -{} {} ".format('SUBTRACT_BACK','Y'))


        args = ['swarp '+inimage+refimage+' -c '+swarp_param+''.join(swarp_args)]
        subprocess.call(args,shell=True)
        return

    def remap_wcs(self,inimage,refimage,outimage,**kwargs):
        from pyraf import iraf
        try:
            sciext = kwargs['sci']
            templext = kwargs['ref']
        except KeyError:
            sciext = 0
            templext = 0
        sci_ext = '{}'.format(sciext)
        ref_ext = '{}'.format(templext)
        return iraf.wregister(input=inimage+sci_ext, output=outimage, wcs="world",
                        reference=refimage+ref_ext, Stdout=1)

    def remap_spalipy(self,incat,refcat,inimage,output,**kwargs):
        from spalipy import spalipy
        try:
            sciext = kwargs['sci']
            templext = kwargs['ref']
        except KeyError:
            sciext = 0
            templext = 0
        sci_ext = '[{}]'.format(sciext)
        s = spalipy.Spalipy(incat, refcat, inimage,output_filename=output)
        s.main()
        return

    def run_hotpants(self,infile,alreffile,subfile,convfile,hpargs=None,**kwargs):
        """
        Additional options:
       [-tu tuthresh]    : upper valid data count, template (25000)
       [-tuk tucthresh]  : upper valid data count for kernel, template (tuthresh)
       [-tl tlthresh]    : lower valid data count, template (0)
       [-tg tgain]       : gain in template (1)
       [-tr trdnoise]    : e- readnoise in template (0)
       [-tp tpedestal]   : ADU pedestal in template (0)
       [-tni fitsfile]   : input template noise array (undef)
       [-tmi fitsfile]   : input template mask image (undef)
       [-iu iuthresh]    : upper valid data count, image (25000)
       [-iuk iucthresh]  : upper valid data count for kernel, image (iuthresh)
       [-il ilthresh]    : lower valid data count, image (0)
       [-ig igain]       : gain in image (1)
       [-ir irdnoise]    : e- readnoise in image (0)
       [-ip ipedestal]   : ADU pedestal in image (0)
       [-ini fitsfile]   : input image noise array (undef)
       [-imi fitsfile]   : input image mask image (undef)

       [-ki fitsfile]    : use kernel table in image header (undef)
       [-r rkernel]      : convolution kernel half width (10)
       [-kcs step]       : size of step for spatial convolution (2*rkernel + 1)
       [-ft fitthresh]   : RMS threshold for good centroid in kernel fit (20.0)
       [-sft scale]      : scale fitthresh by this fraction if... (0.5)
       [-nft fraction]   : this fraction of stamps are not filled (0.1)
       [-mins spread]    : Fraction of kernel half width to spread input mask (1.0)
       [-mous spread]    : Ditto output mask, negative = no diffim masking (1.0)
       [-omi  fitsfile]  : Output bad pixel mask (undef)
       [-gd xmin xmax ymin ymax]
                         : only use subsection of full image (full image)

       [-nrx xregion]    : number of image regions in x dimension (1)
       [-nry yregion]    : number of image regions in y dimension (1)
       -- OR --
       [-rf regionfile]  : ascii file with image regions 'xmin:xmax,ymin:ymax'
       -- OR --
       [-rkw keyword num]: header 'keyword[0->(num-1)]' indicates valid regions

       [-nsx xstamp]     : number of each region's stamps in x dimension (10)
       [-nsy ystamp]     : number of each region's stamps in y dimension (10)
       -- OR --
       [-ssf stampfile]  : ascii file indicating substamp centers 'x y'
       -- OR --
       [-cmp cmpfile]    : .cmp file indicating substamp centers 'x y'

       [-afssc find]     : autofind stamp centers so #=-nss when -ssf,-cmp (1)
       [-nss substamps]  : number of centroids to use for each stamp (3)
       [-rss radius]     : half width substamp to extract around each centroid (15)

       [-savexy file]    : save positions of stamps for convolution kernel (undef)
       [-c  toconvolve]  : force convolution on (t)emplate or (i)mage (undef)
       [-n  normalize]   : normalize to (t)emplate, (i)mage, or (u)nconvolved (t)
       [-fom figmerit]   : (v)ariance, (s)igma or (h)istogram convolution merit (v)
       [-sconv]          : all regions convolved in same direction (0)
       [-ko kernelorder] : spatial order of kernel variation within region (2)
       [-bgo bgorder]    : spatial order of background variation within region (1)
       [-ssig statsig]   : threshold for sigma clipping statistics  (3.0)
       [-ks badkernelsig]: high sigma rejection for bad stamps in kernel fit (2.0)
       [-kfm kerfracmask]: fraction of abs(kernel) sum for ok pixel (0.990)
       [-okn]            : rescale noise for 'ok' pixels (0)
       [-fi fill]        : value for invalid (bad) pixels (1.0e-30)
       [-fin fill]       : noise image only fillvalue (0.0e+00)
       [-convvar]        : convolve variance not noise (0)

       [-oni fitsfile]   : output noise image (undef)
       [-ond fitsfile]   : output noise scaled difference image (undef)
       [-nim]            : add noise image as layer to sub image (0)
       [-ndm]            : add noise-scaled sub image as layer to sub image (0)

       [-oci fitsfile]   : output convolved image (undef)
       [-cim]            : add convolved image as layer to sub image (0)

       [-allm]           : output all possible image layers

       [-nc]             : do not clobber output image (0)
       [-hki]            : print extensive kernel info to output image header (0)

       [-oki fitsfile]   : new fitsfile with kernel info (under)

       [-sht]            : output images 16 bitpix int, vs -32 bitpix float (0)
       [-obs bscale]     : if -sht, output image BSCALE, overrides -inim (1.0)
       [-obz bzero]      : if -sht, output image BZERO , overrides -inim (0.0)
       [-nsht]           : output noise image 16 bitpix int, vs -32 bitpix float (0)
       [-nbs bscale]     : noise image only BSCALE, overrides -obs (1.0)
       [-nbz bzero]      : noise image only BZERO,  overrides -obz (0.0)

       [-ng  ngauss degree0 sigma0 .. degreeN sigmaN]
                         : ngauss = number of gaussians which compose kernel (3)
                         : degree = degree of polynomial associated with gaussian #
                                    (6 4 2)
                         : sigma  = width of gaussian #
                                    (0.70 1.50 3.00)
                         : N = 0 .. ngauss - 1

                         : (3 6 0.70 4 1.50 2 3.00
       [-pca nk k0.fits ... n(k-1).fits]
                         : nk      = number of input basis functions
                         : k?.fits = name of fitsfile holding basis function
                         : Since this uses input basis functions, it will fix :
                         :    hwKernel
                         :
       [-v] verbosity    : level of verbosity, 0-2 (1)
        """
        hpargs = []

        for kw in kwargs.keys():
            if kw not in ['tmpl' ,'sci']:
                hpargs.append(" -{} {} ".format(kw,kwargs[kw]))
        # try:
        #     kwargs['oci']
        #     refcfile = kwargs['oci']
        # except KeyError:
        #     warnings.warn("""Wait a minute! Should run PSF on convolved image.
        #     Making and setting the output convolved image name for HOTPANTS.
        #     """)
        #     refcfile = reffile.replace('.fits','_conv.fits')
        #     hpargs.append(" -{} {} ".format('oci', refcfile))
        #     pass
        print('Running HOTPANTS')
        args = ['hotpants -inim '+infile+kwargs['sci']+' -tmplim '+alreffile+kwargs['tmpl']+' -outim '+
                subfile+ ' -oci '+convfile + ' -n i -sconv '+''.join(hpargs)] #' -oci '+convfile
        print(args,convfile)
        subprocess.call(args,shell=True)
        # VS subprocess.Popen(args).wait() for error handling
        return

    def write_psf_file(self,psf_file,img_nam,cat_file,
            wf_s_coeff_mean, wf_s_coeff_med,wf_s_coeff_gauss,flux_min,cls_min,
            j_max,sub_dim, disk_dim,disk_rad,len_wf_s_coeff_list,
            len_inj_cat_filter,wf_s_coeff_gauss_std,wf_s_coeff_gauss_chi2,wf_s_coeff_mean_std,
            wf_s_coeff_med_std):
        with open(psf_file, 'w') as psfID:
            psfID.write("# Source list compiled from "+img_nam+" using catalog file "+cat_file+".\n")
            psfID.write("# Parameters for source selection were flux_min=%5.3f and cls_min=%1.3f.\n" % (flux_min, cls_min))
            psfID.write("# %4g of %5g catalog objects Z-decomposed up to j_max=%2g.\n" % (len_wf_s_coeff_list, len_inj_cat_filter, j_max))
            psfID.write("# Dimensions of subimages are sub_dim: %2g and disk_dim: %2g.\n" % (sub_dim, disk_dim))
            psfID.write("# %2g\tSUB_DIM\n" % sub_dim)
            psfID.write("# %2g\tDISK_DIM\n" % disk_dim)
            psfID.write("# %3.1f\tDISK_RAD\n" % disk_rad)
            psfID.write("# %2g\tJ_MAX\n" % j_max)
            psfID.write("# 0 \tZ_ORD\tZERNIKE ORDER\n")
            psfID.write("# 1 \tMEAN\tMEAN OF DISTRIBUTION\n")
            psfID.write("# 2 \tMN_SD\tSTANDARD DEVIATION WRT MEAN OF DISTRIBUTION\n")
            psfID.write("# 3 \tMED\t\tMEDIAN OF DISTRIBUTION\n")
            psfID.write("# 4 \tMD_SD\tSTANDARD DEVIATION WRT MEDIAN OF DISTRIBUTION\n")
            psfID.write("# 5 \tGAUSS\tCENTER OF GAUSSIAN FIT\n")
            psfID.write("# 6 \tGS_SD\tSTANDARD DEVIATION OF GAUSSIAN FIT\n")
            psfID.write("# 7 \tCHI2\tCHI^2 OF GAUSSIAN FIT\n")
            for k in range(len(wf_s_coeff_mean)):
                psfID.write("%2g" % k)
                psfID.write("\t%+10.8f" % wf_s_coeff_mean[k])
                psfID.write("\t%+10.8f" % wf_s_coeff_mean_std[k])
                psfID.write("\t%+10.8f" % wf_s_coeff_med[k])
                psfID.write("\t%+10.8f" % wf_s_coeff_med_std[k])
                psfID.write("\t%+10.8f" % wf_s_coeff_gauss[k])
                psfID.write("\t%+10.8f" % wf_s_coeff_gauss_std[k])
                psfID.write("\t%+10.8f" % wf_s_coeff_gauss_chi2[k])
                psfID.write("\n")
            psfID.close()
        return

    def write_trans_file(self,trans_file,trans_list,len_sub_cat,dz_max,filter_trans=0,filter_dz=0):
        print('WRITING TRANS FILE TO',trans_file)
        with open(trans_file, 'w') as transID:
            transID.write("# Zernike Distance: \t%5g\n" %dz_max)
            transID.write("# Number of objects in subtracted image: \t%5g\n" %len_sub_cat)
            transID.write("# Number of transient objects detected: \t%5g\n" %len(trans_list))
            #transID.write("# SEEING: \t%4.2f\n" %seeing)
            transID.write("# 0 \tNUMBER\n")
            transID.write("# 1 \tZ_DIST\n")
            transID.write("# 2 \tX_COORD\n")
            transID.write("# 3 \tY_COORD\n")
            transID.write("# 4 \tRA\n")
            transID.write("# 5 \tDEC\n")
            transID.write("# 6 \tSNR\n")
            transID.write("# 7 \tMAG_INST\n")
            transID.write("# 8 \tMAG_INST_ERR\n")
            transID.write("# 9 \tMAG_APP\n")
            transID.write("# 10 \tMAG_APP_ERR\n")
            transID.write("# 11 \tFLUX\n")
            transID.write("# 12 \tFLUX_ERR\n")
            transID.write("# 13 \tFWHM\n")

            for i,trans_item in enumerate(trans_list):
                x,y,snr,zdist,mag_inst,mag_inst_err,ra,dec,flux,flux_err,\
                                mag_app,mag_app_err,fwhm=trans_item
                if filter_trans and zdist > filter_dz:
                    continue

                transID.write("%4g" % i)
                transID.write("\t%6.2f"  % zdist)
                transID.write("\t%11.6f" % x)
                transID.write("\t%11.6f" % y)
                transID.write("\t%12.6f" % ra)
                transID.write("\t%12.6f" % dec)
                transID.write("\t%6.2f"  % snr)
                transID.write("\t%10.6f" % mag_inst)
                transID.write("\t%10.6f" % mag_inst_err)
                transID.write("\t%10.6f" % mag_app)
                transID.write("\t%10.6f" % mag_app_err)
                transID.write("\t%10.6f" % flux)
                transID.write("\t%10.6f" % flux_err)
                transID.write("\t%6.2f"  % fwhm)
                transID.write("\n")
        transID.close()
        return

    def deg2HMS(self,rain):

       if rain < 0:
          s = -1
          ra = -rain
       else:
          s = 1
          ra = rain

       h = int(ra/15)
       ra -= h*15.
       m = int(ra*4)
       ra -= m/4.
       s = ra*240.

       if s == -1:
          return '-%02d:%02d:%06.3f'%(h,m,s)
       else:
           return '+%02d:%02d:%06.3f'%(h,m,s)

    def deg2DMS(self,decin):

       if decin < 0:
          s = -1
          dec = -decin
       else:
          s = 1
          dec = decin

       d = int(dec)
       dec -= d
       dec *= 100.
       m = int(dec*0.6)
       dec -= m*5./3.
       s = dec*180./5.

       if s == -1:
          return '+%02d:%02d:%06.3f'%(d,m,s)
       else:
           return '-%02d:%02d:%06.3f'%(d,m,s)


    def write_reg(self,reg_file,trans_list,filter_trans=0,filter_dz=0):
        from astropy.coordinates import SkyCoord
        from astropy import units as u

        with open(reg_file, 'w') as iregID:
            iregID.write("# Region file format: DS9 version 4.1\n")
            iregID.write('global color=red dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n')
            iregID.write("fk5\n")
            for i, trans_item in enumerate(trans_list):
                x,y,snr,zdist,mag_inst,mag_inst_err,ra,dec,flux,flux_err,\
                                mag_app,mag_app_err,fwhm=trans_item
                if not filter_trans:
                    iregID.write('circle(%s,%s,17.8744")' % (ra,dec))
                    iregID.write("\n")
                else:
                    if not filter_dz:
                        iregID.write('circle(%s,%s,12.8744")' % (ra,dec))
                        iregID.write("\n")
                    else:
                        if zdist <= filter_dz:
                            iregID.write('circle(%s,%s,12.8744")' % (ra,dec))
                            iregID.write("\n")
        iregID.close()
        return
