
#!/usr/bin/env python3.6
from . import psf
from . import base
from . import test_conda
import argparse
import astropy.io.fits as fits
from astropy.io import ascii
from astropy import wcs
from astropy.table import Table
from astropy.stats import sigma_clip
import warnings
warnings.simplefilter('ignore', UserWarning)
import numpy as np
import os
import glob
import subprocess
from scipy import ndimage
import time
import importlib
from multiprocessing.dummy import Pool as ThreadPool
from multiprocessing import Pool
from configobj import ConfigObj
from astropy.table import Column
import sys
import matplotlib.pyplot as plt

import warnings

zern_poly = psf.ZernikePoly()
fileIO = base.FileIO()
misc = psf.Catalog()
impsf = psf.PSF()
zernobj = base.ZernObj()
trans = base.TransObj()


class Shapelet:

    def __init__(self,inimage,config,outpath=None,catpath=None,nthreads=1,
                verbose=False,del_tmp=False,make_trans=True,plot_psf=True,
                save_as_array=True):

        # if isinstance(ncpus,float):
        #     self.cpupool = Pool(ncpus)

        self.verbose = verbose
        self.del_tmp = del_tmp
        self.make_trans = make_trans
        self.plot_psf = plot_psf
        self.save_as_array = save_as_array
        # self.save_to_file = save_to_file

        if isinstance(config,str):
            try:
                self.config = ConfigObj(config,file_error=True)
            except OSError:
                print("something's wrong")
                # raise SystemExit

            if not len(self.config.sections):
                # TODO Type default values here
                warnings.warn(""""
                No config file found. Using default values
                              """)
                log.critical("""No config file found""")
            else:
                czern = self.config['ZERNIKEDECOMP']
                cpsf = self.config['SCIENCEPSF']

                if self.config['EXTENSIONS']:
                    exts = self.config['EXTENSIONS']
                    try:
                        self.sciext = int(exts['sci_ext'])
                        self.templext = int(exts['templ_ext'])
                        self.subext = int(exts['sub_ext'])
                    except:
                        print('No extensions defined')

            #TODO Don't hardcode this


            if isinstance(inimage, str):
                self.inimage = inimage
            # self.fullinimage = self.inimage+'[{}]'.format(self.sciext)
            incat_file = self.inimage.replace('.fits','_in.cat')
            self.fullincat = os.path.join(self.catpath,incat_file)

            if isinstance(catpath, str):
                self.catpath = catpath
                os.makedirs(self.catpath,exist_ok=True)

            if isinstance(outpath, str):
                self.outpath = outpath
                os.makedirs(self.outpath,exist_ok=True)

            if isinstance(nthreads,float):
                self.iopool = ThreadPool(nthreads)

            subfile = self.inimage
            # self.fullsub = self.inimage+'[{}]'.format(self.subext)
            self.fullsub = subfile
            subcat_file = subfile.replace('.fits','.cat')
            self.fullsubcat = os.path.join(self.catpath,subcat_file)

            reffile = self.inimage
            # self.fullref = self.inimage+'[{}]'.format(self.templ_ext)
            self.fullref = reffile
            ref_cat_file = reffile.replace('.fits','.cat')
            self.fullrefcat = os.path.join(self.catpath,ref_cat_file)

            # TODO TURN ON WHEN CONVOLVED IMAGE OUTPUT FROM HOTPANTS WORKS
            # convfile = self.inimage
            # self.fullconv = self.inimage+'[{}]'.format(self.templ_ext)
            # conv_cat_file = convfile.replace('.fits','.cat')
            # self.fullconvcat = os.path.join(self.catpath,conv_cat_file)

            self.nzmax = int(czern['nz_max'])
            self.disk_dim = int(czern['disk_dim'])
            self.disk_rad = float(czern['disk_rad'])
            self.sub_dim = int(czern['sub_dim'])
            self.dzmax = float(czern['dz_max'])

            global j_vec,zern_list,disk_mask,cov_mat_inv
            j_vec = range(1,self.nzmax+1)
            zern_list,disk_mask = zern_poly.create_fzern_list(j_vec,self.disk_dim,self.disk_rad,0.0,0.0,self.sub_dim)
            cov_mat_inv = zern_poly.inv_cov_mat(zern_list)

            try:
                self.cls_min = float(cpsf['cls_sci_min'])
            except KeyError:
                warnings.warn("No CLS_STAR limits given in param file")

            try:
                self.cls_max = float(cpsf['cls_sci_max'])
            except KeyError:
                warnings.warn("No CLS_STAR limits given in param file")

            try:
                self.flag_req = int(cpsf['flag_sci_req'])
            except KeyError:
                self.flag_req = 0
                warnings.warn("No FLAG_REQ given in param file. Using Flag=0")

            try:
                zp_key = cpsf['zp_keyword']
                if not zp_key:
                    self.zpt = 0
                    self.zpt_err = 0
                    warnings.warn('No ZP Keyword header given in config file. Using ZP=0')
                else:
                    self.zpt = fits.getheader(self.inimage)[zp_key]
                    self.zpt_err = fits.getheader(self.inimage)[zp_key+'*ERR*']
                    if not self.zpt_err:
                        print('Using ZP_ERR=0')
                        self.zpt_err = 0
            except KeyError:
                self.zpt = 0
                self.zpt_err = 0
                warnings.warn('No ZP Keyword header defined at all. Using ZP=0')

    # def check_reference(self,inimage,refimage):
    #     from spherical_geometry import polygon
    #     hdr = fits.getheader(inimage,ext=self.sciext)
    #     img_wcs = wcs.WCS(hdr)
    #
    #     rhdr = fits.getheader(refimage,ext=self.templext)
    #     ref_wcs = wcs.WCS(rhdr)
    #
    #     s1 = polygon.SphericalPolygon.from_wcs(img_wcs)
    #     s2 = polygon.SphericalPolygon.from_wcs(ref_wcs)
    #
    #     return round(s1.overlap(s2),5)

    def filter_catalog(self,cat_file,img_shape,**kwargs):

        try:
            x_key = kwargs['x_key']
            y_key = kwargs['y_key']
        except:
            x_key = 'XWIN_IMAGE'
            y_key = 'YWIN_IMAGE'
        try:
            flux_key = kwargs['flux_key']
            fluxerr_key =kwargs['fluxerr_key']
        except:
            flux_key = 'FLUX_AUTO'
            fluxerr_key ='FLUXERR_AUTO'
        try:
            mag_key = kwargs['mag_key']
            magerr_key = kwargs['mag_key']
        except:
            mag_key = 'MAG_AUTO'
            magerr_key = 'MAGERR_AUTO'

        cat     = fileIO.read_file(cat_file)
        print( str(len(cat))," objects found in image catalog.")
        cat_filter = misc.filter_cat_flag(cat,self.flag_req)
        print( str(len(cat_filter))," objects remain in image catalog after",
        " FLAG={} filter.".format(self.flag_req))
        del(cat)

        cat_filter2 = misc.filter_cat_cls(cat_filter,self.cls_min,self.cls_max)
        print( str(len(cat_filter2))," objects remain in image catalog after",
         " CLASS_STAR filter {} to {}.".format(self.cls_min,self.cls_max))
        del(cat_filter)
        ## TODO: Add SATURATION FLAG
        catflux = cat_filter2.to_pandas()[flux_key]
        catflux_clip = sigma_clip(catflux, sigma=3, maxiters=3)
        self.flux_max = catflux_clip.quantile(0.95)
        self.flux_min = catflux_clip.quantile(0.05)
        del catflux, catflux_clip
        print('FLUX MAX AND MIN',self.flux_max,self.flux_min)
        cat_filter3 = misc.filter_cat_flux(cat_filter2,self.flux_min,self.flux_max)
        print( str(len(cat_filter3))," objects remain in image catalog after flux filter.")
        del cat_filter2
        cat_filter4 = misc.filter_cat_xy(cat_filter3,img_shape,self.sub_dim)
        print( str(len(cat_filter4))," objects remain in image catalog after filtering ",
                "sources {} pixels from edge.".format(self.sub_dim))
        del cat_filter3
        return cat_filter4


    def make_psf(self,cat_file,img_data,**kwargs):
        """
        Calculate the PSF using Zernike polynomical decomposition.

        Parameters:
        :param img_data: the data from the respective FITS HDU, it's defined here
        as the conv file.
        :param cat_filter: the filtered image SExtractor catalog

        """
        print('Making the PSF')
        try:
            x_key = kwargs['x_key']
            y_key = kwargs['y_key']
        except:
            x_key = 'XWIN_IMAGE'
            y_key = 'YWIN_IMAGE'

        cat_filter = self.filter_catalog(cat_file,img_data.shape)
        print(cat_file)
        wf_s_coeff_list = []
        for cat_num, cat_item in enumerate(cat_filter):
            x           = cat_item[x_key]
            y           = cat_item[y_key]
            x_int       = int(round(x))
            y_int       = int(round(y))
            dx          = x-x_int
            dy          = y-y_int

                # psfs.append(impsf.img2wf_lanczos3(img_data,x,y,self.sub_dim,disk_mask))
            # wf_s_img    = impsf.img2wf_lanczos3(img_data,x,y,self.sub_dim,disk_mask)
            sub_img     = impsf.sub_image(img_data,x_int,y_int,self.sub_dim)
            sub_img     = impsf.remove_background(sub_img)
            sub_img_s   = impsf.shift_image(sub_img,dx,dy,5)
            wf_s_img    = impsf.calc_wavefront(sub_img_s,disk_mask)
            if not isinstance(wf_s_img,str):  #???????
                wf_s_coeff  = zern_poly.zern_fit(wf_s_img,zern_list,cov_mat_inv)
                wf_s_rec    = zern_poly.wf_comp_brief(wf_s_coeff,zern_list)
                wf_s_err    = zern_poly.calc_wf_error(wf_s_img,wf_s_rec,disk_mask)
                if wf_s_err < .01:
                    wf_s_coeff_list.append(wf_s_coeff)

            else:
                print(" A wf_s_coeff_list wasn't created")
                print(wf_s_img)
                continue
        print('{} sources found for calculating the model PSF'.format(len(wf_s_coeff_list)))
        if not len(wf_s_coeff_list) > 1:
            #exit
            warnings.warn("""
            Could not write output PSF file. Filtering of sources is too strict.
            """)
        else:
            wf_s_coeff_list     = np.array(wf_s_coeff_list)
            # Calculate mean, median, and stds of coefficients for all sources
            wf_s_coeff_mean     = np.mean(wf_s_coeff_list,axis=0)
            wf_s_coeff_mean_std = np.std(wf_s_coeff_list,axis=0)
            wf_s_coeff_med      = np.median(wf_s_coeff_list,axis=0)
            wf_s_coeff_med_std  = np.zeros(len(wf_s_coeff_med))

            for k in range(len(wf_s_coeff_med)):
                for l in range(len(wf_s_coeff_list)):
                    wf_s_coeff_med_std[k] += (wf_s_coeff_list[l][k]-wf_s_coeff_med[k])**2
                wf_s_coeff_med_std[k]    = np.sqrt(wf_s_coeff_med_std[k])
            n_bins      =   int(np.floor(10*np.log10(len(wf_s_coeff_list))))
            try:
                gfit_res                = misc.gfit_zerndist(wf_s_coeff_list,n_bins)
                wf_s_coeff_gauss        = np.array(gfit_res[0])
                wf_s_coeff_gauss_std    = np.array(gfit_res[1])
                wf_s_coeff_gauss_chi2   = np.array(gfit_res[2])
            except RuntimeError:
                print('RunTimeError encountered when calculating Gaussian Fits')
                wf_s_coeff_gauss        = np.array([-99.0]*self.nzmax)
                wf_s_coeff_gauss_std    = np.array([-99.0]*self.nzmax)
                wf_s_coeff_gauss_chi2   = np.array([-99.0]*self.nzmax)
                pass
            #TODO  Import and run tiling script
            len_wf_s_coeff_list = len(wf_s_coeff_list)
            len_cat_filter = len(cat_filter)
            return len_wf_s_coeff_list, wf_s_coeff_med,wf_s_coeff_med_std
            # fileIO.write_psf_file(self.fullconvpsf,self.fullinterp,self.fullconvcat,
            #     wf_s_coeff_mean,wf_s_coeff_med,wf_s_coeff_gauss,self.flux_min,self.cls_min,
            #     self.nzmax,self.sub_dim, self.disk_dim,self.disk_rad,len_wf_s_coeff_list,
            #     len_cat_filter,wf_s_coeff_gauss_std,wf_s_coeff_gauss_chi2,
            #     wf_s_coeff_mean_std,wf_s_coeff_med_std)
            # return


    def plot_psf_array(self,cat_file,img_data,ext=0,**kwargs):

        try:
            x_key = kwargs['x_key']
            y_key = kwargs['y_key']
        except:
            x_key = 'XWIN_IMAGE'
            y_key = 'YWIN_IMAGE'

        cat_filter = self.filter_catalog(cat_file,img_data.shape)

        sep_val = np.linspace(-.5,.5,4)
        cent_val = np.array([(sep_val[i]+sep_val[i+1])/2.0 for i in range(len(sep_val)-1)])

        tile_cat    = [[[cat_item for cat_item in cat_filter
                if sep_val[i] <= (cat_item[x_key]-int(round(cat_item[x_key]))) < sep_val[i+1]
                and sep_val[j] <= (cat_item[y_key]-int(round(cat_item[y_key]))) < sep_val[j+1]]
                for i in range(len(sep_val)-1)] for j in range(len(sep_val)-1)]

        wf_list_tile = [[[impsf.img2wf(img_data,cat_item[x_key],cat_item[y_key],self.sub_dim,disk_mask)
                    for cat_item in tile_cat[j][i]] for i in range(3)] for j in range(3)]

        psf_wf_tile=np.zeros((3,3,self.disk_dim,self.disk_dim))
        for m in range(0,3):
            for n in range(0,3):
                wf_clip = sigma_clip(wf_list_tile[m][n], sigma=3, maxiters=3)
                wf_weights      = [1/np.var(wf_item) if np.var(wf_item) else 1 for wf_item in wf_clip ]
                try:
                    wf_avg         = np.average(wf_clip,axis=0,weights=wf_weights)*disk_mask
                except ZeroDivisionError:
                    print('ZeroDivisionError encountered when calculating the weighted average')
                    wf_avg         = np.average(wf_clip,axis=0)*disk_mask
                psf_wf_tile[m,n]     = wf_avg/wf_avg.sum()
                print(m,n,len(wf_list_tile[m][n]))

        tile_ref_name = outpath+os.path.splitext(cat_file)[0]+'_ref_tile.png'
        tile_name = outpath+os.path.splitext(cat_file)[0]+'_ref.png'

        psf_inj = np.array(psf_wf_tile,copy=True)

        f, ((ax1, ax2, ax3), (ax4, ax5, ax6), (ax7, ax8, ax9)) = plt.subplots(3,3)

        ax2.set_title(os.path.basename(cat_file))
        ax1.imshow(psf_inj[0][0])
        ax2.imshow(psf_inj[0][1])
        ax3.imshow(psf_inj[0][2])
        ax4.imshow(psf_inj[1][0])
        ax5.imshow(psf_inj[1][1])
        ax6.imshow(psf_inj[1][2])
        ax7.imshow(psf_inj[2][0])
        ax8.imshow(psf_inj[2][1])
        cax = ax9.imshow(psf_inj[2][2])
        f.colorbar(cax)
    #    plt.show()
        plt.savefig(tile_ref_name)
        plt.close()

        f, ax = plt.subplots(1,1)
        ax.set_title(os.path.basename(cat_file))
        cax = ax.imshow(psf_inj[1][1])
        f.colorbar(cax)
        plt.savefig(tile_name)
        plt.close()

        return tile_ref_name, tile_ref

    def make_transient_catalog(self,cat_file,img_data,zern_stats,*kwargs):

        try:
            x_key = kwargs['x_key']
            y_key = kwargs['y_key']
            flux_key = kwargs['flux_key']
            fluxerr_key =kwargs['fluxerr_key']
        except:
            x_key = 'XWIN_IMAGE'
            y_key = 'YWIN_IMAGE'
            flux_key = 'FLUX_AUTO'
            fluxerr_key ='FLUXERR_AUTO'

        cat     = fileIO.read_file(cat_file)
        cat_filter  = misc.filter_cat_xy(cat,img_data.shape,self.sub_dim)
        sys.stdout.write("Done\n")
        sys.stdout.write("Performing Zernike decompositions ...")
        sys.stdout.flush()
        trans_list  = []
        for cand_item in cat_filter:
            if not cand_item[fluxerr_key] == 0:
                x           = cand_item[x_key]
                y           = cand_item[y_key]
                x_int       = int(round(x))
                y_int       = int(round(y))
        #        sub_img     = ppln_misc.sub_image(sub_data,x,y,sub_dim)
        #        sub_img     = ppln_misc.remove_background(sub_img)
        #        sub_img_s   = ppln_misc.shift_image(sub_img,dx,dy,5)
        #        wf_s_img    = ppln_misc.calc_wavefront(sub_img_s,disk_mask)
                wf_s_img = impsf.img2wf_lanczos3(img_data,x,y,self.sub_dim,disk_mask)
                if not isinstance(wf_s_img,str):
                    wf_s_coeff = zern_poly.zern_fit(wf_s_img,zern_list,cov_mat_inv)
                    wf_s_rec = zern_poly.wf_comp_brief(wf_s_coeff,zern_list)
                    wf_s_err = zern_poly.calc_wf_error(wf_s_img,wf_s_rec,disk_mask)
                    wf_s_z_dist = misc.calc_zdist_med_std(wf_s_coeff,zern_stats)

                    flux = cand_item[flux_key]
                    flux_err = cand_item[fluxerr_key]
                    snr = flux/flux_err
                    ra = cand_item['ALPHA_J2000']
                    dec = cand_item['DELTA_J2000']
                    mag_inst = cand_item[mag_key]
                    mag_inst_err = cand_item[magerr_key]
                    mag_app = mag_inst + self.zpt
                    mag_app_err = np.sqrt(mag_inst_err**2+self.zpt_err**2)
                    fwhm_image = cand_item['FWHM_IMAGE']

                    if wf_s_z_dist < self.dzmax:
                        trans_list.append([x,y,snr,wf_s_z_dist,
                                mag_inst,mag_inst_err,ra,dec,flux,flux_err,mag_app,
                                mag_app_err,fwhm_image])
                else:
                    continue
        print('Number of transients recovered',len(trans_list))
        return trans_list

    def make_temp_catalogs(self,inimage,catpath,outpath,imgtype):

        self.sfiles = '/Users/kack0001/astro-tools/goto/ppln_goto/param/goto.sex'
        print('SEXTRACTOR FILE: ',self.sfiles)

        if isinstance(catpath, str):
            self.catpath = catpath
            os.makedirs(self.catpath,exist_ok=True)

        if isinstance(outpath, str):
            self.outpath = outpath
            os.makedirs(self.outpath,exist_ok=True)

        fitsfile = inimage.copy()
        self.fullcat = fitsfile

        if imgtype == 'sci':
            if isinstance(inimage, str):
                self.inimage = inimage
            cat_file = self.inimage.replace('.fits','_in.cat')
            fullsextrcat = os.path.join(self.catpath,cat_file)
            imgext = '[{}]'.format(self.sciext)
            imgextint = int(self.sciext)
            print('Running SEXtractor on the reference image: ',fullsextrcat)
        elif imgtype == 'sub':
            cat_file = fitsfile.replace('.fits','_ref.cat')
            fullsextrcat = os.path.join(self.catpath,cat_file)
            imgext = '[{}]'.format(self.templext)
            imgextint = int(self.templext)
            print('Running SEXtractor on the reference image: ',fullsextrcat)
        elif imgtype == 'templ':
            cat_file = fitsfile.replace('.fits','_sub.cat')
            fullsextrcat = os.path.join(self.catpath,cat_file)
            imgext = '[{}]'.format(self.subext)
            imgextint = int(self.subext)
            print('Running SEXtractor on the subtracted image: ',fullsextrcat)
        else:
            print("Please give imgtype name as 'sci', 'sub', or 'templ'.")

        fileIO.sextractor_script(fitsfile+imgext, fullsextrcat, self.sfiles)
        cat_data,cat_hdr = fits.getdata(fitsfile,ext=imgextint,header=True)

        return fullsextrcat, cat_data

        #This is if you already have your template image definied prior to running the script.
        #FIX THIS!!!! DO NOT USE THE NAME FOR THE REFERENCE IMAGE
        # Make Model PSF


    def make_temp_sub_catalogs(self,inimage,catpath,outpath):

        sub_ext = '[{}]'.format(self.subext)
        print('Running SEXtractor on the subtracted image: ',self.fullsub)
        fileIO.sextractor_script(self.fullsub+sub_ext, self.fullsubcat, self.sfiles)

        print("Making the transient catalog from the subtracted image")
        sub_data = fits.getdata(self.fullsub,ext=self.subext,header=False)
        # sub_wcs     = wcs.WCS(sub_hdr)
        zern_stats = {'MED':wf_s_coeff_med, 'MED_STD': wf_s_coeff_med_std}
        self.make_transient_catalog(self.fullsubcat,sub_data,zern_stats)

            # trans_list.append([x,y,snr,wf_s_z_dist,
            #         mag_inst,mag_inst_err,ra,dec,flux,flux_err,mag_app,
            #         mag_app_err,seeing])

        if self.del_tmp:
            os.remove(self.fullincat)
            os.remove(self.fullrefcat)
            os.remove(self.fullsubcat)



        # # initialise log
        # if log is None:
        #     log = logging.getLogger() #create logger
        #     log.setLevel(logging.INFO) #set level of logger
        #     formatter = logging.Formatter("%(asctime)s %(levelname)s %(message)s") #set format of logger
        #     logging.Formatter.converter = time.gmtime #convert time in logger to UCT
        #     if new_fits is not None:
        #         filehandler = logging.FileHandler(new_fits.replace('.fits','.log'), 'w+') #create log file
        #     elif ref_fits is not None:
        #         filehandler = logging.FileHandler(ref_fits.replace('.fits','.log'), 'w+') #create log file
        #     else:
        #         filehandler = logging.FileHandler('log', 'w+') #create log file
        #     filehandler.setFormatter(formatter) #add format to log file
        #     log.addHandler(filehandler) #link log file to logger
        #     if get_par(C.verbose,tel):
        #         streamhandler = logging.StreamHandler() #create print to screen logging
        #         streamhandler.setFormatter(formatter) #add format to screen logging
        #         log.addHandler(streamhandler) #link logger to screen logging

if __name__ == "__main__":

    ap = argparse.ArgumentParser()
    ap.add_argument('-i', '--inimage', nargs='?', help="List of images", type=str)
    ap.add_argument('-r', '--refimage', nargs='?', help="List of reference images", type=str)
    ap.add_argument('-c', '--config', default=None, help="Configuration parameter file")
    ap.add_argument('-cp', '--catpath', type=str,help='Path to output SEXtractor catalog')
    ap.add_argument('-nt', '--nthreads', default=1, type=int, help='Number of threads for multithreading (IO)')
    # ap.add_argument('-ncpu','--ncpus', default=1, type=int, help='Number of CPUs for multiprocessing (Computational)')
    ap.add_argument('-v', '--verbose', default=False, action='store_true', help='Increase verbosity')
    ap.add_argument('--del_tmp', default=False, action='store_true', help='Delete Temporary files')
    ap.add_argument('--make_trans',default=False, help='Make transient catalog from Subtracted Image using Model PSF Params from Reference Image (default=False)')
    ap.add_argument('--plot_psf',default=False, help='Return the median sigma-clipped PSF for plotting')
    args = vars(ap.parse_args())
    print('CALLING THE SCRIPT')
    r = Shapelet(**args)
    r.main()
