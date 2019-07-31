from . import psf
from . import base
import argparse
import astropy.io.fits as fits
from astropy.io import ascii
from astropy import wcs
from astropy.table import Table
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

import warnings

zern_poly = psf.ZernikePoly()
fileIO = base.FileIO()
misc = psf.Catalog()
impsf = psf.PSF()
zernobj = base.ZernObj()
trans = base.TransObj()

class Shapelet:

    def __init__(self,inimage,refimage,config,inpath,psfpath=None,catpath=None,
                transpath=None,outpath=None,nthreads=1,verbose=False,
                del_tmp=False,wcs_flag=True):

        if isinstance(inimage, str):
            self.inimage = inimage

        if isinstance(refimage, str):
            self.refimage = refimage

        if isinstance(inpath, str):
            self.inpath = inpath
            self.refpath = inpath
            #TODO fix refpath

        if isinstance(psfpath, str):
            self.psfpath = psfpath

        if isinstance(catpath, str):
            self.catpath = catpath

        if isinstance(transpath, str):
            self.transpath = transpath

        if isinstance(outpath, str):
            self.outpath = outpath

        if isinstance(nthreads,float):
            self.iopool = ThreadPool(nthreads)

        # if isinstance(ncpus,float):
        #     self.cpupool = Pool(ncpus)

        self.verbose = verbose
        self.del_tmp = del_tmp
        self.wcs_flag = wcs_flag



        subfile = self.inimage.replace('.fits','_sub.fits')
        self.fullsub = os.path.join(self.outpath,subfile)

        subcat_file = subfile.replace('.fits','.cat')
        self.fullsubcat = os.path.join(self.catpath,subcat_file)

        incat_file = self.inimage.replace('.fits','.cat')
        self.fullincat = os.path.join(self.catpath,incat_file)

        refcat_file = self.refimage.replace('.fits','.cat')
        self.fullrefcat = os.path.join(self.catpath,refcat_file)

        interp_file = self.refimage.replace('.fits','_interp_tmpl.fits')
        self.fullinterp = os.path.join(self.refpath,interp_file)

        interp_cat_file = interp_file.replace('.fits','.cat')
        self.fullinterpcat = os.path.join(self.catpath,interp_cat_file)

        conv_file = self.refimage.replace('.fits','_conv.fits')
        self.fullconv = os.path.join(self.refpath,conv_file)

        conv_cat_file = conv_file.replace('.fits','.cat')
        self.fullconvcat = os.path.join(self.catpath,conv_cat_file)

        conv_psf_file = conv_file.replace('.fits','.psf')
        self.fullconvpsf = os.path.join(self.psfpath,conv_psf_file)

        trans_file = self.fullinterpcat.replace('.cat','.trans')
        self.fulltrans = os.path.join(self.transpath,trans_file)

        if isinstance(config,str):
            try:
                self.config = ConfigObj(config,file_error=True)
            except OSError:
                print("something's wrong")
                # raise SystemExit

            if not len(self.config.sections):
                warnings.warn(""""
                No config file found. Using default values

                              Type default values here
                              """)
                log.critical("""No config file found""")
            else:
                paths = self.config['PATHS']
                czern = self.config['ZERNIKEDECOMP']
                cpsf = self.config['SCIENCEPSF']
                self.chp = self.config['HOTPANTSARGS']

            self.projpath = paths['project_dir']
            if self.refpath is None:
                self.refpath = paths['ref_path']
            if self.inpath is None:
                self.inpath = paths['img_path']

            self.parampath = self.projpath+'/zernike/param/'
            #TODO Don't hardcode this

            self.incatpath = self.catpath
            self.refcatpath = self.catpath
            self.tmplpath = paths['templ_path']
            # if self.subpath is None:
            #     self.subpath = paths['sub_path']
            if self.transpath is None:
                self.transpath = paths['trans_path']

            self.convpath = self.refpath
            self.conv_cat_path = self.refcatpath
            self.conv_psf_path = self.psfpath

            # conv_cat_file = ref_cat.replace('.cat','_conv.cat')
            # conv_psf_path = psf_path
            # conv_psf_file = psf_file

            self.fullin = os.path.join(self.inpath,self.inimage)
            self.fullref = os.path.join(self.refpath,self.refimage)

            self.nzmax = int(czern['nz_max'])
            self.disk_dim = int(czern['disk_dim'])
            self.disk_rad = float(czern['disk_rad'])
            self.sub_dim = int(czern['sub_dim'])
            self.dzmax = float(czern['dz_max'])

            self.overlap_frac = float(self.chp['align_overlap_frac'])


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
                mag_min = float(cpsf['mag_sci_min'])
                mag_max = float(cpsf['mag_sci_max'])
                self.flux_min = 10**(-mag_min/2.5)
                self.flux_max = 10**(-mag_max/2.5)
            except KeyError:
                try:
                    self.flux_min = float(cpsf['flux_sci_min'])
                    self.flux_max = float(cpsf['flux_sci_max'])
                except KeyError:
                    warnings.warn("No Mag or Flux limits given in param file")
                    # log.critical("""No Mag or Flux given in param file""")

            try:
                zp_key = cpsf['zp_keyword']
                if not zp_key:
                    self.zpt = 0
                    self.zpt_err = 0
                    warnings.warn('No ZP Keyword header given. Using ZP=0')
                else:
                    self.zpt = fits.getheader(self.inimage)[zp_key]
                    self.zpt_err = fits.getheader(self.inimage)[zp_key+'*ERR*']
                    if not self.zpt_err:
                        warnings.warn('Using ZP_ERR=0')
                        self.zpt_err = 0
            except KeyError:
                self.zpt = 0
                self.zpt_err = 0
                warnings.warn('No ZP Keyword header given. Using ZP=0')





    def check_reference(self,inimage,refimage):
        from spherical_geometry import polygon
        hdr = fits.getheader(inimage,ext=0)
        img_wcs = wcs.WCS(hdr)

        rhdr = fits.getheader(refimage,ext=0)
        ref_wcs = wcs.WCS(rhdr)

        s1 = polygon.SphericalPolygon.from_wcs(img_wcs)
        s2 = polygon.SphericalPolygon.from_wcs(ref_wcs)

        return s1.overlap(s2)


    def make_psf(self,cat_file,img_data,**kwargs):
        """
        Calculate the PSF using Zernike polynomical decomposition.

        Parameters:
        :param img_data: the data from the respective FITS HDU
        :param cat_filter: the filtered image SExtractor catalog

        """

        try:
            x_key = kwargs['x_key']
            y_key = kwargs['y_key']
        except:
            x_key = 'XWIN_IMAGE'
            y_key = 'YWIN_IMAGE'

        cat     = fileIO.read_file(cat_file)
        cat_filter  = misc.filter_cat_xy(cat,img_data.shape,self.sub_dim)

        wf_s_coeff_list = []

        for cat_num, cat_item in enumerate(cat_filter):
            x           = cat_item[x_key]
            y           = cat_item[y_key]
            wf_s_img    = impsf.img2wf_lanczos3(img_data,x,y,self.sub_dim,disk_mask)
            if not isinstance(wf_s_img,str):

                wf_s_coeff  = zern_poly.zern_fit(wf_s_img,zern_list,cov_mat_inv)
                wf_s_rec    = zern_poly.wf_comp_brief(wf_s_coeff,zern_list)
                wf_s_err    = zern_poly.calc_wf_error(wf_s_img,wf_s_rec,disk_mask)
                if wf_s_err < .01:
                    wf_s_coeff_list.append(wf_s_coeff)
            else:
                continue
        if not len(wf_s_coeff_list) > 1:
            #exit
            warnings.warn("""
            Could not write output PSF file. Filtering of sources is too strict.
            """)
        else:
            wf_s_coeff_list     = np.array(wf_s_coeff_list)
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
                    wf_s_coeff_gauss        = np.array([-99.0]*self.nzmax)
                    wf_s_coeff_gauss_std    = np.array([-99.0]*self.nzmax)
                    wf_s_coeff_gauss_chi2   = np.array([-99.0]*self.nzmax)
                    pass

            len_wf_s_coeff_list = len(wf_s_coeff_list)
            len_cat_filter = len(cat_filter)
            fileIO.write_psf_file(self.fullconvpsf,self.fullinterp,self.fullconvcat,
                wf_s_coeff_mean,wf_s_coeff_med,wf_s_coeff_gauss,self.flux_min,self.cls_min,
                self.nzmax,self.sub_dim, self.disk_dim,self.disk_rad,len_wf_s_coeff_list,
                len_cat_filter,wf_s_coeff_gauss_std,wf_s_coeff_gauss_chi2,
                wf_s_coeff_mean_std,wf_s_coeff_med_std)


    def make_transient_catalog(self,cat_file,img_data,zern_stats,*kwargs):

        try:
            x_key = kwargs['x_key']
            y_key = kwargs['y_key']
            flux_key = kwargs['flux_key']
            fluxerr_key =kwargs['fluxerr_key']
            mag_key = kwargs[mag_key]
            magerr_key = kwargs[mag_key]
        except:
            x_key = 'XWIN_IMAGE'
            y_key = 'YWIN_IMAGE'
            flux_key = 'FLUX_AUTO'
            fluxerr_key ='FLUXERR_AUTO'
            mag_key = 'MAG_AUTO'
            magerr_key = 'MAGERR_AUTO'

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
                    seeing = cand_item['FWHM_IMAGE']

                    if wf_s_z_dist < self.dzmax:
                        trans_list.append([x,y,snr,wf_s_z_dist,
                                mag_inst,mag_inst_err,ra,dec,flux,flux_err,mag_app,
                                mag_app_err,seeing])
                else:
                    continue
        print(len(trans_list))
        if len(trans_list):
            full_reg     = self.fulltrans.replace('.trans','.reg')

            # s_trans_reg10   = cat_file.replace('.cat','dz10.reg')
            len_cat = len(cat)
            fileIO.write_trans_file(self.fulltrans,trans_list,len_cat,self.dzmax)
            fileIO.write_reg(full_reg,trans_list,filter_trans=0)
            # fileIO.write_reg(fullreg,trans_list,filter_trans=1,filter_dz=15)
        else:
            warnings.warn("No transients found in image!")



    def main(self):

        #TODO -- add conversion from magnitude if there is a zp in the header
        global zern_list,disk_mask,cov_mat_inv

        j_vec = range(1,self.nzmax+1)
        zern_list,disk_mask = zern_poly.create_fzern_list(j_vec,self.disk_dim,self.disk_rad,0.0,0.0,21)
        cov_mat_inv = zern_poly.inv_cov_mat(zern_list)

        self.sfiles = self.parampath + 'sex.config'

        if not os.path.isfile(self.fullsub):
            warnings.warn("Running Alignment & HOTPANTS") #logger instead

            overlap = self.check_reference(self.fullin,self.fullref)
            print(overlap)
            if overlap < 0.05:
                warnings.warn("The two images have no overlap in WCS coordinates")
            elif overlap < self.overlap_frac:
                warnings.warn("The two images have less than {}% in WCS coordinates".format(self.overlap_frac*100))
            warnings.warn("The images have {}% overlap".format(overlap*100))


            #TODO: don't hardcode this in

            if self.wcs_flag:
                fileIO.remap_wcs(self.fullin, self.fullref, self.fullinterp)
            else:
                if not os.path.isfile(self.fullincat):
                    fileIO.sextractor_script(self.fullin, self.fullincat, self.sfiles)
                img_cat = fileIO.read_file(self.fullincat)

                if not os.path.isfile(self.fullrefcat):
                    fileIO.sextractor_script(self.fullref, self.fullrefcat, self.sfiles)
                ref_cat = fileIO.read_file(self.fullrefcat)
                fileIO.remap_spalipy(self.fullincat, self.fullrefcat,
                                    self.fullin, self.fullinterp)
                #TODO X_IMAGE, Y_IMAGE, FLUX_BEST, FWHM_IMAGE, FLAGS

                if self.del_tmp:
                    os.remove(self.fullincat)
                    os.remove(self.fullrefcat)


            if not os.path.isfile(self.fullinterpcat):
                fileIO.sextractor_script(self.fullinterp,self.fullinterpcat,self.sfiles)
            ref_aligned_cat = fileIO.read_file(self.fullinterpcat)
            ref_data,ref_hdr = fits.getdata(self.fullinterp,ext=0,header=True)

            cat_filter1 = misc.filter_cat_cls(ref_aligned_cat,self.cls_min,self.cls_max)
            cat_filter2 = misc.filter_cat_flag(cat_filter1,self.flag_req)
            cat_filter3 = misc.filter_cat_flux(cat_filter2,self.flux_min,self.flux_max)
            cat_filter4 = misc.filter_cat_xy(cat_filter3,ref_data.shape,self.sub_dim)
            cat_filter  = cat_filter4.copy()

            del cat_filter1, cat_filter2, cat_filter3, cat_filter4

            # if not ref_hdr[zpt_key]:
            #     warnings.warn("""Please provide the ZP KEYWORD in the header.
            #                     Else provide a calibration catalog""")

            hp_args = self.chp
            fileIO.run_hotpants(self.fullin,self.fullinterp,self.fullsub,self.fullconv,hp_args)
        fileIO.sextractor_script(self.fullconv,self.fullconvcat,self.sfiles)
        fileIO.sextractor_script(self.fullsub,self.fullsubcat,self.sfiles)

        if not os.path.isfile(self.fullconvpsf):
            conv_data = fits.getdata(self.fullinterp,ext=0)
            self.make_psf(self.fullconvcat,conv_data)
        zern_stats = fileIO.read_zern_file(self.fullconvpsf)

        sub_data,sub_hdr = fits.getdata(self.fullsub,ext=0,header=True)
        sub_wcs     = wcs.WCS(sub_hdr)

        self.make_transient_catalog(self.fullsubcat,sub_data,zern_stats)



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





    # run(args.inimage, args.refimage, args.config, args.out_path, args.nthreads)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument('-i', '--inimage', nargs='?', help="List of images", type=str)
    ap.add_argument('-r', '--refimage', nargs='?', help="List of reference images", type=str)
    ap.add_argument('-c', '--config', default=None, help="Configuration parameter file")
    ap.add_argument('-ip', '--inpath', type=str, help='Path to image directory')
    ap.add_argument('-op','--outpath', type=str,help='Output path')
    ap.add_argument('-nt','--nthreads', default=1, type=int, help='Number of threads for multithreading (IO)')
    # ap.add_argument('-ncpu','--ncpus', default=1, type=int, help='Number of CPUs for multiprocessing (Computational)')
    ap.add_argument('-v', '--verbose', default=False, action='store_true', help='Increase verbosity')
    ap.add_argument('--del_tmp', default=False, action='store_true', help='Delete Temporary files')
    ap.add_argument('--wcs_flag', default=True, action='store_true', help='Image Registration with WCS (default=True)')
    # ap.add_argument('-fs', '--flag_shift', nargs='?', default=1, help='Flag for shifts', dest='flag_shift', type=int)
    # ap.add_argument('-fbkg', '--flag_bkg', default=0,
    #                 help='Flag for adding background to image, helpful for subtracting the same image', dest='flag_bkg', type=int)
    # ap.add_argument('-fshot', '--flag_shot', default=0, help='Flag for adding shot noise', dest='flag_shot', type=int)
    # ap.add_argument('-fblur', '--flag_blur', default=0, help='Flag for blurring image with Gaussian', dest='flag_blur', type=int)
    # ap.add_argument('-bfwhm', '--blur_fwhm', default=0.0, help='Flag for the FWHM of the Gaussian blurring (0.05--0.5)', dest='blur_fwhm', type=float)
    # ap.add_argument('-cmin', '--cls_min', default=0.8,  help='Minimum SExtractor CLASS_STAR', dest='cls_min', type=float)
    # ap.add_argument('-cmax', '--cls_max', default=1.0,  help='Maximum SExtractor CLASS_STAR', dest='cls_max', type=float)
    # ap.add_argument('-flag', '--flag_req', default=0,  help='Flag Requirement SExtractor FLAGS', dest='flag_req', type=int)

    args = vars(ap.parse_args())
    print('CALLING THE SCRIPT')
    r = Shapelet(**args)
    r.main()
