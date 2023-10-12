
# -*- coding: utf-8 -*-

import numpy as np
from astropy.table import Table
from scipy.special import factorial as fac

from . import base
#import base
__all__ = ["zernike"]


class ZernikePoly():

    def noll2mn(self,j):
        """
        converts the Noll index j to the index pair (m,n10)
        """

        n = 0
        k = j-1
        while (k > n):
            n += 1
            k -= n
        if n%4==0 or (n-1)%4==0:
            if j%2:
                m=-k
            else:
                m=k+1
        else:
            if j%2:
                m=-k-1
            else:
                m=k
        return m, n

    def pre_fac(self,k,m,n):
        """
        returns the prefactor for the relative weighting of the rho
        terms in R_mn
        """
        return -2*(k%2-0.5)*fac(n-k)/(fac(k)*fac((n+m)/2-k)*fac((n-m)/2-k))


    def R_mn(self,m,n,rho):
        """
        returns the radial part of Z_mn on the grid rho
        """
        return sum(self.pre_fac(k,m,n)*rho**(n-2*k) for k in range(int(((n-m)/2)+1)))

    def Z_mn(self,m,n,rho,phi):
        """
        returns the full Z_mn on the coordinate grid rho,phi
        """
        if (m >= 0):
            return self.R_mn(m,n,rho)*np.cos(m*phi)
        elif (m < 0):
            return self.R_mn(-m,n,rho)*np.sin(-m*phi)

    def avg_subarr(self,arr,n_ds):
        h,w         = arr.shape
        h           = (h//n_ds)*n_ds
        w           = (w//n_ds)*n_ds
        arr         = arr[:h,:w]
        return np.einsum('ijkl->ik',arr.reshape(h//n_ds,n_ds,-1,n_ds))/n_ds**2

    def create_coord_grid(self,dim,r_disc,dx,dy):
        """
        returns the masked coordinate grid rho,phi
        """

        grid        = (np.indices((dim,dim),dtype=np.float)-((dim-1)/2))
        grid[0]     = (grid[0]+dy)/r_disc
        grid[1]     = (grid[1]+dx)/r_disc
        grid_rho    = (grid**2).sum(0)**0.5
        grid_phi    = np.arctan2(grid[0],grid[1])
        grid_mask   = (grid_rho <= 1)*1.0
        return grid_rho,grid_phi,grid_mask

    def create_coord_fgrid(self,dim,r,dx,dy,n_us):
        """
        returns the masked and upsampled coordinate grid rho,phi
        """

        fdim        = dim*n_us
        fgrid       = (np.indices((fdim,fdim),dtype=np.float)-((fdim-1)/2))/n_us
        fgrid[0]    = (fgrid[0]+dy)/r
        fgrid[1]    = (fgrid[1]+dx)/r
        fgrid_rho   = (fgrid**2).sum(0)**0.5
        fgrid_phi   = np.arctan2(fgrid[0],fgrid[1])
        fgrid_mask  = (fgrid_rho <= 1)*1.0
        fgrid_mask_ds   = self.avg_subarr(fgrid_mask,n_us)
        return fgrid_rho,fgrid_phi,fgrid_mask,fgrid_mask_ds

    def create_zern_list(self,z_ord,dim,r_disc,dx,dy):
        """
        returns the Zernike polynomials over the coordinate grid
        """

        grid_rho,grid_phi,grid_mask = self.create_coord_grid(dim,r_disc,dx,dy)
        mn_list = [self.noll2mn(j) for j in z_ord]
        return [self.Z_mn(mn[0],mn[1],grid_rho,grid_phi)*grid_mask for mn in mn_list],grid_mask

    def create_fzern_list(self,z_ord,dim,r_disc,dx,dy,n_us):
        """
        returns the Zernike polynomials over the upsampled coordinate grid
        """

        fgrid_rho,fgrid_phi,fgrid_mask,fgrid_mask_ds = self.create_coord_fgrid(dim,r_disc,dx,dy,n_us)
        mn_list = [self.noll2mn(j) for j in z_ord]
        return [self.avg_subarr(self.Z_mn(mn[0],mn[1],fgrid_rho,fgrid_phi)*fgrid_mask,n_us) for mn in mn_list],fgrid_mask_ds

    def inv_cov_mat(self,zern_list):
        """
        computes the inverse covariance matrix
        """
        cov_mat = np.array([[np.sum(zerni*zernj) for zerni in zern_list] for zernj in zern_list])
        return np.linalg.pinv(cov_mat)

    def zern_fit(self,wf_img,zern_list,cov_mat_inv):
        """
        Calculates the inner product of each Zernike mode with the test surface
        and given the inner product vector of the test wavefront with Zernike
        basis, calculates the Zernike polynomial coefficients
        """
        wf_zern_inprod  = np.array([np.sum(wf_img*zerni) for zerni in zern_list])
        wf_rec_pow      = np.dot(cov_mat_inv,wf_zern_inprod)
        return wf_rec_pow

    def wf_comp_brief(self,z_vec,zern_list):
        wf_temp = sum(val*zern_list[i] for (i, val) in enumerate(z_vec))
        # print('wf_comp_brief',wf_temp.sum())
        # if wf_temp.sum():
        return wf_temp/wf_temp.sum()
        # else:
        #     return wf_temp

    def wf_comp(self,z_vec,dim):
        zern_list,mask   = self.create_zern_list(len(z_vec),dim)
        return self.wf_comp_brief(z_vec,zern_list)

    def calc_wf_error(self,wf_img,wf_rec,mask):
        n_pix   = np.sum(mask)
        return np.sqrt(np.sum(mask*(wf_img-wf_rec)**2)/n_pix)

    def calc_wf_var(self,wf_img,wf_rec,wf_grid):
        import numpy as np

        n_pix   = np.sum(wf_grid)
        return np.sum((wf_img - wf_rec)**2)/n_pix

    def calc_chi2(self,subim,wf,mask):
        return sum(((subim[i,j]-wf[i,j])**2/wf[i,j]) if mask[i,j] else 0 for i in range(len(subim[0])) for j in range(len(subim)))

class PSF():
    """
    Image and subimage manipulation
    """

    def sub_image(self,img,x_int,y_int,sub_dim):
        """
        Defines the stamp cutout around a source

        :param img:
            The FITS image data
        :param x_int:
            The integer x-coordinate on the image pixel grid
        :param y_int:
            The integer y-coordinate on the image pixel grid
        :param sub_dim:
            The total size of the cutout grid in pixels
        """
        return np.array(img[int(y_int-(sub_dim-1)/2):int(y_int+(sub_dim+1)/2),int(x_int-(sub_dim-1)/2):int(x_int+(sub_dim+1)/2)], copy=True)

    def remove_background(self,img):
        """
        Calulates the median as a proxy for the background and subtracts from image.

        :param img:
            The FITS image data
        """

        return img-np.median(img)

    def shift_image(self,img,dx,dy,interp_order):
        """
        Shifts the image by a fractional pixel amount

        :param img:
            The FITS image data
        :param dx:
            The fractional x-coordinate pixel shift
        :param dy:
            The fractional y-coordinate pixel shift
        :interp_order:
            The interpolation order to resample the shifted PSF
        """
        from scipy import ndimage

        return ndimage.interpolation.shift(img, [-dy,-dx], order=interp_order)

    def shift_image_lanczos3(self,img,dx,dy):
        """
        Shifts the image by a fractional pixel amount using Lanczos 3rd
        order resampling

        :param img:
            The FITS image data
        :param dx:
            The fractional x-coordinate pixel shift
        :param dy:
            The fractional y-coordinate pixel shift
        """
        #needs some fine tuning. produces strange psfs
        temp_dim    = img.shape[0]+6
        img_temp    = np.zeros((temp_dim,temp_dim))
        img_temp[3:img.shape[0]+3,3:img.shape[0]+3]=img
        l_krnl      = np.array([[self.lanczos3_2d(dx-k,dy-l) for k in range(-2,4)] for l in range(-2,4)])
        l_krnl      = l_krnl/np.sum(l_krnl)
        img_new     = np.array([[np.sum(img_temp[3+j-2:3+j+4,3+i-2:3+i+4]*l_krnl) for i in range(img.shape[1])] for j in range(img.shape[0])])
        return img_new

    def calc_wavefront(self,sub_img,disk_mask):
        """
        Returns the normalized and masked image on the unit circle

        :param sub_img:
            The FITS sub-image data (from sub_image)
        :param disk_mask:
            The size of the disk_mask in pixel units
        """

        sub_dim     = sub_img.shape[0]
        disk_dim    = disk_mask.shape[0]
        # min_dim = int((sub_dim-1)/2-(disk_dim-1)/2)
        # max_dim = int((sub_dim-1)/2+(disk_dim+1)/2)
        # disk_img    = sub_img[min_dim:max_dim,min_dim:max_dim]*disk_mask
        disk_img = sub_img[int((sub_dim-1)/2-(disk_dim-1)/2):int((sub_dim-1)/2+(disk_dim+1)/2),int((sub_dim-1)/2-(disk_dim-1)/2):int((sub_dim-1)/2+(disk_dim+1)/2)]*disk_mask
        # print('disk_img',disk_img.sum())
        # if disk_img.sum():
        return disk_img/disk_img.sum()
        # else:
        #     return disk_img

    def img2wf(self,img,x,y,sub_dim,disk_mask):

        """
        Returns the normalized and masked image on the unit circle

        :param img:
            The FITS image data
        :param x:
            The pixel-grid x coordinate from SExtractor catalog
        :param y:
            The pixel-grid y coordinate from SExtractor catalog
        :param sub_dim:
            The size of the cutout stamp around source in pixel units
        :param disk_mask:
            The size of the disk_mask in pixel units
        """

        x_int       = int(round(x))
        y_int       = int(round(y))
        # print('SUB_DIM',sub_dim,x_int,y_int,y)
        sub_img     = self.remove_background(self.sub_image(img,x_int,y_int,sub_dim))
        return self.calc_wavefront(sub_img,disk_mask)

    def img2wf_shift(self,img,x_int,y_int,dx,dy,sub_dim,disk_mask):

        """
        Returns the shifted, normalized and masked image on the unit circle

        :param img:
            The FITS image data
        :param x_int:
            The pixel-grid integer x coordinate from SExtractor catalog
        :param y_int:
            The pixel-grid y integer coordinate from SExtractor catalog
        :param dx:
            The fractional x-coordinate pixel shift
        :param dy:
            The fractional y-coordinate pixel shift
        :param sub_dim:
            The size of the cutout stamp around source in pixel units
        :param disk_mask:
            The size of the disk_mask in pixel units
        """

        sub_img         = self.remove_background(self.sub_image(img,x_int,y_int,sub_dim))
        sub_img_shift   = self.shift_image(sub_img,dx,dy,5)
        return self.calc_wavefront(sub_img_shift,disk_mask)

    def img2wf_lanczos3(self,img,x,y,sub_dim,disk_mask,**kwargs):
        """
        Returns the  normalized and masked image on the unit circle using
        Lanczos 3rd order resampling

        :param img:
            The FITS image data
        :param x:
            The pixel-grid x coordinate from SExtractor catalog
        :param y:
            The y integer coordinate from SExtractor catalog
        :param sub_dim:
            The size of the cutout stamp around source in pixel units
        :param disk_mask:
            The size of the disk_mask in pixel units
        """

        x_int = int(x)
        y_int = int(y)
        dx          = x-x_int
        dy          = x-y_int
        sub_img     = self.remove_background(self.sub_image(img,x_int,y_int,sub_dim))
        disk_dim    = disk_mask.shape
        disk_img    = np.zeros(disk_dim)
        pix_off     = (sub_dim-disk_dim[0])/2

        l_krnl      = np.array([[self.lanczos3_2d(dx-k,dy-l) for k in range(-2,4)] for l in range(-2,4)])

        l_krnl      = l_krnl/np.sum(l_krnl)
        try:
            disk_img    = np.array([[np.sum(sub_img[int(pix_off+j-2):int(pix_off+j+4), \
                                    int(pix_off+i-2):int(pix_off+i+4)]*l_krnl) \
                                    for i in range(disk_mask.shape[1])] \
                                    for j in range(disk_mask.shape[0])])
            disk_img    = disk_img*disk_mask
            if np.sum(disk_img) != 0.0:
                return disk_img/np.sum(disk_img)
            else:
                return 'ERROR'
        except ValueError:
            return 'ERROR'
    #    for i in xrange(disk_dim[1]):
    #        for j in xrange(disk_dim[0]):
    #            for k in xrange(-2,4):
    #                for l in xrange(-2,4):
    #                    disk_img[j,i] += sub_img[pix_off+j+l,pix_off+i+k]*lanczos3_2d(dx-k,dy-l)


    def lanczos3_2d(self,x,y):
        import numpy as np

        return np.sinc(x)*np.sinc(x/3)*np.sinc(y)*np.sinc(y/3)

    def calc_wf_err(self,wf_img,wf_rec,disk_mask):
        import numpy as np

        n_pix   = np.sum(disk_mask)
        return np.sqrt(np.sum(disk_mask*(wf_img - wf_rec)**2)/n_pix)

class Catalog():

    psf = PSF()
    zernpoly = ZernikePoly()

    def filter_cat_xy(self,cat_list,img_shape,sub_dim,**kwargs):
        # print(img_shape)
        # print(img_shape[0])

        x_min       = (sub_dim+1)/2
        x_max       = img_shape[1]-x_min
        y_min       = (sub_dim+1)/2
        y_max       = img_shape[0]-y_min
        print('X AND Y: ',x_min,x_max,y_min,y_max)

        try:
            x_key = kwargs['x_key']
            y_key = kwargs['y_key']
        except:
            x_key = 'XWIN_IMAGE'
            y_key = 'YWIN_IMAGE'
        # return(cat_list[[x_min < cat_list[x_key]] and [cat_list[x_key] < x_max] \
        #         and [y_min < cat_list[y_key]] and [cat_list[y_key] < y_max]])

        return cat_list[(cat_list[x_key] < x_max) & (cat_list[x_key] > x_min) & \
                        (cat_list[y_key] < y_max) & (cat_list[y_key] > y_min)]



    def filter_cat_flux(self,cat_list,flux_min,flux_max,**kwargs):
        try:
            flux_key = kwargs['flux_key']
        except:
            flux_key = 'FLUX_AUTO'
        print('FLUX_MIN = {}, FLUX_MAX = {}'.format(flux_min,flux_max))

        return cat_list[(cat_list[flux_key] <= flux_max) & (cat_list[flux_key] > flux_min)]

    def filter_cat_mag(self,cat_list,mag_min,mag_max, **kwargs):
        try:
            mag_key = kwargs['mag_key']
        except:
            mag_key = 'MAG_AUTO'
        return cat_list[(cat_list[mag_key] <= mag_max) & (cat_list[mag_key] > mag_min)]
        # return(cat_list[[mag_min < cat_list[mag_key]] and [cat_list[mag_key] <= mag_max]])

    def filter_cat_cls(self,cat_list,cls_min,cls_max, **kwargs):
        return cat_list[(cat_list['CLASS_STAR'] <= cls_max) & (cat_list['CLASS_STAR'] > cls_min)]
        # return(cat_list[[cat_list['CLASS_STAR'] <= cls_max] and [cat_list['CLASS_STAR'] > cls_min ]])

    def filter_cat_flag(self,cat_list,flag_val):
        return(cat_list[cat_list['FLAGS'] == flag_val])

    def filter_cat_radec(self,cat_list,ra,dec):
        return(cat_list[[cat_list['ALPHA_J2000'] != ra] and [cat_list['DELTA_J2000'] != dec]])

    ### Zernike calculations

    def cat2zern(self,cat_list,img_data,surr_dim,disc_dim,chi2_thr,zern_list,grid_mask,cov_mat_inv):

        zc_list = []
        for cat_item in cat_list:
            wf_img      = psf.img2wf(img_data,cat_item.x,cat_item.y,surr_dim,disc_dim)
            zc_item     = zernpoly.zern_fit(wf_img,zern_list,cov_mat_inv)
            wf_rec      = zernpoly.wf_comp_brief(zc_item,zern_list,grid_mask)
            if zernpoly.calc_chi2(wf_img,wf_rec,grid_mask)<chi2_thr:
                zc_list.append(zc_item)
        return(zc_list)

    def cat2zern_wf_shift(self,cat_list,img_data,surr_dim,disc_dim,chi2_thr,zern_list,grid_mask,cov_mat_inv):


        zc_list     =   []
        for i, cat_item in enumerate(cat_list):
            dx          = cat_item.x - int(round(cat_item.x))
            dy          = cat_item.y - int(round(cat_item.y))
            wf_img      = psf.img2wf_shift(img_data,cat_item.x,cat_item.y,dx,dy,surr_dim,disc_dim)
            zc_item     = zernpoly.zern_fit(wf_img,zern_list,cov_mat_inv)
            wf_rec      = zernpoly.wf_comp_brief(zc_item,zern_list,grid_mask)
            if zernpoly.calc_chi2(wf_img,wf_rec,grid_mask)<chi2_thr:
                zc_list.append(zc_item)
        return(zc_list)

    ### Gauss fitting

    def calc_zdist_gauss(self,zz,ref_zz):
        import numpy as np
        return np.sqrt(np.sum([((zz[i]-ref_zz['GAUSS'][i])/ref_zz['GAUSS_STD'][i])**2 for i in range(len(zz))]))

    def calc_zdist_med_std(self,zz,ref_zz):
        import numpy as np
        return np.sqrt(np.sum([((zz[i]-ref_zz['MED'][i])/ref_zz['MED_STD'][i])**2 for i in range(len(zz))]))

    def gfit_zerndist(self,zc_lst,n_bins):
        import numpy as np

        z_med   = np.median(zc_lst,axis=0)
        z_gpk   = []
        z_gsd   = []
        z_chi2  = []
        hists   = []
        cents   = []

        for zern_ord, z_init in enumerate(z_med):

            hist, bin_edges = np.histogram(zc_lst[:,zern_ord],bins=n_bins,range=(np.min(zc_lst[:,zern_ord]),np.max(zc_lst[:,zern_ord])))
            bin_centers     = (bin_edges[:-1] + bin_edges[1:])/2
            p_init          = [np.amax(hist), z_init, .01]
            [amp,gpk,gsd]   = self.gauss_fit(hist,bin_centers,p_init)
    	#print amp
            gsd             = np.abs(gsd)
            g_dist          = amp*np.exp((-(bin_centers-gpk)**2)/(2*gsd**2))
            chi2            = 0
            for i in range(len(hist)):
                    if g_dist[i]>1:
                        chi2            = chi2 + ((hist[i]-g_dist[i])**2)/g_dist[i]
            z_gpk.append(gpk)
            z_gsd.append(gsd)
            z_chi2.append(chi2)
            hists.append(hist)
            cents.append(bin_centers)
        return(z_gpk,z_gsd,z_chi2,hists,cents)

    def cat_comp_radec(self,cat1_list,cat2_list):
        import numpy as np

        cat1_cand   = []
        for cat2_item in cat2_list:
            for cat1_item in cat1_list:
                if np.abs(cat1_item['ALPHA_J2000']-cat2_item['ALPHA_J2000'])<0.0005 \
                        and np.abs(cat1_item['DELTA_J2000']-cat2_item['DELTA_J2000']) < 0.0005:
                    cat1_cand.append(cat1_item)
        return(cat1_cand)

    def cat_comp_xy(self,cat1_list,cat2_list):
        import numpy as np

        cat1_cand   = []
        for cat2_item in cat2_list:
            for cat1_item in cat1_list:
                if np.abs(cat1_item.x-cat2_item.x)<4.0 and np.abs(cat1_item.y-cat2_item.y)<4.0:
                    cat1_cand.append(cat1_item)
        return(cat1_cand)

    def gauss(self,x, *p):
        import numpy as np
        A, mu, sigma = p
        return A*np.exp(-(x-mu)**2/(2.*sigma**2))

    def gauss_fit(self,hist,bin_centers,p0):
        from scipy.optimize import curve_fit
        coeff, var_matrix = curve_fit(self.gauss, bin_centers, hist, p0=p0)
        return coeff
