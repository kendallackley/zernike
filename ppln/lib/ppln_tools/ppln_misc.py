import numpy as np

### Image and subimage manipulation

def sub_image(img,x_int,y_int,sub_dim):

    return np.array(img[int(y_int-(sub_dim-1)/2):int(y_int+(sub_dim+1)/2),int(x_int-(sub_dim-1)/2):int(x_int+(sub_dim+1)/2)], copy=True)

def remove_background(img):

    return img-np.median(img)

def shift_image(img,dx,dy,interp_order):
    from scipy import ndimage

    return ndimage.interpolation.shift(img, [-dy,-dx], order=interp_order)

def shift_image_lanczos3(img,dx,dy):
    #needs some fine tuning. produces strange psfs
    temp_dim    = img.shape[0]+6
    img_temp    = np.zeros((temp_dim,temp_dim))
    img_temp[3:img.shape[0]+3,3:img.shape[0]+3]=img
    l_krnl      = np.array([[lanczos3_2d(dx-k,dy-l) for k in range(-2,4)] for l in range(-2,4)])
    l_krnl      = l_krnl/np.sum(l_krnl)
    img_new     = np.array([[np.sum(img_temp[3+j-2:3+j+4,3+i-2:3+i+4]*l_krnl) for i in range(img.shape[1])] for j in range(img.shape[0])])
    return img_new

def calc_wavefront(sub_img,disk_mask):

    sub_dim     = sub_img.shape[0]
    disk_dim    = disk_mask.shape[0]
    disk_img    = sub_img[int((sub_dim-1)/2-(disk_dim-1)/2):int((sub_dim-1)/2+(disk_dim+1)/2),int((sub_dim-1)/2-(disk_dim-1)/2):int((sub_dim-1)/2+(disk_dim+1)/2)]*disk_mask
    return disk_img/disk_img.sum()

def img2wf(img,x,y,sub_dim,disk_mask):

    x_int       = int(round(x))
    y_int       = int(round(y))
    sub_img     = remove_background(sub_image(img,x_int,y_int,sub_dim))
    return calc_wavefront(sub_img,disk_mask)

def img2wf_shift(img,x_int,y_int,dx,dy,sub_dim,disk_mask):

    sub_img         = remove_background(sub_image(img,x_int,y_int,sub_dim))
    sub_img_shift   = shift_image(sub_img,dx,dy,5)
    return calc_wavefront(sub_img_shift,disk_mask)

def img2wf_lanczos3(img,x,y,sub_dim,disk_mask):
    import numpy as np

    x_int       = int(np.floor(x))
    y_int       = int(np.floor(y))
    dx          = x-x_int
    dy          = y-y_int
    sub_img     = remove_background(sub_image(img,x_int,y_int,sub_dim))
    disk_dim    = disk_mask.shape
    disk_img    = np.zeros(disk_dim)
    pix_off     = (sub_dim-disk_dim[0])/2

    l_krnl      = np.array([[lanczos3_2d(dx-k,dy-l) for k in range(-2,4)] for l in range(-2,4)])

    l_krnl      = l_krnl/np.sum(l_krnl)
    disk_img    = np.array([[np.sum(sub_img[int(pix_off+j-2):int(pix_off+j+4),int(pix_off+i-2):int(pix_off+i+4)]*l_krnl) for i in range(disk_mask.shape[1])] for j in range(disk_mask.shape[0])])

    disk_img    = disk_img*disk_mask
    return disk_img/np.sum(disk_img)

def lanczos3_2d(x,y):
    import numpy as np

    return np.sinc(x)*np.sinc(x/3)*np.sinc(y)*np.sinc(y/3)

def calc_wf_err(wf_img,wf_rec,disk_mask):
    import numpy as np

    n_pix   = np.sum(disk_mask)
    return np.sqrt(np.sum(disk_mask*(wf_img - wf_rec)**2)/n_pix)

### Catalog filtering

def filter_cat_xy(cat_list,img_shape,sub_dim):

    x_min       = (sub_dim+1)/2
    x_max       = img_shape[1]-x_min
    y_min       = (sub_dim+1)/2
    y_max       = img_shape[0]-y_min
    print('X AND Y: ',x_min,x_max,y_min,y_max)
    filter_list = []
    for cat_item in cat_list:
        if (x_min < cat_item.x < x_max and y_min < cat_item.y < y_max):
            filter_list.append(cat_item)
    return(filter_list)

def filter_cat_flux(cat_list,flux_min,flux_max):

    return([cat_item for cat_item in cat_list if (flux_min < cat_item.flux <= flux_max)])


def filter_cat_mag(cat_list,mag_min,mag_max):

    return([cat_item for cat_item in cat_list if (mag_min < cat_item.mag <= mag_max)])


def filter_cat_cls(cat_list,cls_min,cls_max):

    return([cat_item for cat_item in cat_list if (cls_min < cat_item.cls <= cls_max)])


def filter_cat_flag(cat_list,flag_val):

    return([cat_item for cat_item in cat_list if (cat_item.flag == flag_val)])

def filter_cat_radec(cat_list,ra,dec):

    return([cat_item for cat_item in cat_list if (cat_item.ra != ra and cat_item.dec != dec)])

### Zernike calculations

def cat2zern(cat_list,img_data,surr_dim,disc_dim,chi2_thr,zern_list,grid_mask,cov_mat_inv):
    import ppln_zernike

    zc_list = []
    for cat_item in cat_list:
        wf_img      = img2wf(img_data,cat_item.x,cat_item.y,surr_dim,disc_dim)
        zc_item     = ppln_zernike.zern_fit(wf_img,zern_list,cov_mat_inv)
        wf_rec      = ppln_zernike.wf_comp_brief(zc_item,zern_list,grid_mask)
        if ppln_zernike.calc_chi2(wf_img,wf_rec,grid_mask)<chi2_thr:
            zc_list.append(zc_item)
    return(zc_list)


def cat2zern_wf_shift(cat_list,img_data,surr_dim,disc_dim,chi2_thr,zern_list,grid_mask,cov_mat_inv):
    import ppln_zernike

    zc_list     =   []
    for i, cat_item in enumerate(cat_list):
        dx          = cat_item.x - int(round(cat_item.x))
        dy          = cat_item.y - int(round(cat_item.y))
        wf_img      = img2wf_shift(img_data,cat_item.x,cat_item.y,dx,dy,surr_dim,disc_dim)
        zc_item     = ppln_zernike.zern_fit(wf_img,zern_list,cov_mat_inv)
        wf_rec      = ppln_zernike.wf_comp_brief(zc_item,zern_list,grid_mask)
        if ppln_zernike.calc_chi2(wf_img,wf_rec,grid_mask)<chi2_thr:
            zc_list.append(zc_item)
    return(zc_list)

### Gauss fitting

def calc_zdist_gauss(zz,ref_zz):
    import numpy as np
    return np.sqrt(np.sum([((zz[i]-ref_zz.z_gauss[i])/ref_zz.z_gauss_std[i])**2 for i in range(len(zz))]))


def calc_zdist_med_std(zz,ref_zz):
    import numpy as np
    return np.sqrt(np.sum([((zz[i]-ref_zz.z_med[i])/ref_zz.z_mean_std[i])**2 for i in range(len(zz))]))


def gfit_zerndist(zc_lst,n_bins):
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
        [amp,gpk,gsd]   = gauss_fit(hist,bin_centers,p_init)
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

def gauss(x, *p):
    import numpy as np

    A, mu, sigma = p
    return A*np.exp(-(x-mu)**2/(2.*sigma**2))


def gauss_fit(hist,bin_centers,p0):
    from scipy.optimize import curve_fit
    coeff, var_matrix = curve_fit(gauss, bin_centers, hist, p0=p0)
    return coeff
