import numpy as np
from scipy.special import factorial as fac

def noll2mn(j):
    """
    Convert to Noll index ordering
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

def pre_fac(k,m,n):
    """
    Returns the prefactor for the relative weighting of the rho terms in R_mn
    """
    return -2*(k%2-0.5)*fac(n-k)/(fac(k)*fac((n+m)/2-k)*fac((n-m)/2-k))

def R_nm(m,n,rho):
    """
    Returns the radial part of Z_mn on the grid rho
    """
    return sum(pre_fac(k,m,n)*rho**(n-2*k) for k in range(int(((n-m)/2)+1)))

def Z_nm(m,n,rho,phi):
    """
    Returns the full Z_mn on the coordinate grid rho,phi
    """
    if (m >= 0): return R_nm(m,n,rho)*np.cos(m*phi)
    if (m < 0): return R_nm(-m,n,rho)*np.sin(-m*phi)

def avg_subarr(arr,n_ds):
    h,w         = arr.shape
    h           = (h//n_ds)*n_ds
    w           = (w//n_ds)*n_ds
    arr         = arr[:h,:w]
    return np.einsum('ijkl->ik',arr.reshape(h//n_ds,n_ds,-1,n_ds))/n_ds**2

def create_coord_fgrid(dim,r,dx,dy,n_us):
    """
    Creates a finer grid for interpolation
    """
    fdim        = dim*n_us
    fgrid       = (np.indices((fdim,fdim),dtype=np.float)-((fdim-1)/2))/n_us
    fgrid[0]    = (fgrid[0]+dy)/r
    fgrid[1]    = (fgrid[1]+dx)/r
    fgrid_rho   = (fgrid**2).sum(0)**0.5
    fgrid_phi   = np.arctan2(fgrid[0],fgrid[1])
    fgrid_mask  = (fgrid_rho <= 1)*1.0
    fgrid_mask_ds   = avg_subarr(fgrid_mask,n_us)
    return fgrid_rho,fgrid_phi,fgrid_mask,fgrid_mask_ds

def create_zern_list(z_ord,dim,r_disc,dx,dy):
    grid_rho,grid_phi,grid_mask = create_coord_grid(dim,r_disc,dx,dy)
    mn_list = [noll2mn(j) for j in z_ord]
    return [Z_nm(mn[0],mn[1],grid_rho,grid_phi)*grid_mask for mn in mn_list],grid_mask

def create_fzern_list(z_ord,dim,r_disc,dx,dy,n_us):
    fgrid_rho,fgrid_phi,fgrid_mask,fgrid_mask_ds = create_coord_fgrid(dim,r_disc,dx,dy,n_us)
    mn_list = [noll2mn(j) for j in z_ord]
    return [avg_subarr(Z_nm(mn[0],mn[1],fgrid_rho,fgrid_phi)*fgrid_mask,n_us) for mn in mn_list],fgrid_mask_ds

def inv_cov_mat(zern_list):
    cov_mat = np.array([[np.sum(zerni*zernj) for zerni in zern_list] for zernj in zern_list])
    return np.linalg.pinv(cov_mat)

def zern_fit(wf_img,zern_list,cov_mat_inv):
    """
    Returns the reconstructed PSF
    """
    # Calculate the inner product of each Zernike mode with the test surface
    wf_zern_inprod  = np.array([np.sum(wf_img*zerni) for zerni in zern_list])
    # Given the inner product vector of the test wavefront with Zernike basis,
    # calculate the Zernike polynomial coefficients
    wf_rec_pow      = np.dot(cov_mat_inv,wf_zern_inprod)
    return wf_rec_pow

def wf_comp_brief(z_vec,zern_list):
    wf_temp = sum(val*zern_list[i] for (i, val) in enumerate(z_vec))
    return wf_temp/wf_temp.sum()

def wf_comp(z_vec,dim):
    zern_list,mask   = create_zern_list(len(z_vec),dim)
    return wf_comp_brief(z_vec,zern_list)

def calc_wf_error(wf_img,wf_rec,mask):
    n_pix   = np.sum(mask)
    return np.sqrt(np.sum(mask*(wf_img-wf_rec)**2)/n_pix)

def calc_wf_var(wf_img,wf_rec,wf_grid):
    import numpy as np

    n_pix   = np.sum(wf_grid)
    return np.sum((wf_img - wf_rec)**2)/n_pix

def calc_chi2(subim,wf,mask):
    return sum(((subim[i,j]-wf[i,j])**2/wf[i,j]) if mask[i,j] else 0 for i in range(int(len(subim[0]))) for j in range(len(subim)))
