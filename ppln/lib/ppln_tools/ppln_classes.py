class zern_obj:
    def __init__ (self,z_ord,z_mean,z_mean_std,z_med,z_med_std,z_gauss,z_gauss_std,z_gauss_chi2):
        self.z_ord          = [z_ord]
        self.z_mean         = [z_mean]
        self.z_mean_std     = [z_mean_std]
        self.z_med          = [z_med]
        self.z_med_std      = [z_med_std]
        self.z_gauss        = [z_gauss]
        self.z_gauss_std    = [z_gauss_std]
        self.z_gauss_chi2   = [z_gauss_chi2]

class trans_obj:
    def __init__ (self,x,y,snr,z_dist,mag_inst,mag_inst_err,ra,dec,flux,flux_err,mag_app,mag_app_err):
        self.x              = x
        self.y              = y
        self.snr            = snr
        self.z_dist         = z_dist
        self.mag_inst       = mag_inst
        self.mag_inst_err   = mag_inst_err
        self.ra             = ra
        self.dec            = dec
        self.flux           = flux
        self.flux_err       = flux_err
        self.mag_app        = mag_app
        self.mag_app_err    = mag_app_err

class wf_obj:
    def __init__ (self,wf,snr,var):
        self.wf     = wf
        self.snr    = snr
        self.var    = var
