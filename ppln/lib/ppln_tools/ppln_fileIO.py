def read_zern_file(zern_file):
    import sys
    sys.path.append('/Users/kack0001/astro-tools/goto/ppln_goto/lib/')
    from ppln_tools.ppln_classes import zern_obj

    zern_list   =   []
    with open(zern_file,'r') as zernID:
        for zern_line in zernID:
            line    = zern_line.split()
            if line[2] == "SUB_DIM":
                sub_dim = int(line[1])
            if line[2] == "SURR_DIM":
                sub_dim = int(line[1])
            if line[2] == "DISK_DIM":
                disk_dim = int(line[1])
            if line[2] == "DISK_RAD":
                disk_rad = float(line[1])
            if line[2] == "J_MAX":
                j_max = int(line[1])
            if line[0] == "#":
                if line[2] == 'Z_ORD':
                    ord_ind         = int(line[1])
                elif line[2] == 'MEAN':
                    mean_ind        = int(line[1])
                elif line[2] == 'MN_SD':
                    mean_std_ind    = int(line[1])
                elif line[2] == 'MED':
                    med_ind         = int(line[1])
                elif line[2] == 'MD_SD':
                    med_std_ind     = int(line[1])
                elif line[2] == 'GAUSS':
                    gauss_ind       = int(line[1])
                elif line[2] == 'GS_SD':
                    gauss_std_ind   = int(line[1])
                elif line[2] == 'CHI2':
                    gauss_chi2_ind  = int(line[1])
            else:
                z_ord           = int(line[ord_ind])
                z_mean          = float(line[mean_ind])
                z_mean_std      = float(line[mean_std_ind])
                z_med           = float(line[med_ind])
                z_med_std       = float(line[med_std_ind])
                z_gauss         = float(line[gauss_ind])
                z_gauss_std     = float(line[gauss_std_ind])
                z_gauss_chi2    = float(line[gauss_chi2_ind])
                if zern_list:
                    zern_list.z_ord.append(z_ord)
                    zern_list.z_mean.append(z_mean)
                    zern_list.z_mean_std.append(z_mean_std)
                    zern_list.z_med.append(z_med)
                    zern_list.z_med_std.append(z_med_std)
                    zern_list.z_gauss.append(z_gauss)
                    zern_list.z_gauss_std.append(z_gauss_std)
                    zern_list.z_gauss_chi2.append(z_gauss_chi2)
                else:
                    zern_list=zern_obj(z_ord, z_mean, z_mean_std, z_med, z_med_std, z_gauss, z_gauss_std, z_gauss_chi2)
    return zern_list,sub_dim,disk_dim,disk_rad,j_max
