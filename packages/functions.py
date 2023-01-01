# import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import math
import statistics
import lasio as lasio
from scipy.stats import gaussian_kde

# Rock properties
K_QUARTZ = 36.8  # GPa
MU_QUARTZ = 44  # GPa
K_CLAY = 15  # GPa
MU_CLAY = 5  # GPa


phi=np.linspace(0.01, 0.4)


# define basic styles for plotting log curves (sty0), sand (sty1) and shale (sty2)
sty0 = {'lw':1, 'color':'k', 'ls':'-'}
sty1 = {'marker':'o', 'color':'g', 'ls':'none', 'ms':6, 'mec':'none', 'alpha':0.5}
sty2 = {'marker':'o', 'color':'r', 'ls':'none', 'ms':6, 'mec':'none', 'alpha':0.5}


# import warnings
# warnings.filterwarnings("ignore")

import platform
import glob
import os
my_os = platform.system()
if my_os == 'Windows':
    path = ".\data"
elif my_os == 'Linux':
    path = './data'
elif my_os == 'Darwin':
    path = './data'

# get all paths and alphabetically ordered
paths = sorted(glob.glob(os.path.join(path, "*.las")))
print(paths)


def load():
    well_df = [0] * len(paths)

    for i in range(len(paths)):
        # read with lasio
        well = lasio.read(paths[i])
        # convert to dataframe
        df = well.df()
        # print(df.head())  # ok
        # in this dataframe, depth is positioned as index, not as column
        # so better to change depth index to column
        well_df[i] = df.reset_index()
        well_df[i].rename(columns={'DEPT': 'DEPTH'}, inplace=True)
        # replace null values (-9999.0) with NaN
        well_df[i] = well_df[i].replace(-9999.0, np.NaN)

    return well_df




def well_add_features(df):


    # well2 adding features

    df = vshale_from_gr(df)

    df["PHIE"] = (2.65 - df.RHOB) / (2.65 - 1.05)
    df["IP"] = df.VP * 1000 * df.RHOB
    df["IS"] = df.VS * 1000 * df.RHOB
    df["VPVS"] = df.VP / df.VS
    # df['K'] = df.DEN*(df.VP**2 - 4/3.*(df.VS)**2)  # ?? what is this bulk modulus ??

    df['sandy-shaly'] = np.where(df['VSH'] >= 0.35, 'shaly', 'sandy')

    shale = df.VSH.values  #
    # shale = df.VWCL.values

    sand = 1 - shale - df.PHIE.values
    shaleN = shale / (shale + sand)  # normalized shale and sand volumes
    sandN = sand / (shale + sand)

    # mineral mixture bulk and shear moduli, k0 and mu0  # ?? A mineral mixture BUT NOT dry rock bulk modulus??
    k_u, k_l, mu_u, mu_l, k0, mu0 = vrh([shaleN, sandN], [K_CLAY, K_QUARTZ], [MU_CLAY, MU_QUARTZ])

    df['K0'] = k0

    ## well2 should be
    # Facies I: Gravels => not used
    # Facies II: Thick bedded sandstones; IIa, IIb, IIc, IId
    # Facies III: Interbedded sandstone-shale
    # Facies IV: Silty shales and silt-laminated shale
    # Facies V: Pure shales
    # Facies VI: Chaotic deposits => not used

    conditions = [
        (df["DEPTH"].ge(2078.0) & df["DEPTH"].lt(2105.0)),
        (df["DEPTH"].ge(2143.2) & df["DEPTH"].lt(2154.1)),
        (df["DEPTH"].ge(2154.1) & df["DEPTH"].lt(2164.1)),
        (df["DEPTH"].ge(2168.1) & df["DEPTH"].lt(2184.1)),
        (df["DEPTH"].ge(2186.1) & df["DEPTH"].lt(2200.1)),
        (df["DEPTH"].ge(2254.0) & df["DEPTH"].lt(2300.1)),
    ]
    # facies = ["shale", "sltShale", "clnSand", "sltSand1", "sltSand2", "cemSand"]
    facies = [1, 2, 3, 4, 5, 6]  # == FCODES[1:]
    df["FACIES"] = np.select(conditions, facies)

    reservoir = [0, 0, 1, 1, 0, 1]
    df["RESERVOIR"] = np.select(conditions, reservoir)

    facies_labels = ["shale", "sltShale", "clnSand", "sltSand1", "sltSand2", "cemSand"]
    df["LABELS"] = np.select(conditions, facies_labels)
    # df["FACIES"] = np.select(conditions, facies)  # can I use FCODES[1:] ???

    facies_codes = [6, 0, 6, 1, 2, 6, 3, 6, 4, 6, 5, 6]  # for well log plot
    conditions = [
        (df["DEPTH"].ge(df.DEPTH.min()) & df["DEPTH"].lt(2078.0)),  # undef=0    6
        (df["DEPTH"].ge(2078.0) & df["DEPTH"].lt(2105.0)),  # shale=1    0
        (df["DEPTH"].ge(2105.0) & df["DEPTH"].lt(2143.2)),  # undef=0    6
        (df["DEPTH"].ge(2143.2) & df["DEPTH"].lt(2154.1)),  # sltShale=2 1
        (df["DEPTH"].ge(2154.1) & df["DEPTH"].lt(2164.1)),  # clnSand=3  2
        (df["DEPTH"].ge(2164.1) & df["DEPTH"].lt(2168.1)),  # undef=0    6
        (df["DEPTH"].ge(2168.1) & df["DEPTH"].lt(2184.1)),  # sltSand1=4 3
        (df["DEPTH"].ge(2184.1) & df["DEPTH"].lt(2186.1)),  # undef=0    6
        (df["DEPTH"].ge(2186.1) & df["DEPTH"].lt(2200.1)),  # sltSand2=5 4
        (df["DEPTH"].ge(2200.1) & df["DEPTH"].lt(2254.0)),  # undef=0    6
        (df["DEPTH"].ge(2254.0) & df["DEPTH"].lt(2300.1)),  # cemSand=6  5
        (df["DEPTH"].ge(2300.1) & df["DEPTH"].lt(df.DEPTH.max()))  # undef=0    6
    ]
    df["FCODES"] = np.select(conditions, facies_codes)

    return df



def vshale_from_gr(df):
    """
    Creates Clavier, Larionov old, Larionov new, Steiber VSH
    """    
    GR_min = df.GR.min()
    GR_max = df.GR.max()
    df.loc[:, 'IGR'] = (df.GR - GR_min) / (GR_max - GR_min)
    df.loc[:, 'VSH_clavier'] = 1.7 - ((3.38 - (df.IGR + 0.7)**2)**0.5)
    df.loc[:, 'VSH_larionovO'] = 0.33 * (2**(2*df.IGR)-1)
    df.loc[:, 'VSH_steiber'] = df.IGR / (3 - 2*df.IGR)
    df.loc[:, 'VSH_larionovT'] = 0.083*(2**(3.7*df.IGR)-1)
    # Pick one to be "main" VSH
    df['VSH'] = df.VSH_larionovO
    return df


def kde_resample(df, col, lithologies, logs,  num_samples=1000):
    """
    Resample from each litho-facies a set number of times
    Inputs:
    df: well data in pandas dataframe with column of litho-facies and columns of Vp, Vs and rho logs
    col: column containing the litho-facies categories
    lithologies: the categories (could use column and .unique() instead of input list)
    logs: list of column names of Vp, Vs, rho (e.g. 'RHOB' or 'DEN', etc)
    num_samples: number of resamples from gaussian_kde.resample()
    """
    
    kde_resample_array = [0] * len(lithologies)
    
    for i in range(len(lithologies)):
        lith_logs = df[df[col]==lithologies[i]][logs[0]],\
                        df[df[col]==lithologies[i]][logs[1]],\
                        df[df[col]==lithologies[i]][logs[2]]
        # reset_index() ??
        kde_lith_logs = gaussian_kde(lith_logs)
        kde_resample_array[i] = kde_lith_logs.resample(num_samples)
            
    return kde_resample_array


def r0g(vp1, vs1, rho1, vp2, vs2, rho2):
    
    R0 = []
    G = []
    R_theta = []

    for i in range(len(vp1)):

        delta_rho = rho2[i] - rho1[i]
        delta_vp = vp2[i] - vp1[i]
        delta_vs = vs2[i] - vs1[i]

        vp = (vp2[i] + vp1[i]) / 2
        vs = (vs2[i] + vs1[i]) / 2
        rho = (rho2[i] + rho1[i]) / 2

        R0_temp = 1/2*((delta_vp / vp) + (delta_rho / rho))
        R0.append(R0_temp)

        G_temp = 1/2 * (delta_vp / vp) - 2 * vs**2/vp**2 * ((delta_rho / rho) + (2*delta_vs / vs))
        G.append(G_temp)

        R_theta_i = []
        for theta in range(41):
            R_theta_i.append(R0_temp + G_temp*(np.sin(math.radians(theta))**2))
        R_theta.append([R_theta_i])        

    median = []

    for j in range(41):
        angle = []
        for k in range(len(vp2)):
            angle.append(R_theta[k][0][j])
        median.append(statistics.median(angle))
        
    return R0, G, R_theta, median


def vrh(volumes,k,mu):
    f = np.array(volumes).T
    k = np.resize(np.array(k),np.shape(f))
    mu = np.resize(np.array(mu),np.shape(f))

    k_u = np.sum(f*k, axis=1)
    k_l = 1. / np.sum(f/k, axis=1)
    mu_u = np.sum(f*mu, axis=1)
    mu_l = 1. / np.sum(f/mu, axis=1)
    k0 = (k_u+k_l) / 2.
    mu0 = (mu_u+mu_l) / 2.
    return k_u, k_l, mu_u, mu_l, k0, mu0




# Fluid Replacement Modelling

def xfrm(vp1, vs1, rho1, rho_f1, k_f1, rho_f2, k_f2, k0, phi):
    """
    INPUT
    vp1, vs1, rho1: (vector) Measured Vp, Vs and density saturated with fluid 1
    rho_f1, k_f1:   (vector) Density and bulk modulus of fluid 1 (requires Sw)
    rho_f2, k_f2:   (scalar) Density and bulk modulus of fluid 2 (rho_o, k_o, rho_g, k_g, etc)
    k0:             (scalar) mineral bulk modulus - (k_u - k_l)/2
    phi:            (vector) porosity
    
    RETURN 
    vp2, vs2, rho2, k_s2: Vp, Vs, density and bulk modulus of rock with fluid 2. 
    
    Velocities are in m/s and densities in g/cm3.
    
    USAGE:
    vp1ox, vs1ox, rho1ox, k1ox = frm(vp1, vs1, rho1, RHO_WATER, K_WATER, RHO_GAS, K_GAS, K0, phiex1)
    """
    # convert Vp,Vs from m/s to km/s for calculation
    vp1  = vp1 / 1000.
    vs1  = vs1 / 1000.
    
    mu1  = rho1 * vs1**2.
    k_s1 = rho1 * vp1**2 - (4./3.)*mu1  # mu1 = rho1 * vs1**2
    
    # The dry rock bulk modulus
    kdry = (k_s1 * ((phi*k0)/k_f1 + 1 - phi) - k0) / ((phi*k0)/k_f1 + (k_s1/k0) - 1 - phi)
    
    # Now we can apply Gassmann to get the new values
    k_s2 = kdry + (1- (kdry/k0))**2 / ( (phi/k_f2) + ((1-phi)/k0) - (kdry/k0**2) )
    rho2 = rho1 - phi*rho_f1 + phi*rho_f2
    mu2  = mu1
    vp2  = np.sqrt(((k_s2 + (4./3)*mu2))/rho2)
    vs2  = np.sqrt((mu2/rho2))

    # return Vp,Vs as m/s
    return vp2*1000, vs2*1000, rho2, k_s2


def frm(vp_1, vs_1, rho_1, rho_f1, k_f1, rho_f2, k_f2, K0, phi_):
    """
    INPUT
    vp1, vs1, rho1: (vector) Measured Vp, Vs and density saturated with fluid 1
    rho_f1, k_f1:   (vector) Density and bulk modulus of fluid 1 (requires Sw)
    rho_f2, k_f2:   (scalar) Density and bulk modulus of fluid 2 (rho_o, k_o, rho_g, k_g, etc)
    k0:             (scalar) mineral bulk modulus - (k_u - k_l)/2
    phi:            (vector) porosity

    RETURN
    vp2, vs2, rho2, k_s2: Vp, Vs, density and bulk modulus of rock with fluid 2.

    Velocities are in m/s and densities in g/cm3.

    USAGE:
    vp1ox, vs1ox, rho1ox, k1ox = frm(vp1, vs1, rho1, RHO_WATER, K_WATER, RHO_GAS, K_GAS, K0, phiex1)
    """

    vp_frm = []
    vs_frm = []
    rho_frm = []

    for i in range(0, len(vp_1)):

        # convert Vp,Vs from m/s to km/s for calculation
        vp1 = vp_1[i] / 1000.
        vs1 = vs_1[i] / 1000.
        rho1 = rho_1[i]
        phi = phi_[i]
        k0 = K0[i]

        mu1 = rho1 * vs1 ** 2.
        rho2 = rho1 - phi * rho_f1 + phi * rho_f2
        rho_frm.append(rho2)
        mu2 = mu1
        vs2 = np.sqrt((mu2 / rho2))
        vs_frm.append(vs2*1000)

        k_s1 = rho1 * vp1 ** 2 - (4. / 3.) * mu1  # mu1 = rho1 * vs1**2    

        # The dry rock bulk modulus
        kdry = (k_s1 * ((phi*k0)/k_f1+1-phi)-k0) / ((phi*k0)/k_f1+(k_s1/k0)-1-phi)

        # kdry values were negative when k0=36 was too high -> use vrh()
        # Vp2 values were getting > 6000 when kdry was > 70 (why so high?)
        if (0 < kdry < 36):

            # Now we can apply Gassmann to get the new values
            k_s2 = kdry + (1 - (kdry / k0)) ** 2 / ((phi / k_f2) + ((1 - phi) / k0) - (kdry / k0 ** 2))
            
            if k_s2 > 0:
                
            # print(f"k_s2: {k_s2}")
                vp2 = np.sqrt(((k_s2 + (4. / 3.) * mu2)) / rho2)
                # vp_frm.append(np.round(vp2*1000, 0))
                vp_frm.append(int(vp2*1000))
            else: 
                vp2 = vp1
                vp_frm.append(vp2*1000)

        else:
            vp2 = vp1
            vp_frm.append(vp2*1000)
            
    # return as arrays
    vp_frm = np.array(vp_frm)
    vs_frm = np.array(vs_frm)
    rho_frm = np.array(rho_frm)
    
    return vp_frm, vs_frm, rho_frm, k_s2


