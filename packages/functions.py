# import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import math
import statistics
import lasio as lasio
from scipy.stats import gaussian_kde
import platform
import glob
import os
import warnings
warnings.filterwarnings("ignore")


# Rock properties
K_QUARTZ = 36.8  # GPa
MU_QUARTZ = 44  # GPa
K_CLAY = 15  # GPa
MU_CLAY = 5  # GPa

phi=np.linspace(0.01, 0.4)



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

