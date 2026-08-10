# -*- coding: utf-8 -*-
"""EAM/alloy potential reader and analytical cohesive-energy statistics.

Equation references
-------------------
[Vacancy manuscript] A. Baski et al., "An Analytical Method for
Quantifying Vacancy Energetics and Vacancy Transport Behavior in
Concentrated Solid Solutions" (revised manuscript supplied with this code).

[CMS 2022] R. Jagatramka, C. Wang, and M. Daly, Computational Materials
Science 214 (2022) 111763, doi:10.1016/j.commatsci.2022.111763.

The calculation follows the coordination-shell reparameterization of the
embedded-atom method (EAM). Warren-Cowley short-range-order parameters are
introduced through P_zeta^(XY) = C_X C_Y (1-alpha_zeta^(XY)); see vacancy
manuscript Eqs. (5)-(12). Numerical expressions are intentionally unchanged.
"""

def potential_read(fname='NiCo-lammps-2014.alloy'):
    """Read a LAMMPS EAM/alloy ``setfl`` file.

    Returns distance/density grids and tabulated electron-density, embedding,
    and pair-potential functions. Pair data are expanded to a symmetric
    ``(Nr, n_elements, n_elements)`` array.
    """
    import numpy as np
    with open(fname) as f:
        lines=f.readlines()
    chem=lines[3]
    chem=int(str.split(chem)[0])
    Nrho=int(str.split(lines[4])[0])
    drho=float(str.split(lines[4])[1])
    Nr=int(str.split(lines[4])[2])
    dr=float(str.split(lines[4])[3])
    cols=np.shape(str.split(lines[6]))[0]
    
    Nrho_rows=int(Nrho/cols)
    Nr_rows=int(Nr/cols)
    
    rho=np.zeros((Nr,chem))
    Fr=np.zeros((Nrho,chem))
    Pp=np.zeros((Nr,sum(np.arange(0,chem+1))))
    
    k=6
    for i in np.arange(0,chem):
        m=0
        for j in np.arange(0,Nrho):
            
            Fr[j,i]=float(str.split(lines[i*(Nr_rows-Nrho_rows)+i*Nrho_rows+k+i])[m])            
            m=m+1                
            if m==cols:
                m=0
                k=k+1    
       
    k=6
    # check1=[]
    for i in np.arange(0,chem):
        m=0
        for j in np.arange(0,Nr): 
            rho[j,i]=float(str.split(lines[(i+1)*Nrho_rows+k+i])[m])
            m=m+1                
            if m==cols:
                m=0
                k=k+1
    
    
    k=chem*Nr_rows+chem*Nrho_rows+5+chem  
    
    for i in np.arange(0,np.shape(Pp)[1]):
        m=0
        for j in np.arange(0,Nr):
            Pp[j,i]=float(str.split(lines[k])[m])
            m=m+1                
            if m==cols:
                m=0
                k=k+1
    
    lists=[]
    m=0
    for i in np.arange(0,chem):
        for j in np.arange(0,chem):
            if i < j:
                continue
            else:
                lists.append([i,j,m])
                m=m+1
    
    Pp_new=np.zeros((Nr,chem,chem))  
    for i in np.arange(0,np.shape(lists)[0]):
        ind1=lists[i][0]
        ind2=lists[i][1]
        ind3=lists[i][2]
        Pp_new[:,ind1,ind2]=Pp[:,ind3]
        
        if ind2<ind1:
            Pp_new[:,ind2,ind1]=Pp[:,ind3]
    
    Pp=Pp_new
    
    rrange=np.arange(0,Nr)*dr
    rhorange=np.arange(0,Nrho)*drho
#%%    
    return rrange,rhorange,rho,Fr,Pp
#%%
def potential_fetch(rrange,rhorange,rho,Fr,Pp,atoms,Neighbors,ind):
    """Evaluate the conventional atom-resolved EAM terms for one atom.

    This is the direct site form of vacancy-manuscript Eq. (1), with the
    local electron density defined by Eq. (2).
    """
#%%
    import numpy as np
    
    rho_a=0
    Fr_a=0
    Pp_a=0    
    
    for i in np.arange(1,np.shape(Neighbors[1])[1]):    
        #print(i)
        r=Neighbors[0][ind,i]
        cent=int(atoms[ind,3]-1)
        pair=int(atoms[Neighbors[1][ind,i],3]-1)        
        rho_a=rho_a+np.interp(r,rrange,rho[:,pair])        
        Pp_a=Pp_a+(0.5*np.interp(r,rrange,Pp[:,cent,pair]))/r
        
    Fr_a=np.interp(rho_a,rhorange,Fr[:,cent])
 #%%
    return rho_a, Fr_a, Pp_a
# #%%
# def potential_stats(rrange, rhorange, rho, Fr, Pp, comp, cn, alpha):
#     """Return mean/std cohesive-energy statistics for one environment.

#     ``cn[k] = (r_k, N_k)`` defines the coordination-shell structure factor.
#     ``alpha[k,i,j]`` is the Warren-Cowley parameter for shell k and pair i-j.

#     Main references: vacancy-manuscript Eqs. (3)-(12); the random-alloy
#     limit is equivalent to CMS-2022 Eqs. (3)-(10).
#     """
#     import numpy as np
#     import pandas as pd
    
#     if alpha.shape[0] < cn.shape[0]:
#         pad_len = cn.shape[0] - alpha.shape[0]
#         pad = np.zeros((pad_len, alpha.shape[1], alpha.shape[2]))
#         alpha = np.concatenate((alpha, pad), axis=0)

#     rho_bar = 0
#     rho_cn = np.zeros((np.shape(cn)[0], np.shape(comp)[0]))

#     for j in np.arange(0, np.shape(comp)[0]):
#         for i in np.arange(0, np.shape(cn)[0]):
#             c = comp[j]
#             r = cn[i, 0]
#             m = cn[i, 1]
#             Rho = np.interp(r, rrange, rho[:, j])
#             # Mean local electron density, vacancy manuscript Eq. (3)
#             # (and the charge-density contribution underlying Eq. (1)).
#             rho_bar = rho_bar + c * m * Rho
#             rho_cn[i, j] = np.interp(r, rrange, rho[:, j])

#     rho_std = np.zeros((np.shape(rho_cn)[0]))
#     for i in np.arange(0, np.shape(rho_std)[0]):
#         rho_std[i] = np.sum((rho_cn[i, :] - np.sum(rho_cn[i, :] * comp))**2 * comp) * cn[i, 1]
#         rho_std[i] = np.sqrt(rho_std[i])

#     F_bar = 0
#     Fs = np.zeros((np.shape(comp)[0]))
#     for j in np.arange(0, np.shape(comp)[0]):
#         c = comp[j]
#         Fs[j] = np.interp(rho_bar, rhorange, Fr[:, j])
#         F_bar = F_bar + Fs[j] * c

#     # Standard deviation of embedding energy: vacancy manuscript Eq. (10)
#     # embedding contribution; CMS-2022 Eq. (5).
#     F_std = np.sqrt(np.sum((Fs - F_bar)**2 * comp))

#     Pp_bar = np.zeros((np.shape(comp)[0], np.shape(comp)[0]))
#     Pp_cn = np.zeros((np.shape(cn)[0], np.shape(comp)[0], np.shape(comp)[0]))

#     for i in np.arange(0, np.shape(Pp_bar)[0]):
#         for j in np.arange(0, np.shape(Pp_bar)[1]):
#             Pp_ind = 0
#             for k in np.arange(0, np.shape(cn)[0]):
#                 Pp_cn[k, i, j] = np.interp(cn[k, 0], rrange, Pp[:, i, j]) / cn[k, 0]
#                 # Ordered pair probability C_i C_j (1-alpha_kij):
#                 # vacancy manuscript Eqs. (5)-(7).
#                 Pp_ind += comp[i]*comp[j]*cn[k,1]*Pp_cn[k,i,j]*(1-alpha[k,i,j])
#             Pp_bar[i, j] = Pp_ind

#     form_E = np.zeros((2, 4))
#     form_E[0, 0] = rho_bar
#     form_E[1, 0] = np.sqrt(np.sum(rho_std**2))
#     form_E[0, 1] = F_bar
#     form_E[1, 1] = F_std
#     # Mean pair contribution (1/2 avoids double counting):
#     # vacancy manuscript Eq. (9); CMS-2022 Eq. (3).
#     form_E[0, 2] = np.sum(Pp_bar) * 0.5

#     Pp_std_cn = np.zeros((np.shape(cn)[0], np.shape(comp)[0]))
#     Pp_std_avg = np.zeros((np.shape(cn)[0], np.shape(comp)[0]))
#     for k in np.arange(0, np.shape(cn)[0]):
#         for j in np.arange(0, np.shape(comp)[0]):
#             # Central-species/shell mean pair interaction. This is the SRO
#             # extension of CMS-2022 Appendix-B Eq. (B.2).
#             avg = np.sum(Pp_cn[k, j, :] * comp*(1-alpha[k,j,:]))
#             Pp_std_avg[k, j] = avg * cn[k, 1]
#             Pp_std_cn[k, j] = np.sum((Pp_cn[k, j, :] - avg)**2 * comp*(1-alpha[k,i,j])) * cn[k, 1]

#     Pp_std_cn2 = np.sqrt(np.sum(Pp_std_cn, axis=0))
#     Pp_std_avg2 = np.sum(Pp_std_avg, axis=0)

#     # Pair-interaction standard deviation: vacancy manuscript Eqs. (11)-(12);
#     # random-alloy counterpart CMS-2022 Eqs. (6)-(8).
#     Pp_std = np.sqrt(np.sum(comp * (Pp_std_cn2**2 + (Pp_std_avg2 - np.sum(Pp_bar))**2)))

#     form_E[1, 2] = Pp_std * 0.5
#     # Mean per-atom cohesive/site energy: vacancy manuscript Eq. (9);
#     # CMS-2022 Eq. (3).
#     form_E[0, 3] = F_bar + np.sum(Pp_bar) * 0.5

#     # Covariance cov(F, V/2): vacancy manuscript Eq. (10);
#     # CMS-2022 Eq. (10).
#     exp_FV = np.sum(comp * Fs * Pp_std_avg2 * 0.5) - F_bar * np.sum(Pp_bar) * 0.5
#     covar2 = exp_FV
#     # Total site-energy standard deviation from variance propagation:
#     # vacancy manuscript Eq. (10); CMS-2022 Eq. (9).
#     form_E[1, 3] = np.sqrt(F_std**2 + (Pp_std**2 / 4) + 2 * covar2)

#     form_E = pd.DataFrame(form_E, columns=["rho", "F", "Pp", "E"])
#     form_E.index = ["Mean", "Std"]

#     covars = np.array([covar2])
#     test = [Fs + Pp_std_avg2 * 0.5]

#     return form_E, test, covars


#%%
def potential_stats(rrange, rhorange, rho, Fr, Pp, comp, cn, alpha):
    """Return mean/std cohesive-energy statistics for one environment.

    ``cn[k] = (r_k, N_k)`` defines the coordination-shell structure factor.
    ``alpha[k,i,j]`` is the Warren-Cowley parameter for shell k and pair i-j.

    Main references: vacancy-manuscript Eqs. (3)-(12); the random-alloy
    limit is equivalent to CMS-2022 Eqs. (3)-(10).
    """
    import numpy as np
    import pandas as pd
    
    if alpha.shape[0] < cn.shape[0]:
        pad_len = cn.shape[0] - alpha.shape[0]
        pad = np.zeros((pad_len, alpha.shape[1], alpha.shape[2]))
        alpha = np.concatenate((alpha, pad), axis=0)

    # ------------------------------------------------------------------
    # Mean electron density with SRO
    #
    # rho_bar[j] is the mean electron density for central species j:
    #
    # rho_bar_X = sum_eta sum_Y
    #             N_eta*C_Y*(1-alpha_eta_XY)*rho_eta_Y
    # ------------------------------------------------------------------
    rho_bar = np.zeros((np.shape(comp)[0]))
    rho_cn = np.zeros((np.shape(cn)[0], np.shape(comp)[0]))

    # Electron-density contribution of each neighboring species
    # at each coordination-shell distance
    for j in np.arange(0, np.shape(comp)[0]):
        for i in np.arange(0, np.shape(cn)[0]):
            r = cn[i, 0]
            rho_cn[i, j] = np.interp(r, rrange, rho[:, j])

    # Mean rho for each possible central species
    # j = central species
    # k = neighboring species
    # i = coordination shell
    for j in np.arange(0, np.shape(comp)[0]):
        for i in np.arange(0, np.shape(cn)[0]):
            m = cn[i, 1]

            for k in np.arange(0, np.shape(comp)[0]):
                c = comp[k]
                Rho = rho_cn[i, k]

                rho_bar[j] = rho_bar[j] + c * m * Rho * (1 - alpha[i, j, k])

    # ------------------------------------------------------------------
    # Standard deviation of rho with SRO
    #
    # rho_std[j] = sigma_(rho|X)
    # ------------------------------------------------------------------
    rho_std = np.zeros((np.shape(comp)[0]))

    for j in np.arange(0, np.shape(comp)[0]):

        rho_var = 0

        for i in np.arange(0, np.shape(cn)[0]):

            rho_avg = np.sum(
                rho_cn[i, :] * comp * (1 - alpha[i, j, :])
            )

            rho_var = rho_var + np.sum(
                (rho_cn[i, :] - rho_avg)
                * (rho_cn[i, :] - rho_avg)
                * comp
                * (1 - alpha[i, j, :])
            ) * cn[i, 1]

        rho_std[j] = np.sqrt(rho_var)

    # Overall mean and standard deviation of rho
    #
    # Var(rho) = E_X[Var(rho|X)] + Var_X(E[rho|X])
    rho_bar_avg = np.sum(rho_bar * comp)

    rho_std_avg = np.sqrt(
        np.sum(
            comp * (
                rho_std**2
                + (rho_bar - rho_bar_avg)**2
            )
        )
    )

    # ------------------------------------------------------------------
    # Mean and standard deviation of embedding energy F
    #
    # Each central species uses its own SRO-dependent rho_bar[j].
    # sigma_(F|X) = |dF_X/drho| * sigma_(rho|X)
    # ------------------------------------------------------------------
    F_bar = 0
    Fs = np.zeros((np.shape(comp)[0]))
    F_std_rho = np.zeros((np.shape(comp)[0]))

    for j in np.arange(0, np.shape(comp)[0]):
        c = comp[j]

        Fs[j] = np.interp(rho_bar[j], rhorange, Fr[:, j])
        F_bar = F_bar + Fs[j] * c

        # Derivative dF/drho from the tabulated embedding function
        dF = np.gradient(Fr[:, j], rhorange)

        F_prime = np.interp(
            rho_bar[j],
            rhorange,
            dF
        )

        # Embedding-energy standard deviation caused by rho fluctuations
        F_std_rho[j] = np.abs(F_prime) * rho_std[j]

    # Overall standard deviation of F
    #
    # Var(F) = E_X[Var(F|X)] + Var_X(E[F|X])
    F_std = np.sqrt(
        np.sum(
            comp * (
                F_std_rho**2
                + (Fs - F_bar)**2
            )
        )
    )

    Pp_bar = np.zeros((np.shape(comp)[0], np.shape(comp)[0]))
    Pp_cn = np.zeros((np.shape(cn)[0], np.shape(comp)[0], np.shape(comp)[0]))

    for i in np.arange(0, np.shape(Pp_bar)[0]):
        for j in np.arange(0, np.shape(Pp_bar)[1]):
            Pp_ind = 0
            for k in np.arange(0, np.shape(cn)[0]):
                Pp_cn[k, i, j] = np.interp(cn[k, 0], rrange, Pp[:, i, j]) / cn[k, 0]
                # Ordered pair probability C_i C_j (1-alpha_kij):
                # vacancy manuscript Eqs. (5)-(7).
                Pp_ind += comp[i]*comp[j]*cn[k,1]*Pp_cn[k,i,j]*(1-alpha[k,i,j])
            Pp_bar[i, j] = Pp_ind

    form_E = np.zeros((2, 4))
    form_E[0, 0] = rho_bar_avg
    form_E[1, 0] = rho_std_avg
    form_E[0, 1] = F_bar
    form_E[1, 1] = F_std
    # Mean pair contribution (1/2 avoids double counting):
    # vacancy manuscript Eq. (9); CMS-2022 Eq. (3).
    form_E[0, 2] = np.sum(Pp_bar) * 0.5

    Pp_std_cn = np.zeros((np.shape(cn)[0], np.shape(comp)[0]))
    Pp_std_avg = np.zeros((np.shape(cn)[0], np.shape(comp)[0]))
    for k in np.arange(0, np.shape(cn)[0]):
        for j in np.arange(0, np.shape(comp)[0]):
            # Central-species/shell mean pair interaction. This is the SRO
            # extension of CMS-2022 Appendix-B Eq. (B.2).
            avg = np.sum(Pp_cn[k, j, :] * comp*(1-alpha[k,j,:]))
            Pp_std_avg[k, j] = avg * cn[k, 1]
            Pp_std_cn[k, j] = np.sum((Pp_cn[k, j, :] - avg)**2 * comp*(1-alpha[k,i,j])) * cn[k, 1]

    Pp_std_cn2 = np.sqrt(np.sum(Pp_std_cn, axis=0))
    Pp_std_avg2 = np.sum(Pp_std_avg, axis=0)

    # Pair-interaction standard deviation: vacancy manuscript Eqs. (11)-(12);
    # random-alloy counterpart CMS-2022 Eqs. (6)-(8).
    Pp_std = np.sqrt(np.sum(comp * (Pp_std_cn2**2 + (Pp_std_avg2 - np.sum(Pp_bar))**2)))

    form_E[1, 2] = Pp_std * 0.5
    # Mean per-atom cohesive/site energy: vacancy manuscript Eq. (9);
    # CMS-2022 Eq. (3).
    form_E[0, 3] = F_bar + np.sum(Pp_bar) * 0.5

    # Covariance cov(F, V/2): vacancy manuscript Eq. (10);
    # CMS-2022 Eq. (10).
    exp_FV = np.sum(comp * Fs * Pp_std_avg2 * 0.5) - F_bar * np.sum(Pp_bar) * 0.5
    covar2 = exp_FV
    # Total site-energy standard deviation from variance propagation:
    # vacancy manuscript Eq. (10); CMS-2022 Eq. (9).
    form_E[1, 3] = np.sqrt(F_std**2 + (Pp_std**2 / 4) + 2 * covar2)

    form_E = pd.DataFrame(form_E, columns=["rho", "F", "Pp", "E"])
    form_E.index = ["Mean", "Std"]

    covars = np.array([covar2])
    test = [Fs + Pp_std_avg2 * 0.5]

    return form_E, test, covars
