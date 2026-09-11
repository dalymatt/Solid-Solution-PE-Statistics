# -*- coding: utf-8 -*-
"""EAM/alloy potential reader and analytical cohesive-energy statistics.

Equation references
-------------------
[Vacancy manuscript] A. Baski et al., "A mechanistic model for vacancy 
    energetics in concentrated solid solutions with short-range order".

[CMS 2022] R. Jagatramka, C. Wang, and M. Daly, Computational Materials
Science 214 (2022) 111763, doi:10.1016/j.commatsci.2022.111763.

The calculation follows the coordination-shell reparameterization
of the embedded-atom method (EAM), with Warren-Cowley short-range-order
parameters entering through the conditional pair probabilities p(Y|X)
= c_Y (1 - alpha_XY); see vacancy manuscript Eqs. (1)-(9). 
Shell-resolved alpha values are mapped onto the peaks of each coordination
environment by distance.
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
def potential_stats(rrange, rhorange, rho, Fr, Pp, comp, cn, alpha):
    """Return mean/std cohesive-energy statistics for one environment.

    ``cn[k] = (r_k, N_k)`` defines the coordination-shell structure factor.
    ``alpha[k,i,j]`` is the Warren-Cowley parameter for shell k and pair i-j.
    """
    import numpy as np
    import pandas as pd
    
    alpha = np.asarray(alpha, dtype=float)
    cn = np.asarray(cn, dtype=float)
    
    # Validate array dimensions and Warren-Cowley conditional probabilities.
    nsp = np.shape(comp)[0]
    if alpha.shape[1:] != (nsp, nsp):
        raise ValueError(
            f"alpha must have shape (n_shells, {nsp}, {nsp}); got {alpha.shape}")
    if alpha.ndim != 3:
        raise ValueError(f"alpha must be three-dimensional; got {alpha.shape}")
    if not np.isclose(np.sum(comp), 1.0, atol=1e-10):
        raise ValueError(f"Composition must sum to 1.0; got {np.sum(comp):.12g}")
    if np.any(comp < 0):
        raise ValueError("Composition entries must be non-negative")

    # P(Y|X,zeta) = C_Y * (1-alpha_zeta^(XY)). For a physically valid
    # WC description each conditional distribution should be normalized and
    # non-negative. We raise an error instead of silently renormalizing.
    pcond_shell = comp[None, None, :] * (1.0 - alpha)
    if np.any(pcond_shell < -1e-12):
        bad = np.argwhere(pcond_shell < -1e-12)[0]
        raise ValueError(
            "Invalid Warren-Cowley parameters: negative conditional "
            f"probability at shell={bad[0]}, X={bad[1]}, Y={bad[2]}."
        )
    norm_shell = np.sum(pcond_shell, axis=2)
    max_norm_err = float(np.max(np.abs(norm_shell - 1.0)))
    if max_norm_err > 1e-6:
        raise ValueError(
            "Warren-Cowley parameters violate conditional-probability "
            f"normalization; max |sum_Y P(Y|X)-1| = {max_norm_err:.3e}."
        )    
        
    # Map shell-level alpha to every coordination peak.  The first peak in
    # the shipped FCC/vacancy/TS environments is a first-neighbour peak, so
    # it provides a robust local estimate of a0.  Peaks farther than the
    # supplied WC-shell range receive alpha=0 rather than inheriting an
    # unrelated shell by array position.
    
    fcc_ratios = np.array([
        1.0 / np.sqrt(2.0),
        1.0,
        np.sqrt(3.0 / 2.0),
        np.sqrt(2.0),
        np.sqrt(5.0 / 2.0),
        np.sqrt(3.0),
        np.sqrt(7.0 / 2.0),
        2.0,
    ])
    if alpha.shape[0] > fcc_ratios.size:
        raise ValueError(
            f"Distance-based WC mapping currently supports up to "
            f"{fcc_ratios.size} FCC shells; got {alpha.shape[0]}."
        )

    a0_est = cn[0, 0] / fcc_ratios[0]
    ratios = cn[:, 0] / a0_est
    mapped_alpha = np.zeros((cn.shape[0], nsp, nsp), dtype=float)
    shell_ratios = fcc_ratios[:alpha.shape[0]]
    for peak_index, ratio in enumerate(ratios):
        shell_index = int(np.argmin(np.abs(shell_ratios - ratio)))
        if abs(shell_ratios[shell_index] - ratio) < 0.06:
            mapped_alpha[peak_index] = alpha[shell_index]
    alpha = mapped_alpha

    pcond = comp[None, None, :] * (1.0 - alpha)
    if np.any(pcond < -1e-12):
        bad = np.argwhere(pcond < -1e-12)[0]
        raise ValueError(
            "Invalid Warren-Cowley parameters: negative conditional "
            f"probability at shell={bad[0]}, X={bad[1]}, Y={bad[2]}."
        )

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

            rho_avg = np.sum(rho_cn[i, :] * comp * (1 - alpha[i, j, :]))

            rho_var = rho_var + np.sum((rho_cn[i, :] - rho_avg) * (rho_cn[i, :] - rho_avg) * comp * (1 - alpha[i, j, :])) * cn[i, 1]

        rho_std[j] = np.sqrt(rho_var)

    # Overall mean and standard deviation of rho
    #
    # Var(rho) = E_X[Var(rho|X)] + Var_X(E[rho|X])
    rho_bar_avg = np.sum(rho_bar * comp)

    rho_std_avg = np.sqrt(np.sum( comp * (rho_std**2 + (rho_bar - rho_bar_avg)**2 )))

    # ------------------------------------------------------------------
    # Mean and standard deviation of embedding energy F
    #
    # Each central species uses its own SRO-dependent rho_bar[j].
    # Manuscript Eq. (8b) is used for the embedding-energy variance.
    # ------------------------------------------------------------------
    F_bar = 0
    Fs = np.zeros((np.shape(comp)[0]))

    for j in np.arange(0, np.shape(comp)[0]):
        c = comp[j]

        Fs[j] = np.interp(rho_bar[j], rhorange, Fr[:, j])
        F_bar = F_bar + Fs[j] * c

    # Manuscript Eq. (8b): neglect the intra-species embedding-energy
    # contribution and retain only the inter-species spread of
    # F^X(rho_bar^X). This is equivalent to
    # sum_{X<Y} C_X C_Y [F^X - F^Y]^2.
    # Equation (8b): Standard deviation of the embedding energy
    bb_F = np.zeros((np.shape(comp)[0], np.shape(comp)[0]))

    for i in np.arange(0, np.shape(bb_F)[0]):
        for j in np.arange(0, np.shape(bb_F)[0]):
            if i <= j:
                continue

            bb_F[i, j] = (comp[i] * comp[j] * np.square(Fs[i] - Fs[j]) )

    F_std = np.sqrt(np.sum(bb_F))
    
    #Pair Interactions

    Pp_bar = np.zeros((np.shape(comp)[0], np.shape(comp)[0]))
    Pp_cn = np.zeros((np.shape(cn)[0], np.shape(comp)[0], np.shape(comp)[0]))

    for i in np.arange(0, np.shape(Pp_bar)[0]):
        for j in np.arange(0, np.shape(Pp_bar)[1]):
            Pp_ind = 0
            for k in np.arange(0, np.shape(cn)[0]):
                Pp_cn[k, i, j] = np.interp(cn[k, 0], rrange, Pp[:, i, j]) / cn[k, 0]
                # Ordered pair probability C_i C_j (1-alpha_kij):
                # vacancy manuscript Eqs. (3)-(4).
                Pp_ind += comp[i]*comp[j]*cn[k,1]*Pp_cn[k,i,j]*(1-alpha[k,i,j])
            Pp_bar[i, j] = Pp_ind

    form_E = np.zeros((2, 4))
    form_E[0, 0] = rho_bar_avg
    form_E[1, 0] = rho_std_avg
    form_E[0, 1] = F_bar
    form_E[1, 1] = F_std
    # Mean pair contribution (1/2 avoids double counting):
    # vacancy manuscript Eq. (4); CMS-2022 Eq. (3).
    form_E[0, 2] = np.sum(Pp_bar) * 0.5

    Pp_std_cn = np.zeros((np.shape(cn)[0], np.shape(comp)[0]))
    Pp_std_avg = np.zeros((np.shape(cn)[0], np.shape(comp)[0]))
    for k in np.arange(0, np.shape(cn)[0]):
        for j in np.arange(0, np.shape(comp)[0]):
            # Central-species/shell mean pair interaction. This is the SRO
            # extension of CMS-2022 Appendix-B Eq. (B.2).
            avg = np.sum(Pp_cn[k, j, :] * comp*(1-alpha[k,j,:]))
            Pp_std_avg[k, j] = avg * cn[k, 1]
            Pp_std_cn[k, j] = np.sum((Pp_cn[k, j, :] - avg)**2 * comp * (1 - alpha[k, j, :])) * cn[k, 1]

    Pp_std_cn2 = np.sqrt(np.sum(Pp_std_cn, axis=0))
    Pp_std_avg2 = np.sum(Pp_std_avg, axis=0)

    # Pair-interaction standard deviation:
    # First term = intra-species contribution
    # Second term = inter-species contribution, evaluated explicitly
    # over unique X-Y pairs.
    bb = np.zeros((np.shape(comp)[0], np.shape(comp)[0]))
    for i in np.arange(0, np.shape(bb)[0]):
        for j in np.arange(0, np.shape(bb)[0]):
            if i <= j:
                continue
            bb[i, j] = (comp[i] * comp[j] * np.square(Pp_std_avg2[i] - Pp_std_avg2[j]))

    # Standard deviation of the pair-interaction energy
    Pp_std = np.sqrt(np.sum(comp * np.square(Pp_std_cn2)) + np.sum(bb))

    form_E[1, 2] = Pp_std * 0.5
    # Mean per-atom cohesive/site energy: vacancy manuscript Eq. (4);
    # CMS-2022 Eq. (3).
    form_E[0, 3] = F_bar + np.sum(Pp_bar) * 0.5

    # Covariance cov(F, V/2): vacancy manuscript Eq. (9);
    # CMS-2022 Eq. (10).
    exp_FV = np.sum(comp * Fs * Pp_std_avg2 * 0.5) - np.sum(F_bar) * np.sum(Pp_bar) * 0.5
    covar2 = exp_FV
    # Total site-energy standard deviation from variance propagation:
    # vacancy manuscript Eq. (9); CMS-2022 Eq. (9).
    form_E[1, 3] = np.sqrt(F_std**2 + (Pp_std**2 / 4) + 2 * covar2)

    form_E = pd.DataFrame(form_E, columns=["rho", "F", "Pp", "E"])
    form_E.index = ["Mean", "Std"]

    covars = np.array([covar2])
    test = [Fs + Pp_std_avg2 * 0.5]

    return form_E, test, covars
