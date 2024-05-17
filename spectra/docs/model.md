## **Modeling the molecular rotational transitions of gas with Gaussian velocity profiles**

**By Yue Cao**

**Created Sep 23, 2020**

**Modified May 13, 2024**

### **1. Introduction**

​    This essay describes how to model the molecular rotational spectra of a gas tracer given its basic physical parameters (column density, excitation temperature, velocity dispersion, etc.). Most of the theoretical framework in this essay can be referred to [Mangum & Shirley (2015)](https://iopscience.iop.org/article/10.1086/680323/meta) (hereafter [MS15](https://iopscience.iop.org/article/10.1086/680323/meta)). My Python package `spectra` is dedicated to model spectra data and derive the relevant physical quantities based on the models described in this essay.

### **2. Symbols**

​    The symbols used throughout the essay are listed below.

| Symbol          | Type         | Description                                                  |
| --------------- | ------------ | ------------------------------------------------------------ |
| $J$             | Transitional | Total angular momentum quantum number.                       |
| $K$             | Transitional | Projected angular momentum quantum number.                   |
| $E_{\rm u}$     | Transitional | Upper energy level.                                          |
| $\nu_0$         | Transitional | Rest frequency of the transition.                            |
| $\delta v_j$    | Transitional | Velocity shift of the $j$th hyperfine line relative to $\nu_0$. Radio velocity convention used. |
| $R_j$           | Transitional | Relative strength of the $j$th hyperfine line s.t. $\Sigma_jR_j=1$. |
| $S$             | Transitional | Line strength.                                               |
| $g_J, g_K, g_I$ | Transitional | Degeneracies related to the $J$, $K$, $I$ (nuclear spin) quantum numbers. |
| $g_{\rm u}$     | Transitional | Degeneracy of the upper energy level, $=g_{J_{\rm u}}g_{K_{\rm u}}g_{I_{\rm u}}$. |
| $Q_{\rm rot}$   | Transitional | Rotational partition function.                               |
| $\mu$           | Molecular    | Electric dipole moment of a molecule.                        |
| $A_0,B_0,C_0$   | Molecular    | Rigid rotor constants of a molecule.                         |
| $\sigma$        | Molecular    | Symmetry number of a molecule.                               |
| $N$             | Macroscopic  | Total column density of a gas tracer.                        |
| $T_{\rm ex}$    | Macroscopic  | Excitation temperature of a gas tracer.                      |
| $v_0$           | Macroscopic  | Central line-of-sight velocity of a Gaussian velocity component. |
| $\sigma_v$      | Macroscopic  | 1-$\sigma$ velocity dispersion of a Gaussian velocity component. |
| $T_{\rm bg}$    | Macroscopic  | Background brightness temperature.                           |

### **3. Assumptions**

​    The following assumptions and restrictions are made for deriving the rotational transition model.

1. Only linear or symmetric top molecules are considered in this essay.
2. There are no interactions between the different kinds of quantum energy states of the molecule, and the molecule is in ground states of electronic and vibrational energies. In this way, $Q_{\rm e}=Q_{\rm v}=1$ and $Q_{\rm rot}=Q_{\rm r}Q_{\rm n}$.

3. Rigid rotor approximation.

4. McDowell's approximation for calculating $Q_{\rm rot}$ (see below).

5. For symmetric top molecules, we are only interested in one nuclear spin state (e.g. para-$\rm NH_3$), i.e. $g_I=1$.

6. The gas materials along one line-of-sight (los) can be decomposed into multiple velocity components with Gaussian velocity profile
   $$
   f_i(v)=\frac{1}{\sigma_{vi}\sqrt{2\pi}}\exp\left(-\frac{(v-v_{0i})^2}{2\sigma_{vi}^2}\right)
   $$

7. All the velocity components share the same $T_{\rm ex}$ (Otherwise the information of the 3D gas distribution is needed to integrate the radiative transfer equation).

### **4. Theory**

​    In this section we derive the modeled spectral intensity as a function of los velocity, i.e. $\Delta I_\nu(v)$, given the parameters {$T_{\rm ex}$, $N_i$, $v_{0i}$, $\sigma_{vi}$} of the velocity components and the constants listed in Sect 2. 

​    According to Eq. 32 of [MS15](https://iopscience.iop.org/article/10.1086/680323/meta) (see their derivations), the opacity of the $j$th hyperfine line of the $i$th velocity component as a function of $v$​ can be written as
$$
\tau_{vij}(v)=\frac{4\sqrt2\pi^\frac{5}{2}S\mu^2R_jg_{\rm u}}{3h\sigma_{vi}Q_{\rm rot}(T_{\rm ex})}N_i\exp\left(-\frac{(v-v_{0i}-\delta v_j)^2}{2\sigma_{vi}^2}-\frac{E_{\rm u}}{k_{\rm B}T_{\rm ex}}\right)\left[\exp\left(\frac{h\nu_0}{k_{\rm B}T_{\rm ex}}\right)-1\right]\label{eq:tau_ij}
$$
​    The total opacity
$$
\tau(v)=\sum_{ij}\tau_{\nu ij}(v)\label{eq:tau}
$$
​    The intensity
$$
\Delta I_\nu(v;T_{\rm ex},N_i,v_{0i},\sigma_{vi},T_{\rm bg})=[B_\nu(T_{\rm ex})-B_\nu(T_{\rm bg})]\left(1-e^{-\tau_\nu(v)}\right)\label{eq:intensity}
$$
​    Eq $\ref{eq:tau_ij}$, $\ref{eq:tau}$, $\ref{eq:intensity}$ are essentially the spectral model. Now we calculate the constants in Eq $\ref{eq:tau_ij}$​ for **linear and symmetric top** molecules. 

#### **4.1 Line strength $S$** 

1. For $(J,K)\rightarrow(J-1,K)$ transitions
   $$
   S=\frac{J^2-K^2}{J(2J+1)}
   $$

2. For $(J,K)\rightarrow(J,K)$​ transitions
   $$
   S=\frac{K^2}{J(J+1)}
   $$

3. For $(J,K)\rightarrow(J+1,K)$​ transitions
   $$
   S=\frac{(J+1)^2-K^2}{(J+1)(2J+1)}
   $$

​    For linear molecules, set $K=0$ in the above equations.

#### **4.2 Degeneracy $g_{\rm u}$ **

​    Degeneracy of the upper energy state $g_{\rm u}=g_{J_{\rm u}}g_{K_{\rm u}}g_{I_{\rm u}}$. For linear and symmetric top molecules, $g_{J_{\rm u}}=2J_{\rm u}+1$. For $K=0$ or linear molecules, $g_{K_{\rm u}}=1$; for $K\neq0$ symmetric top molecules, $g_{K_{\rm u}}=2$. If you are only interested in one nuclear spin state (e.g. para-$\rm NH_3$), $g_{I_{\rm u}}=1$.

#### **4.3 Rotational partition function $Q_{\rm rot}$**

​    With assumption 2 in Sect. 3,
$$
Q_{\rm rot}=Q_{\rm r}Q_{\rm n}=\sum_{J,K,I}g_Jg_Kg_I\exp\left(-\frac{E_{JK}}{k_{\rm B}T_{\rm ex}}\right)
$$
​    Here we use McDowell's approximation (see the references in [MS15](https://iopscience.iop.org/article/10.1086/680323/meta)) to evaluate $Q_{\rm rot}$​. For linear molecules,
$$
Q_{\rm rot,linear}\approx\frac{k_{\rm B}T_{\rm ex}}{hB_0}\exp\left(\frac{hB_0}{3k_{\rm B}T_{\rm ex}}\right)\label{eq:Q_rot_lin}
$$
​    This equation is more accurate at higher temperatures. For $\frac{hB_0}{k_{\rm B}T_{\rm ex}}<0.2$, the error is less than 0.01%.

​    For symmetric top molecules,
$$
Q_{\rm rot,sym}\approx\frac{\sqrt{m\pi}}{\sigma}\exp\left(\frac{hB_0(4-m)}{12k_{\rm B}T_{\rm ex}}\right)\left(\frac{k_{\rm B}T_{\rm ex}}{hB_0}\right)^\frac{3}2\left[1+\frac{1}{90}\left(\frac{hB_0(1-m)}{k_{\rm B}T_{\rm ex}}\right)^2\right]\label{eq:Q_rot_sym}
$$
where $m=B_0/A_0$ for prolate molecules and $m=B_0/C_0$ for oblate molecules. For the calculation of $\sigma$, see Sect 2 of [this paper](https://link.springer.com/article/10.1007/s00214-007-0328-0).

### **5. Numerical implementation**

​    This section describes how the spectral model (essentially Eq. $\ref{eq:tau_ij}$, $\ref{eq:intensity}$, $\ref{eq:Q_rot_lin}$, $\ref{eq:Q_rot_sym}$​) is numerically evaluated in `spectra`. We introduce the following dimensionless constants, which can be evaluated first.
$$
\begin{flalign}
C_N&=\frac{4\sqrt2\pi^\frac{5}{2}S\mu^2R_jg_{\rm u}}{3h}\rm\frac{cm^{-2}}{km\ s^{-1}}\\
C_{T1}&=\frac{E_{\rm u}}{k_{\rm B}}\rm\frac{1}{K}\\
C_{T2}&=\frac{h\nu_0}{k_{\rm B}}\rm\frac{1}{K}\\
C_I&=\frac{2h\nu_0^3}{c^2}\rm\frac{1}{Jy\ sr^{-1}}\\
C_{Q1}&=\frac{hB_0}{k_{\rm B}}\rm\frac{1}{K}\\
C_{Q2}&=\frac{\sqrt{m\pi}}{\sigma}\left(\frac{k_{\rm B}{\rm K}}{hB_0}\right)^\frac{3}2\\
C_{Q3}&=\frac{hB_0(4-m)}{12k_{\rm B}}\rm\frac{1}{K}\\
C_{Q4}&=\frac{1}{90}\left(\frac{hB_0(1-m)}{k_{\rm B}{\rm K}}\right)^2
\end{flalign}
$$
​    Eq. $\ref{eq:tau_ij}$, $\ref{eq:intensity}$, $\ref{eq:Q_rot_lin}$, $\ref{eq:Q_rot_sym}$​ can be expressed numerically as followings, respectively.
$$
\begin{flalign}
\tau_{\nu ij}(v)&=C_N\left(\frac{N_i}{\rm cm^{-2}}\right)\left(\frac{\sigma_{vi}}{\rm km\ s^{-1}}\right)^{-1}\frac{1}{Q_{\rm rot}(T_{\rm ex})}e^{-\frac{(v-v_{0i}-\delta v_j)^2}{2\sigma_{vi}^2}-C_{T1}\left(\frac{T_{\rm ex}}{\rm K}\right)^{-1}}\left(e^{C_{T2}\left(\frac{T_{\rm ex}}{\rm K}\right)^{-1}}-1\right)\\
\frac{\Delta I_\nu(v)}{\rm Jy\ sr^{-1}}&=C_I\left[\left(e^{C_{T2}\left(\frac{T_{\rm ex}}{\rm K}\right)^{-1}}-1\right)^{-1}-\left(e^{C_{T2}\left(\frac{T_{\rm bg}}{\rm K}\right)^{-1}}-1\right)^{-1}\right]\left(1-e^{-\tau_\nu(v)}\right)\label{eq:intensity_num}\\
Q_{\rm rot,lin}&=\frac{1}{C_{Q1}}\exp\left(\frac{C_{Q1}}3\left(\frac{T_{\rm ex}}{\rm K}\right)^{-1}\right)\\
Q_{\rm rot,sym}&=C_{Q2}\left(\frac{T_{\rm ex}}{\rm K}\right)^\frac{3}2\exp\left(C_{Q3}\left(\frac{T_{\rm ex}}{\rm K}\right)^{-1}\right)\left[1+C_{Q4}\left(\frac{T_{\rm ex}}{\rm K}\right)^{-2}\right]
\end{flalign}
$$
​    We can also calculate the brightness temperature of the spectrum 
$$
T_{\rm b}=C_{T2}\ln^{-1}\left[1+\left(\left(e^{C_{T2}\left(\frac{T_{\rm ex}}{\rm K}\right)^{-1}}-1\right)^{-1}-\left(e^{C_{T2}\left(\frac{T_{\rm bg}}{\rm K}\right)^{-1}}-1\right)^{-1}\right)^{-1}\left(1-e^{-\tau_\nu}\right)^{-1}\right]\label{eq:T_b_num}
$$
​    For the optically thin limit, replace $(1-e^{-\tau_\nu})$ with $\tau_\nu$ in Eq. $\ref{eq:intensity}$, $\ref{eq:intensity_num}$, $\ref{eq:T_b_num}$​. 

​    The mean opacity of a velocity component, $\bar\tau_i$​, can be given as
$$
\begin{flalign*}
\bar\tau_i&=\frac{1}{2\sqrt{2\ln2}\sigma_{vi}}\int\sum_j\tau_{\nu ij}(v){\rm d}v\\
&=\frac{\sqrt\pi}{2\sqrt{\ln2}}\left(\sum_jC_{Nj}\right)\left(\frac{N_i}{\rm cm^{-2}}\right)\left(\frac{\sigma_{vi}}{\rm km\ s^{-1}}\right)^{-1}\frac{1}{Q_{\rm rot}(T_{\rm ex})}e^{-C_{T1}\left(\frac{T_{\rm ex}}{\rm K}\right)^{-1}}\left(e^{C_{T2}\left(\frac{T_{\rm ex}}{\rm K}\right)^{-1}}-1\right) 
\end{flalign*}\label{eq:mean_tau}
$$
​    In this way Eq $\ref{eq:tau}$​ can be written as
$$
\tau_\nu(v)=2\sqrt{\frac{\ln2}\pi}\sum_{ij}R_j\bar\tau_i\exp\left(-\frac{(v-v_{0i}-\delta v_j)^2}{2\sigma_{vi}^2}\right)
$$
​    If you want to calculate the mean opacity within a velocity range $\Delta v$ (e.g. a velocity channel), replace $\sigma_v$ in Eq $\ref{eq:mean_tau}$ with $\Delta v/(2\sqrt{2\ln2})$.

### **6. Flowchart of** `spectra`

<img src="/Users/yuecao/Documents/astro/spectra/model.assets/Screenshot 2024-05-14 at 12.22.02.png" alt="Screenshot 2024-05-14 at 12.22.02" style="zoom:50%;" />

### **7. Some discussions**

​    Here we discuss the behavior of the peak intensity $\Delta I_{\nu,\rm p}$ as a function of $N$ and $T_{\rm ex}$. For simplicity we make the following restrictions.

1. Only consider linear molecules.
2. One velocity component and one hyperfine line. 
3. $T_{\rm bg}=0$.

of the modeled spectrum as a function of

and .

In this way is similar to the function

, which saturates at high . is similar to the

In this way function close to

is similar to the function , which saturates at high . is similar to the . When (high ), it is close to ; when (low ), it is

. These functions are demonstrated below.



































