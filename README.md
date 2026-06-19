# `CRASHE` - Collisional Radiative modelling for AStrophysics with Homologous Expansion
This is a code for solving collisional radiative problems, with atomic data input through the adf04 standard. Eventually, this will be a python package with compiled fortran doing the numerical work.
Essentially, the equation of motion is something like

$$
n_i \left(\sum_{j\neq i} n_e q_{i\to j} +\sum_{j < i} A^{\text{eff}}_{i \to j}\right) = \sum_{j\neq i} n_j n_e q_{j\to i} + \sum_{j>i}n_j A^{\text{eff}}_{j\to i},
$$

at some electron temperature $T_e$ and electron density $n_e$. Collisional excitation/dexcitation rate coefficients $q_{i\to j}$ and radiative rates $A^{\text{eff}}_{i \to j}$. These radiative rates are either the Einstein A-coefficients, or these coefficients multiplied by a Sobolev escape factor,

$$
    \tau_{j\to i} = \frac{\lambda_{j \to i} ^3}{8\pi}A_{j \to i} t ~n_i~\Big(\frac{g_j}{g_i}-\frac{n_j}{n_i}\Big),\\
    \beta_{j\to i} = \frac{1}{\tau_{j\to i}} \left(1-e^{-\tau_{j\to i}} \right).
$$

The solution of the problem is coded in this repository. Basically, define a matrix $C_{ij}$ as

$$
C_{ij} =   \begin{cases}
      q_{i\to j}, & \text{if}\ i < j \\
      -\sum_{i>j} A_{i \to j} - \sum_{i \neq j} q_{i\to j}, &  \text{if}\ i = j \\
      q_{i\to j} + A_{i \to j}, & \text{if}\ i > j \\
    \end{cases}
$$

And the steady state is obtained as the matrix equation $\sum_j C_{ij} n_j = 0$. For $N$ levels in the system, we have $N$ equations that are linearly dependent. Therefore we additionally need to enforce some other condition to have a unique solution. The particle conservation $\sum_j n_j =1$ is adequate. This is typically used in conjunction with the diagonalization of $C_{ij}$, and the selection of the eigenvector with eigenvalue zero. 

For an inversion - we choose to write all levels relative to $n_1$, e.g 

$$
n_i/n_1 \to n_i, i>1. 
$$

We therefore have for the steady state condition,

$$
\sum_{j>1} C_{ij} n_j = -C_{i1}.
$$

Define the matrix $A$ = $( C_{ij} ;~i,j>1)$, and the column vectors $x =  (n_j; ~j>2)$, $b =  (-C_{i1}; ~i>2)$. The resulting $N-1 \times N-1$ matrix equation is then solved by

$$
x=A^{-1} b,
$$

or

$$
n_i/n_1 = -\sum_{j>2 }(A^{-1})_{ij} C_{j1}
$$

which is solved by the lapack routines dgesv or dgesvx. The matrix $A$ is stored in the $N-1 \times N-1$ fortran array `crm`, and the vector of mathematical size $N-1$ $C_{i1}$ is stored in the fortran array `col1` - that is actually allocated the size $N$. This is because the resultant levels are re-normalized such that

$$
1 = \sum_i n_i,
$$

and all of the populations including the ground (total $N$) are stored in `col1`. See the routine `solve_cr_populations_axb` (where these variables are masked as `Q` and `qcol1` respectively).

Should the user request Sobolev opacity with the namelist variable `sobolev=.true.`, they are required to additionally specify the time since explosion in days `timeSinceExplosionDays`. They should also specify one of:
- `velocityExpansionC` - the rough expansion velocity in units of the speed of light, as well as `massElementSolar` the ionic mass in solar masses. From these two variables, combined with `timeSinceExplosionDays` - an estimate of the absolute number density of the element is calculated. Provided the code has been given an atomic number (which is gauranteed in the adf04), it makes an estimate of the number of particles as `massElementSolar / A(Z)` where `A(Z)` is the atomic mass.
- or, `fractionOverride`- which sets the total ionic number density as `fractionOverride * electronDensity`. Note for total luminosity calculations the code will still need velocity (broadening) and mass (Sets the magnitude). 