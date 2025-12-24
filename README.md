# HFGN: A Python code for Hartree-Fock calculations for Graphene Nanoribbons
Hartree-Fock (HF) calculations for zigzag graphene nanoribbons (ZGNR). This is the Hartree-Fock approximation for the Mott-Anderson Hamiltonian 

$$
H = -t \sum_{\left< ij \right>, \sigma} 
\left( c^\dagger_{i, \sigma} c_{j, \sigma} + c^\dagger_{j, \sigma} c_{i, \sigma}   \right) +
U \sum_{i} n_{i, \uparrow} n_{i, \downarrow} +
\sum_{i} V_i n_{i, \sigma}.
$$


The results obtained from this code are presented in these following publications:
1) Phase diagram and crossover phases of topologically ordered graphene zigzag nanoribbons: role of localization effects
Hoang-Anh Le, In-Hwan Lee, Young-Heon Kim and S.-R. Eric Yang, J. Phys.: Condens. Matter 36 265604 (2024)
https://dx.doi.org/10.1088/1361-648X/ad38f9
2) Exploring Topological Order of Zigzag Graphene Nanoribbons: Phase Diagram and Crossover Phases, Hoang Anh Le, PhD Thesis, Korea University (2024)
http://dcollection.korea.ac.kr/srch/srchDetail/000000278801

Please cite this software if you use it for your scientific research. 

Feel free to reach me via hoanganhle at qns.science or lehoanganh1112 at gmail.com if you need assistance.

# Usage

`HF_spinspliting.py` and `HF_zgnr_doped.py` are the main codes, respectively used for half-filling and doped systems. They return the eigensystems of the system: `val_up, vec_up, val_dn, vec_dn` via the following function

```
def HF_vec(Lx, Ly):
    '''
    Return Hartree-Fock eigensystem.
    < Input >
    Lx: (int) length of ribbon.
    Ly: (int) width of ribbon.
    < Output >
    val, vec
    val: (2*Lx*Ly) column vector of eigenvalues.
    vec[:, i]: (2*Lx*Ly) column-eigenvector corresponding to val[i].
    '''
    the_system = self_consistent_solution(Lx, Ly)
    val_up, vec_up = cp.linalg.eigh(the_system[0])
    val_dn, vec_dn = cp.linalg.eigh(the_system[1])
    return val_up, vec_up, val_dn, vec_dn
```

# Examples
1. Compute topological entanglement entropy (Refer to Section 3.5 of my thesis)

2. Visualize eigenstates (Section 3.6)

3. Plot $q_A$-energy diagram (Section 3.6)

4. Plot site spins of the system (Section 3.4)

5. Visualize the geometry of zigzag graphene nanoribbon








