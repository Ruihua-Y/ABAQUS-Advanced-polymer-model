# ABAQUS-Advanced-polymer-model
# The source code may have some potential issues and is for communication and learning purposes only
1. All files(umat/vumat) were written and compiled by Ruihua Yin. If you are willing to collaborate with us, we will share the source code
2. BB_UMAT is the Bergström–Boyce dual-network finite viscoelastic model
3. BB_UMAT2 is the BB model (F·Fv decomposition)
4. BB_Mullins_UMAT is the BB model + Ogden–Roxburgh Mullins damage
5. BB_Mullins_UMAT2 is the BB model (F·Fv decomposition) + Mullins damage
6. BB_UMAT4 is the thermo-mechanical model for double-network glassy polymers
7. HCM70_UMAT_VUMAT is a platform for 70 types of hyperelastic constitutive models
8. The compilation environment is: ABAQUS (2021); Intel Visual Fortran 19.1.0; Visual Studio 2019
9. Before compiling the HCM70 file containing VUMAT, you must prepare ‘vaba_param.inc’  and rename it;
10. if the Explicit analysis uses `double=explicit`, please switch to `vaba_param_dp.inc` and recompile.
11. BB_UMAT material parameters：
- props(1)  Bulk modulus κ
- props(2)  A-network (equilibrium) shear modulus μ_eq
- props(3)  B-network (viscous) shear modulus μ_v
- props(4)  Number of segments in A-network (locking parameter N_eq)
- props(5)  Number of network chain segments B, N_v
- props(6)  Reference flow rate coefficient γ_ref
- props(7)  Creep (chain elongation) exponent C 
- props(8)  Stress exponent m
- props(9)  Chain elongation regularization coefficient ε_reg
- props(10) Newton convergence tolerance tol
- props(11)  Newton maximum iteration count maxit
- Status variables (*DEPVAR = 6)
- SDV(1–6)  Inverse of the right Cauchy-Green tensor of the viscous intermediate configuration C_v⁻¹
12. BB_Mullins_UMAT material parameters：
- props(1)  Bulk modulus κ
- props(2)  A-network (equilibrium) shear modulus μ_eq
- props(3)  B-network (viscous) shear modulus μ_v
- props(4)  Number of segments in A-network (locking parameter N_eq)
- props(5)  Number of network chain segments B, N_v
- props(6)  Reference flow rate coefficient γ_ref
- props(7)  Creep (chain elongation) exponent C 
- props(8)  Stress exponent m
- props(9)  Regularization coefficient ε_reg
- props(10)  tol Newton tolerance tol 
- props(11)  Newton maximum iterations maxit
- props(12)  Damage parameter r 
- props(13)  Damage energy scale û 
- props(14)  Damage parameter β 
- props(15)  Minimum recovery factor a_min
- props(16)  Recovery exponent b
- Status variables (*DEPVAR = 10)
- SDV(1–6)  Inverse of the right Cauchy-Green tensor of the viscous intermediate configuration C_v⁻¹
- SDV(7)  U_max (maximum energy in damage history)
- SDV(8)  U_min (minimum energy in damage history)
- SDV(9)  η (current damage reduction factor, 1 = no damage) 
- SDV(10)  α_d (recovery factor) 
14. BB_UMAT2 material parameters：
- props(1)  μ  Initial shear modulus (Network A) 
- props(2)  λ_L Locking stretch 
- props(3)  κ  Bulk modulus
- props(4)  s  Stiffness ratio of Network B
- props(5)  ξ (xi) Flow regularization paramete
- props(6)  C (cexp) Creep exponent
- props(7)  τ_base Reference flow stress
- props(8)  m (mexp) Stress exponent
- props(9)  τ_cut Stress threshold coefficient
- props(10)  ε_reg Regularization coefficient
- props(11) tol Newton tolerance
- props(12)  maxit Newton maximum iterations
- Status variables (*DEPVAR = 13)
- SDV(1)  Current total time time(2)
- SDV(2)  Cumulative equivalent viscous flow Δγ
- SDV(3)  Chain strain ½ln(tr(b)/3)
- SDV(4)  Reserved
- SDV(5–13)  Viscous deformation gradient Fᵛ
15. BB_Mullins_UMAT2 material parameters：
- props(1)  μ  Initial shear modulus (Network A) 
- props(2)  λ_L Locking stretch 
- props(3)  κ  Bulk modulus
- props(4)  s  Stiffness ratio of Network B
- props(5)  ξ (xi) Flow regularization paramete
- props(6)  C (cexp) Creep exponent
- props(7)  τ_base Reference flow stress
- props(8)  m (mexp) Stress exponent
- props(9)  τ_cut Stress threshold coefficient
- props(10)  ε_reg Regularization coefficient
- props(11) tol Newton tolerance
- props(12)  maxit Newton maximum iterations
- props(13)  Damage parameter r 
- props(14)  Damage energy scale û 
- props(15)  Damage parameter β 
- props(16)  Minimum recovery factor a_min
- props(17)  Recovery exponent b
- Status variables (*DEPVAR = 17)
- SDV(1)  Current total time time(2)
- SDV(2)  Cumulative equivalent viscous flow Δγ
- SDV(3)  Chain strain ½ln(tr(b)/3)
- SDV(4)  Reserved
- SDV(5–13)  Viscous deformation gradient Fᵛ
- SDV(14)  U_max (maximum energy in damage history)
- SDV(15)  U_min (minimum energy in damage history)
- SDV(16)  η (current damage reduction factor, 1 = no damage) 
- SDV(17)  α_d (recovery factor)
- 16. BB_UMAT4 material parameters：
- props(1)  R   Universal gas constant 
- props(2)  k_B   Boltzmann constant (used with N_S)  
- props(3)  G_B   B network shear modulus 
- props(4)  κ   B network bulk modulus  
- props(5)  μ*   Structural shear modulus at temperature T*  
- props(6)  C   Vogel constant for viscosity  
- props(7)  ΔH   Activation enthalpy  
- props(8)  Ω_s   Shear activation volume  
- props(9)  Ω_p   Pressure activation volume  
- props(10)  N_S   Conformational chain density  
- props(11)  γ0    Rate pre-factor  
- props(12)  C_s   Vogel constant for γ relaxation  
- props(13)  T*   Reference temperature 
- props(14)  T∞   Vogel limit temperature  
- props(15)  α    Conformational network non-Gaussian parameter 
- props(16)  η    Conformational network parameter 
- props(17)  a_crit   Critical chain stretch constant term 
- props(18)  b_crit   Critical chain stretch temperature coefficient
- Status variables (*DEPVAR = 18)
- SDV(1-6)  Total logarithmic strain ε
- SDV(7-12)  B-network back stress deviator τ_B^dev (Voigt)
- SDV(13-15)  Conformational network principal direction strain ε_N
- SDV(16-18)  Total principal strain from the previous increment step ε_p
17. HCM70_UMAT_VUMAT
18. In this file, props(1) is model number (1–70), selected from the model table below;
- props(2 … nprops-1) are material parameters for that model (in the order given in the model table);
- props(nprops) is bulk modulus for the nearly incompressible penalty function W_vol = ½·KBULK·(J−1)²;
- 1  3-Param Gent  α, μ, I_m
- 2  Alexander  C₁, k, C₂, γ, C₃
- 3  Amin  C₁, C₂, N, C₃, M, C₄
- 4  Arruda-Boyce  μ, N
- 5  Bahreman-Darijani  A₂, B₂, A₄, A₆
- 6  Biderman  C₁₀, C₀₁, C₂₀, C₃₀
- 7  Carroll  A, B, C
- 8  Chevalier-Marco  N(≤6), a₀…a₆, b₀…b₆ 
- 9 Davis-De-Thomas A, n, C, k
- 10 Exp-Ln A, a, b
- 11 Fung-Demiray μ, b
- 12 Gent E, I_m
- 13  Gent-Thomas  C₁, C₂
- 14  Gen-Yeoh  K₁, m, K₂, p, K₃, q
- 15  Gornet-Desmorat  h₁, h₂, h₃
- 16  Gregory  A, B, C, m, n
- 17  Haines-Wilson  C₁₀, C₀₁, C₁₁, C₀₂, C₂₀, C₃₀
- 18  Hartmann-Neff  N(≤6), α, C₁₀…C₆₀, C₀₁…C₀₆
- 19  Hart-Smith  G, k₁, k₂
- 20  Haupt-Sedlan  C₁₀, C₀₁, C₁₁, C₀₂, C₃₀
- 21  Hoss-Marczak-1  α, β, μ, b, n
- 22  Hoss-Marczak-2  α, β, μ, b, n, C₂
- 23  Humphrey-Yin  C₁, C₂
- 24  Isihara  C₁₀, C₂₀, C₀₁
- 25  James-Green-Simpson  C₁₀, C₀₁, C₁₁, C₂₀, C₃₀
- 26  Lambert-Diani-Rey  a₀, a₁, a₂, b₀, b₁
- 27  Lion  C₁₀, C₀₁, C₅₀
- 28  Li-Zhang  A, α, M, B, β, N
- 29  Lopez-Pamies  N(≤10), (μ₁,α₁)…(μ₁₀,α₁₀)
- 30  Mansouri-Darijani  A₁, m₁, B₁, n₁
- 31  Modified-Yeoh  C₁₀, C₂₀, C₃₀, α, β
- 32  Polynomial_N1 (Mooney-Rivlin)  C₁₀, C₀₁
- 33  Polynomial_N2  C₁₀, C₀₁, C₁₁, C₂₀, C₀₂
- 34  Pucci-Saccomandi  K, μ, J_m
- 35  RedPoly_N1 (Neo-Hookean)  C₁₀
- 36  RedPoly_N2  C₁₀, C₂₀
- 37  RedPoly_N3 (Yeoh)  C₁₀, C₂₀, C₃₀
- 38  RedPoly_N4  C₁₀, C₂₀, C₃₀, C₄₀
- 39  RedPoly_N5  C₁₀, C₂₀, C₃₀, C₄₀, C₅₀
- 40  RedPoly_N6  C₁₀, C₂₀, C₃₀, C₄₀, C₅₀, C₆₀
- 41  Swanson  N(≤10), (A₁,a₁)…(A₁₀,a₁₀), (B₁,b₁)…(B₁₀,b₁₀)
- 42  Takamizawa-Hayashi  c, J_m
- 43  Van der Waals  μ, λ_m, β, a
- 44  Veronda-Westmann  C₁, C₂, α
- 45  Vito  μ, b, α
- 46  Warner  μ, I_m
- 47  Yamashita-Kawabata  C₁, C₂, C₃, N
- 48  Yeoh-Fleming  A, B, I_m, C₁₀
- 49  Attard  N(≤10), A₁…A₁₀, B₁…B₁₀
- 50  Bechir-4term  C₁₁, C₂₁, C₁₂, C₂₂
- 51  Bootstrapped-8Chn  μ, N
- 52  Constr-Junctions  G_c, νkT, κ
- 53  Continuum-Hybrid  K₁, K₂, α, μ
- 54  Eight-Chain  μ, N, C₁…C₅ 
- 55  Edwards-Vilgis  G_c, G_e, α, η
- 56  Extended-Tube  G_c, G_e, δ, β
- 57  Horgan-Murphy  μ, c, I_m
- 58  Micro-Sphere  μ (optional N, p, ≤0 uses default)
- 59  Ogden_N1  μ₁, α₁
- 60  Ogden_N2  μ₁, α₁, μ₂, α₂
- 61  Ogden_N3  μ₁…μ₃, α₁…α₃
- 62  Ogden_N4  μ₁…μ₄, α₁…α₄
- 63  Ogden_N5  μ₁…μ₅, α₁…α₅
- 64  Ogden_N6  μ₁…μ₆, α₁…α₆
- 65 Peng-Landel E
- 66 Shariff E, α₂…α₆ 
- 67 Slip-Link G_e, η, α
- 68 SpT G_c, G_e, N
- 69  Tube  G_c, G_e, β
- 70  Valanis-Landel  μ
19. standardU.dll is used for Abaqus/Standard, while explicitU.dll is used for Abaqus/Explicit; both are compiled from the same source file, but their names cannot be interchanged.
20. The number of *DEPVAR for each model must be ≥ the requirements in the table above (BB_UMAT2 /BB_Mullins_UMAT2 at least 13; BB_UMAT4 must be 18; HCM70 at least 1).
21. The number of CONSTANTS in *USER MATERIAL must be provided in full according to the order in the table above;BB_UMAT4 must have exactly 18.
22. The analysis must include a temperature field (e.g., *TEMPERATURE or a coupled thermal-mechanical analysis), with temperature in units of K.
23. HCM70 does not support plane stress elements; the BB series UMAT only supports 3D (ndi=3/nshr=3).
24. All viscoelastic/damage models use implicit Newton iteration, and the time increment is automatically reduced when convergence fails; if convergence still fails, optional numerical parameters can be adjusted.
25. It is recommended that the bulk modulus be 1000–5000 times the shear modulus 
26. explicitU.dll is compiled with the default single precision (vaba_param_sp.inc); If the job uses double=explicit, the library must be recompiled with vaba_param_dp.inc.
