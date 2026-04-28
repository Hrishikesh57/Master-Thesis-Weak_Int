# Weak Energy Corrections Basis-Sum Runlog

- Generated at: 2026-04-21 09:31:07 UTC
- Repository commit: `0da3419`
- Compiler:
```text
GNU Fortran (GCC) 15.2.1 20260123 (Red Hat 15.2.1-7)
Copyright (C) 2025 Free Software Foundation, Inc.
This is free software; see the source for copying conditions.  There is NO
warranty; not even for MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
```

## Build And Run Commands

```text
(cd /home/hrishikeshsriman-narayan/Documents/Masters/Master-Thesis--Weak_Int && git rev-parse --short HEAD)
```
- returncode = 0
- stdout tail:
```text
0da3419
```

```text
gfortran --version
```
- returncode = 0
- stdout tail:
```text
GNU Fortran (GCC) 15.2.1 20260123 (Red Hat 15.2.1-7)
Copyright (C) 2025 Free Software Foundation, Inc.
This is free software; see the source for copying conditions.  There is NO
warranty; not even for MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
```

```text
make -C /home/hrishikeshsriman-narayan/Documents/Masters/Master-Thesis--Weak_Int/workspace_f2py/programs/rad practise something
```
- returncode = 0
- stdout tail:
```text
make: Entering directory '/home/hrishikeshsriman-narayan/Documents/Masters/Master-Thesis--Weak_Int/workspace_f2py/programs/rad'
make: Nothing to be done for 'practise'.
make: Nothing to be done for 'something'.
make: Leaving directory '/home/hrishikeshsriman-narayan/Documents/Masters/Master-Thesis--Weak_Int/workspace_f2py/programs/rad'
```

```text
(cd /home/hrishikeshsriman-narayan/Documents/Masters/Master-Thesis--Weak_Int/workspace_f2py && bash /home/hrishikeshsriman-narayan/Documents/Masters/Master-Thesis--Weak_Int/workspace_f2py/programs/rad_f2py/build_f2py.sh)
```
- returncode = 0
- stdout tail:
```text
              e_pot,r_rel,ierr = first_order_weak_shift(znucl,nneutrons,r_fm,r,g,f)
		Creating wrapper for Fortran subroutine "second_order_shift_sampled"("second_order_shift_sampled")...
            Constructing wrapper function "rad_kernel_f2py.second_order_shift_sampled"...
              m_ns,delta_e,r_rel,ierr = second_order_shift_sampled(znucl,nneutrons,r_fm,e_initial,e_intermediate,r_intermediate,g_intermediate,f_intermediate,r_initial,g_initial,f_initial)
		Creating wrapper for Fortran subroutine "second_order_basis_sum_rel"("second_order_basis_sum_rel")...
            Constructing wrapper function "rad_kernel_f2py.second_order_basis_sum_rel"...
              delta_e_total,num_positive_available,num_states_used,excluded_basis_index,energies_rel,matrix_elements,contributions,r_rel,ierr = second_order_basis_sum_rel(znucl,nneutrons,r_fm,ref_nmain,ref_kappa,e_initial,basis_kappa,basis_order,basis_rcavity,max_positive_states,knots,spectrum,coefmat,r_initial,g_initial,f_initial)
		Creating wrapper for Fortran subroutine "discrete_amp_sampled"("discrete_amp_sampled")...
            Constructing wrapper function "rad_kernel_f2py.discrete_amp_sampled"...
              raw_a_s,raw_a_w,a_s_final,a_w_final,eta,r_fk,r_ki,r_sk_weak,r_kf_weak,r_rel,ierr = discrete_amp_sampled(znucl,nneutrons,r_fm,e_initial,e_final,e_intermediate,e0_field,state_factor,weak_enabled,r_initial,g_initial,f_initial,r_intermediate,g_intermediate,f_intermediate,r_final,g_final,f_final)
    Wrote C/API module "rad_kernel_suite" to file "/home/hrishikeshsriman-narayan/Documents/Masters/Master-Thesis--Weak_Int/workspace_f2py/programs/rad_f2py/build/rad_kernel_suitemodule.c"
    Fortran 90 wrappers are saved to "/home/hrishikeshsriman-narayan/Documents/Masters/Master-Thesis--Weak_Int/workspace_f2py/programs/rad_f2py/build/rad_kernel_suite-f2pywrappers2.f90"
```
- stderr tail:
```text
f951: Warning: ‘&’ not allowed by itself in line 206
f951: Warning: ‘&’ not allowed by itself in line 248
```

```text
(cd /home/hrishikeshsriman-narayan/Documents/Masters/Master-Thesis--Weak_Int/workspace_f2py/programs/rad && /home/hrishikeshsriman-narayan/Documents/Masters/Master-Thesis--Weak_Int/workspace_f2py/programs/rad/bin/practise)
stdin:
```text
1
1
55
1s
4.8041000000000000
```
```
- returncode = 0
- stdout tail:
```text
  197.797798   -8.8228877843E-015   -2.5861738049E-015
  197.997998   -8.8920467320E-015   -2.6119985034E-015
  198.198198   -9.0135519991E-015   -2.5803234391E-015
  198.398398   -9.1959177023E-015   -2.4807041974E-015
  198.598599   -9.4397801268E-015   -2.3075871074E-015
  198.798799   -9.7382706438E-015   -2.0612835751E-015
  198.998999   -1.0078098831E-014   -1.7484944000E-015
  199.199199   -1.0441259934E-014   -1.3823060089E-015
  199.399399   -1.0807228854E-014   -9.8159718916E-016
  199.599600   -1.1155446699E-014   -5.6983285667E-016
  199.799800   -1.1467847678E-014   -1.7333255238E-016
  200.000000   -1.1731162335E-014    1.8079814318E-016
```

```text
(cd /home/hrishikeshsriman-narayan/Documents/Masters/Master-Thesis--Weak_Int/workspace_f2py/programs/rad && /home/hrishikeshsriman-narayan/Documents/Masters/Master-Thesis--Weak_Int/workspace_f2py/programs/rad/bin/practise)
stdin:
```text
1
2
55
p
4.8041000000000000
```
```
- returncode = 0
- stdout tail:
```text
 2 - Basis data
 Enter your choice (1 or 2):
 Single-electron basis-data export selected.
 Enter the value for znucl (atomic number Z):
 Enter the channel (choose among: s, p, p3, d, d5, f, f7):
 Enter the value for rms radius (in fm):
 nspectrum =         2010
 File opened successfully: outputs/wavefunctions/single_electron/55_spectrum_p.txt
 # ======================================================================================
 # Single-electron basis-data export completed successfully!
# Basis data written to: outputs/wavefunctions/single_electron/55_spectrum_p.txt
 # ======================================================================================
```

```text
(cd /home/hrishikeshsriman-narayan/Documents/Masters/Master-Thesis--Weak_Int/workspace_f2py/programs/rad && /home/hrishikeshsriman-narayan/Documents/Masters/Master-Thesis--Weak_Int/workspace_f2py/programs/rad/bin/practise)
stdin:
```text
1
2
55
p3
4.8041000000000000
```
```
- returncode = 0
- stdout tail:
```text
 2 - Basis data
 Enter your choice (1 or 2):
 Single-electron basis-data export selected.
 Enter the value for znucl (atomic number Z):
 Enter the channel (choose among: s, p, p3, d, d5, f, f7):
 Enter the value for rms radius (in fm):
 nspectrum =         2010
 File opened successfully: outputs/wavefunctions/single_electron/55_spectrum_p3.txt
 # ======================================================================================
 # Single-electron basis-data export completed successfully!
# Basis data written to: outputs/wavefunctions/single_electron/55_spectrum_p3.txt
 # ======================================================================================
```

```text
(cd /home/hrishikeshsriman-narayan/Documents/Masters/Master-Thesis--Weak_Int/workspace_f2py/programs/rad && /home/hrishikeshsriman-narayan/Documents/Masters/Master-Thesis--Weak_Int/workspace_f2py/programs/rad/bin/practise)
stdin:
```text
1
1
55
2s
4.8041000000000000
```
```
- returncode = 0
- stdout tail:
```text
  197.797798   -1.4907614290E-014   -3.2698913065E-015
  197.997998   -1.5852066540E-014   -3.2729911549E-015
  198.198198   -1.6800753637E-014   -3.2178299803E-015
  198.398398   -1.7733755524E-014   -3.0981644421E-015
  198.598599   -1.8628869473E-014   -2.9102265756E-015
  198.798799   -1.9462463646E-014   -2.6532947634E-015
  198.998999   -2.0210568377E-014   -2.3301036150E-015
  199.199199   -2.0850144652E-014   -1.9470449451E-015
  199.399399   -2.1360447901E-014   -1.5141141627E-015
  199.599600   -2.1724384398E-014   -1.0445681912E-015
  199.799800   -2.1929744700E-014   -5.5430765719E-016
  200.000000   -2.1970206636E-014   -6.1043256578E-017
```

```text
(cd /home/hrishikeshsriman-narayan/Documents/Masters/Master-Thesis--Weak_Int/workspace_f2py/programs/rad && /home/hrishikeshsriman-narayan/Documents/Masters/Master-Thesis--Weak_Int/workspace_f2py/programs/rad/bin/practise)
stdin:
```text
1
1
87
1s
5.5914999999999999
```
```
- returncode = 0
- stdout tail:
```text
  197.797798   -9.4805748793E-015    1.2224000006E-016
  197.997998   -7.9713808816E-015    1.5397095055E-015
  198.198198   -6.0835063947E-015    2.8646263911E-015
  198.398398   -3.9601806943E-015    3.9930653526E-015
  198.598599   -1.7647413879E-015    4.8379756852E-015
  198.798799    3.3292814399E-016    5.3367288588E-015
  198.998999    2.1704794190E-015    5.4563875882E-015
  199.199199    3.6063993357E-015    5.1960893965E-015
  199.399399    4.5314178693E-015    4.5865431690E-015
  199.599600    4.8771846994E-015    3.6869493947E-015
  199.799800    4.6216280017E-015    2.5797390366E-015
  200.000000    3.7906560550E-015    1.3635885864E-015
```

```text
(cd /home/hrishikeshsriman-narayan/Documents/Masters/Master-Thesis--Weak_Int/workspace_f2py/programs/rad && /home/hrishikeshsriman-narayan/Documents/Masters/Master-Thesis--Weak_Int/workspace_f2py/programs/rad/bin/practise)
stdin:
```text
1
2
87
p
5.5914999999999999
```
```
- returncode = 0
- stdout tail:
```text
 2 - Basis data
 Enter your choice (1 or 2):
 Single-electron basis-data export selected.
 Enter the value for znucl (atomic number Z):
 Enter the channel (choose among: s, p, p3, d, d5, f, f7):
 Enter the value for rms radius (in fm):
 nspectrum =         2010
 File opened successfully: outputs/wavefunctions/single_electron/87_spectrum_p.txt
 # ======================================================================================
 # Single-electron basis-data export completed successfully!
# Basis data written to: outputs/wavefunctions/single_electron/87_spectrum_p.txt
 # ======================================================================================
```

```text
(cd /home/hrishikeshsriman-narayan/Documents/Masters/Master-Thesis--Weak_Int/workspace_f2py/programs/rad && /home/hrishikeshsriman-narayan/Documents/Masters/Master-Thesis--Weak_Int/workspace_f2py/programs/rad/bin/practise)
stdin:
```text
1
2
87
p3
5.5914999999999999
```
```
- returncode = 0
- stdout tail:
```text
 2 - Basis data
 Enter your choice (1 or 2):
 Single-electron basis-data export selected.
 Enter the value for znucl (atomic number Z):
 Enter the channel (choose among: s, p, p3, d, d5, f, f7):
 Enter the value for rms radius (in fm):
 nspectrum =         2010
 File opened successfully: outputs/wavefunctions/single_electron/87_spectrum_p3.txt
 # ======================================================================================
 # Single-electron basis-data export completed successfully!
# Basis data written to: outputs/wavefunctions/single_electron/87_spectrum_p3.txt
 # ======================================================================================
```

```text
(cd /home/hrishikeshsriman-narayan/Documents/Masters/Master-Thesis--Weak_Int/workspace_f2py/programs/rad && /home/hrishikeshsriman-narayan/Documents/Masters/Master-Thesis--Weak_Int/workspace_f2py/programs/rad/bin/practise)
stdin:
```text
1
1
87
2s
5.5914999999999999
```
```
- returncode = 0
- stdout tail:
```text
  197.797798   -1.4891638950E-014   -5.0509009170E-016
  197.997998   -1.3738585185E-014    8.4521835220E-017
  198.198198   -1.2480380079E-014    7.0148575503E-016
  198.398398   -1.1194517984E-014    1.2942268411E-015
  198.598599   -9.9629393964E-015    1.8145036142E-015
  198.798799   -8.8653328040E-015    2.2216732648E-015
  198.998999   -7.9725344765E-015    2.4861353212E-015
  199.199199   -7.3406379451E-015    2.5915783266E-015
  199.399399   -7.0062856383E-015    2.5359243562E-015
  199.599600   -6.9834666444E-015    2.3310213263E-015
  199.799800   -7.2620276078E-015    2.0011899730E-015
  200.000000   -7.8079906036E-015    1.5807797858E-015
```

```text
(cd /home/hrishikeshsriman-narayan/Documents/Masters/Master-Thesis--Weak_Int/workspace_f2py/programs/rad && /home/hrishikeshsriman-narayan/Documents/Masters/Master-Thesis--Weak_Int/workspace_f2py/programs/rad/bin/practise)
stdin:
```text
1
1
92
1s
5.8571000000000000
```
```
- returncode = 0
- stdout tail:
```text
  197.797798    5.8172646237E-015   -7.0412647346E-016
  197.997998    4.6591209539E-015   -1.3069813104E-015
  198.198198    2.9263244558E-015   -2.1393489865E-015
  198.398398    7.1519968810E-016   -3.1489640085E-015
  198.598599   -1.8413624750E-015   -4.2643725340E-015
  198.798799   -4.5836580611E-015   -5.4017323873E-015
  198.998999   -7.3370311380E-015   -6.4727562341E-015
  199.199199   -9.9249708813E-015   -7.3930577281E-015
  199.399399   -1.2182373669E-014   -8.0900107588E-015
  199.599600   -1.3967777278E-014   -8.5091876779E-015
  199.799800   -1.5173554047E-014   -8.6187722239E-015
  200.000000   -1.5733399483E-014   -8.4117731935E-015
```

```text
(cd /home/hrishikeshsriman-narayan/Documents/Masters/Master-Thesis--Weak_Int/workspace_f2py/programs/rad && /home/hrishikeshsriman-narayan/Documents/Masters/Master-Thesis--Weak_Int/workspace_f2py/programs/rad/bin/practise)
stdin:
```text
1
2
92
p
5.8571000000000000
```
```
- returncode = 0
- stdout tail:
```text
 2 - Basis data
 Enter your choice (1 or 2):
 Single-electron basis-data export selected.
 Enter the value for znucl (atomic number Z):
 Enter the channel (choose among: s, p, p3, d, d5, f, f7):
 Enter the value for rms radius (in fm):
 nspectrum =         2010
 File opened successfully: outputs/wavefunctions/single_electron/92_spectrum_p.txt
 # ======================================================================================
 # Single-electron basis-data export completed successfully!
# Basis data written to: outputs/wavefunctions/single_electron/92_spectrum_p.txt
 # ======================================================================================
```

```text
(cd /home/hrishikeshsriman-narayan/Documents/Masters/Master-Thesis--Weak_Int/workspace_f2py/programs/rad && /home/hrishikeshsriman-narayan/Documents/Masters/Master-Thesis--Weak_Int/workspace_f2py/programs/rad/bin/practise)
stdin:
```text
1
2
92
p3
5.8571000000000000
```
```
- returncode = 0
- stdout tail:
```text
 2 - Basis data
 Enter your choice (1 or 2):
 Single-electron basis-data export selected.
 Enter the value for znucl (atomic number Z):
 Enter the channel (choose among: s, p, p3, d, d5, f, f7):
 Enter the value for rms radius (in fm):
 nspectrum =         2010
 File opened successfully: outputs/wavefunctions/single_electron/92_spectrum_p3.txt
 # ======================================================================================
 # Single-electron basis-data export completed successfully!
# Basis data written to: outputs/wavefunctions/single_electron/92_spectrum_p3.txt
 # ======================================================================================
```

```text
(cd /home/hrishikeshsriman-narayan/Documents/Masters/Master-Thesis--Weak_Int/workspace_f2py/programs/rad && /home/hrishikeshsriman-narayan/Documents/Masters/Master-Thesis--Weak_Int/workspace_f2py/programs/rad/bin/practise)
stdin:
```text
1
1
92
2s
5.8571000000000000
```
```
- returncode = 0
- stdout tail:
```text
  197.797798    4.5925526175E-016   -6.6751064472E-016
  197.997998   -6.5914062286E-016   -1.0371165763E-015
  198.198198   -2.0367712005E-015   -1.5077664688E-015
  198.398398   -3.6083059858E-015   -2.0400316725E-015
  198.598599   -5.2907350489E-015   -2.5865099922E-015
  198.798799   -6.9890386916E-015   -3.0961282712E-015
  198.998999   -8.6029271906E-015   -3.5188876157E-015
  199.199199   -1.0034166002E-014   -3.8106277041E-015
  199.399399   -1.1193892427E-014   -3.9373115375E-015
  199.599600   -1.2009262106E-014   -3.8783165402E-015
  199.799800   -1.2428861439E-014   -3.6284131382E-015
  200.000000   -1.2426516616E-014   -3.1983606453E-015
```

## Basis Parameters

- Single-electron generation path: `workspace_f2py/programs/rad/src/practise.f90`
- Reference wavefunction mode: `practise` single-electron sampled export (`wavefunction_type=1`, `output_data_type=1`)
- Basis-spectrum mode: `practise` single-electron basis export (`wavefunction_type=1`, `output_data_type=2`)
- Spline order: `9`
- Number of knots: `1016`
- Number of B-splines: `1007`
- Total basis eigenstates per channel: `2010`
- Basis cavity radius from spectrum metadata [rel]: `3.4258999784750000e+02`
- Knot-generation convention from source: `kntseq1(rnucl, rcavity, npoints/10, order, knots)` with `npoints=1000`, so `ninside=100`
- Weak-integration radius convention: `R_rel = rms_fm / electron_reduced_compton` and denominators use binding energies in `m_e = c = \hbar = 1` units

## Per-Case State Counts

- `Z=55`, `1s`: np1/2 used = `1005` of `1005` positive states; np3/2 used = `0` of `1005` positive states
- `Z=55`, `2s`: np1/2 used = `1005` of `1005` positive states; np3/2 used = `0` of `1005` positive states
- `Z=87`, `1s`: np1/2 used = `1005` of `1005` positive states; np3/2 used = `0` of `1005` positive states
- `Z=87`, `2s`: np1/2 used = `1005` of `1005` positive states; np3/2 used = `0` of `1005` positive states
- `Z=92`, `1s`: np1/2 used = `1005` of `1005` positive states; np3/2 used = `0` of `1005` positive states
- `Z=92`, `2s`: np1/2 used = `1005` of `1005` positive states; np3/2 used = `0` of `1005` positive states

## Validation

- Deterministic rerun for `Z=92`, `1s`: `True`
- Rerun `DeltaE2_np12_rel`: `-1.1917855324955772e-18`
- Rerun `DeltaE2_np32_rel`: `0.0000000000000000e+00`
- Rerun total correction: `-2.4649880334212528e-12`
- Discrete `2p`-only diagnostic contribution [rel]: `-3.6950304449584790e-22`
- Full `np1/2` basis-sum contribution [rel]: `-1.1917855324955772e-18`
- Truncated `np1/2` cutoff diagnostic with first `905` positive states [rel]: `-1.1719588389773858e-18`
- Largest `Z=92`, `1s`, `np1/2` intermediate-state contributions:
  1. basis_index=1515, E_n=+6.8316135556263589e+01, M_nref=+6.7384205919073839e-10, contribution=-6.6214464985534243e-21
  2. basis_index=1513, E_n=+6.6825437571951397e+01, M_nref=+6.6638081566089688e-10, contribution=-6.6195212634519458e-21
  3. basis_index=1517, E_n=+6.9839652040861566e+01, M_nref=+6.8112285389746267e-10, contribution=-6.6182700724198288e-21
  4. basis_index=1511, E_n=+6.5366850976190321e+01, M_nref=+6.5875787054277620e-10, contribution=-6.6127203145647575e-21
  5. basis_index=1519, E_n=+7.1396709698049094e+01, M_nref=+6.8820333264692854e-10, contribution=-6.6097634566902277e-21

## Warnings

- np3/2 basis file contained 1005 positive-energy pseudostates, but the existing code-level weak angular rule suppresses kappa=-2 contributions exactly
- np3/2 basis file contained 1005 positive-energy pseudostates, but the existing code-level weak angular rule suppresses kappa=-2 contributions exactly
- np3/2 basis file contained 1005 positive-energy pseudostates, but the existing code-level weak angular rule suppresses kappa=-2 contributions exactly
- np3/2 basis file contained 1005 positive-energy pseudostates, but the existing code-level weak angular rule suppresses kappa=-2 contributions exactly
- np3/2 basis file contained 1005 positive-energy pseudostates, but the existing code-level weak angular rule suppresses kappa=-2 contributions exactly
- np3/2 basis file contained 1005 positive-energy pseudostates, but the existing code-level weak angular rule suppresses kappa=-2 contributions exactly
