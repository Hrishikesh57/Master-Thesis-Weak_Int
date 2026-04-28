# Literature novelty and positioning check

Search date: 2026-04-27.

Project target: spin-independent, parity-conserving Standard-Model weak contribution to the H-like `1s` bound-electron `g` factor, evaluated through a central-potential mass-derivative relation and finite proton/neutron Fermi densities.

Internal benchmark values to position:

```text
Pb-208, H-like Pb^{81+}: delta_g_W(1s, zero neutron skin) = -1.3328870635281932e-11
U-238,  H-like U^{91+}:  delta_g_W(1s, zero neutron skin) = -3.1233277218115705e-11
```

For neutron skins from `0.00` to `0.20 fm`, the magnitude decreases by about `0.665%` in Pb and `0.777%` in U.

## 1. Search Scope And Keywords

Searches used combinations of:

```text
weak interaction correction bound electron g factor
parity conserving weak correction bound-electron g factor
electroweak correction bound electron g factor
Karshenboim Lee Milstein weak correction g factor
weak interaction contribution to g factor highly charged ions
bound-electron g factor weak interaction H-like uranium
parity nonconserving bound-electron g factor
weak energy shift mass derivative g factor bound electron
neutron skin bound-electron g factor
finite nuclear size bound-electron g factor neutron density
```

Sources checked included APS journal pages, Nature articles, ScienceDirect abstracts, OSTI metadata, arXiv-indexed results, PubMed, MDPI review/context papers, and accessible full-text thesis/review material.

## 2. Direct Prior Work On Weak Corrections To Bound-Electron g Factors

### No direct peer-reviewed Pb/U finite p/n-density weak-g benchmark found

This search did not find a peer-reviewed paper that directly tabulates the Standard-Model spin-independent weak correction to the H-like `1s` bound-electron `g` factor for Pb/U using finite proton and neutron densities, nor a paper that studies neutron-skin sensitivity of that correction.

That statement is deliberately limited: it means **no direct prior calculation was found in this search**, not that none exists.

### Shabaev, Glazov, Plunien, Volotka review, 2015

Citation: V. M. Shabaev, D. A. Glazov, G. Plunien, and A. V. Volotka, "Theory of Bound-Electron g Factor in Highly Charged Ions," *J. Phys. Chem. Ref. Data* **44**, 031205 (2015). DOI: https://doi.org/10.1063/1.4921299 . OSTI page: https://www.osti.gov/biblio/22413365 .

What they computed/reviewed: the standard theory of bound-electron `g` factors in highly charged ions, including relativistic, QED, nuclear recoil, nuclear structure, and interelectronic-interaction corrections. The review's H-like decomposition is:

```text
g = g_D + Delta g_QED + Delta g_NR + Delta g_NS + Delta g_NP
```

where the last three terms are nuclear recoil, nuclear size, and nuclear polarization. The review also notes free-electron hadronic and weak contributions in the free-QED part, but the accessible summary and text do not show a finite-nuclear-density Standard-Model weak source calculation of the kind used in this project.

Relevance: must cite for the conventional bound-electron `g` correction budget and for showing where this work fits: it is not replacing QED/recoil/nuclear-size theory, but adding a tiny Standard-Model weak central-source estimate.

### Karshenboim, Lee, Milstein, 2005

Citation: S. G. Karshenboim, R. N. Lee, and A. I. Milstein, "g factor of an electron or muon bound by an arbitrary central potential," *Phys. Rev. A* **72**, 042101 (2005). DOI: https://doi.org/10.1103/PhysRevA.72.042101 .

What they computed: a general method expressing the bound-particle `g` factor and central-potential corrections in terms of binding-energy information. Applications include finite nuclear size and vacuum-polarization-type potential corrections.

H-like / 1s: yes, the formalism is directly relevant to H-like states and arbitrary central potentials.

Pb/U finite p/n weak densities: not found in the accessible abstract/metadata. The paper provides the formula-level machinery, not the finite proton/neutron weak-density benchmark.

Relation to this project: this is the main conceptual ancestor of the mass-derivative route. The project is best framed as applying this central-potential machinery to the Standard-Model weak source with finite p/n density modeling.

### Non-peer-reviewed thesis material

Citation: H. Sriman-Narayan, "Weak Interaction Effects in Highly Charged Ions and Heavy Atoms," Master's thesis, Heidelberg/MPIK, August 1, 2025. Public PDF: https://pure.mpg.de/rest/items/item_3665358/component/file_3665359/content .

What it contains: a section titled "Weak Interaction Correction to the Bound-Electron g factor" and analytic/methodological material on weak effects in H-like systems. It is not a peer-reviewed journal calculation and appears to be project-provenance material rather than an independent external prior.

Relation to this project: this is a red flag only for internal provenance and citation hygiene. It does not remove novelty for a later peer-reviewed finite-density numerical benchmark, but the paper should not pretend the mass-derivative idea was discovered here.

## 3. Related Mass-Derivative / Virial g-Factor Formulas

The mass-derivative or virial route is not novel by itself.

Must-cite formula/theory papers:

- Karshenboim, Lee, Milstein, *Phys. Rev. A* **72**, 042101 (2005), DOI: https://doi.org/10.1103/PhysRevA.72.042101 .
- V. M. Shabaev, "Virial relations for the Dirac equation and their applications to calculations of hydrogen-like atoms," in *Precision Physics of Simple Atomic Systems*, Springer (2003). This is repeatedly cited in modern bound-electron `g` reviews as the virial-relations foundation.
- D. A. Glazov and V. M. Shabaev, "Finite nuclear size correction to the bound-electron g factor in a hydrogenlike atom," *Phys. Lett. A* **297**, 408-411 (2002), DOI: https://doi.org/10.1016/S0375-9601(02)00021-X .

Positioning consequence: the paper should not claim novelty for using a mass derivative. The novelty, if retained, is the specific Standard-Model weak-source application with finite proton/neutron densities and neutron-skin sensitivity.

## 4. Related Finite Nuclear-Size And Nuclear-Density g-Factor Work

### Finite nuclear size and charge-density modeling

Key source: I. A. Valuev, Z. Harman, C. H. Keitel, and N. S. Oreshkina, "Skyrme-type nuclear interaction as a tool for calculating the finite-nuclear-size correction to atomic energy levels and the bound-electron g factor," *Phys. Rev. A* **101**, 062502 (2020). DOI: https://doi.org/10.1103/PhysRevA.101.062502 .

What they computed: finite-nuclear-size corrections to energy levels and bound-electron `g` factors for highly charged H-like ions using realistic proton charge distributions, compared with homogeneous sphere, two-parameter Fermi, and experimental charge distributions.

H-like ions: yes.

Pb/U: likely within heavy-ion scope, but this is a charge-distribution/FNS `g` calculation, not a Standard-Model weak-source calculation.

Finite p/n weak densities or neutron skin: not the target; the focus is proton charge distribution and FNS correction.

Relation to this project: must cite as the closest nuclear-density modeling context for bound-electron `g` factors. It strengthens the case that nuclear-density modeling matters, but it does not appear to compute the weak correction studied here.

### Nuclear shape and nuclear polarization

Modern bound-electron `g` correction budgets include nuclear shape, nuclear polarization, recoil, and QED effects. Examples appearing in recent review/reference lists include:

- J. Zatorski, N. S. Oreshkina, C. H. Keitel, and Z. Harman, "Nuclear shape effect on the g factor of hydrogenlike ions," *Phys. Rev. Lett.* **108**, 063005 (2012). DOI: https://doi.org/10.1103/PhysRevLett.108.063005 .
- A. V. Nefiodov, G. Plunien, and G. Soff, "Nuclear-polarization correction to the bound-electron g factor in heavy hydrogenlike ions," *Phys. Rev. Lett.* **89**, 081802 (2002). DOI: https://doi.org/10.1103/PhysRevLett.89.081802 .
- N. Michel, J. Zatorski, N. S. Oreshkina, and C. H. Keitel, "Nonperturbative analysis of nuclear shape effects on the bound electron g factor," *Phys. Rev. A* **99**, 012505 (2019). DOI: https://doi.org/10.1103/PhysRevA.99.012505 .

Relation to this project: these are not weak-interaction calculations, but they are essential background because the weak correction is much smaller than many nuclear-structure corrections and must be framed as part of a hierarchy.

## 5. Related Precision HCI g-Factor Context

### H-like tin and current precision frontier

J. Morgner et al., "Stringent test of QED with hydrogen-like tin," *Nature* **622**, 53-57 (2023). DOI: https://doi.org/10.1038/s41586-023-06453-2 .

Relevance: demonstrates current high-precision H-like `g`-factor measurement and theory in a heavy-ish H-like ion, with modern correction budgets. It is context, not a prior weak-density calculation.

### Isotope shifts and new physics with bound-electron g factors

T. Sailer et al., "Measurement of the bound-electron g-factor difference in coupled ions," *Nature* **606**, 479-483 (2022). DOI: https://doi.org/10.1038/s41586-022-04807-w .

What they did: measured an isotope shift in H-like Ne with high precision and used it to constrain new scalar-boson physics. Their theory discussion includes finite nuclear size, recoil, and new-physics scalar contributions.

Relation: important context for why isotope shifts of bound-electron `g` factors are interesting. It is not a calculation of the Standard-Model weak neutral-current correction with finite p/n densities.

### 2025 isotope-shift theory limits

D. S. Akulov, R. R. Abdullin, D. V. Chubukov, D. A. Glazov, and A. V. Volotka, "g-Factor Isotopic Shifts: Theoretical Limits on New Physics Search," *Atoms* **13**, 52 (2025). DOI: https://doi.org/10.3390/atoms13060052 .

What they did: analyzed isotope-shift sensitivity of H-like HCI `g` factors to hypothetical scalar-boson interactions and discussed limiting uncertainties from finite nuclear size, recoil, and nuclear polarization.

Relation: close positioning context for `g`-factor sensitivity to short-range electron-nucleus interactions. It is BSM scalar-boson focused, not a finite-density Standard-Model weak-neutral-current benchmark.

## 6. Related Weak-Density And Neutron-Skin Work

The neutron-density part of this project is better connected to atomic parity violation than to existing bound-electron `g` work.

Key source:

A. V. Viatkina, D. Antypas, M. G. Kozlov, D. Budker, and V. V. Flambaum, "Dependence of atomic parity-violation effects on neutron skins and new physics," *Phys. Rev. C* **100**, 034318 (2019). DOI: https://doi.org/10.1103/PhysRevC.100.034318 .

What they did: estimated neutron-skin and new-physics contributions to spin-independent atomic parity nonconservation, using Fermi proton/neutron density distributions for multiple heavy atoms, including Pb.

Relation: must cite for the nuclear weak-density and neutron-skin logic. However, their observable is APV, not a bound-electron `g` correction.

Recent related preprint:

A. V. Viatkina, Ch. G. Mertens, B. Ohayon, V. A. Yerokhin, and A. Surzhykov, "Atomic parity violation in highly charged 40,48Ca and 208Pb ions," arXiv:2602.19962 (2026). Search result/public abstract: https://papers.cool/arxiv/2602.19962 .

Relation: very relevant to H-like highly charged ions and neutron-skin effects, but for parity-violating E1 amplitudes, not the parity-conserving weak `g` correction.

## 7. Novelty Assessment

Category selected: **C. Formula exists but Pb/U finite-density numerical benchmarks appear absent.**

Reasoning:

- The mass-derivative/central-potential `g` machinery exists in the literature.
- Standard bound-electron `g` budgets are highly developed and include QED, recoil, finite nuclear size, nuclear polarization, and sometimes broad mention of weak/hadronic contributions.
- New-physics scalar-boson contributions to H-like `g`-factor isotope shifts have been studied.
- Finite nuclear charge-density effects on bound-electron `g` factors have been studied in detail.
- Neutron-skin and finite weak-density effects have been studied in atomic parity violation and HCI APV amplitudes.
- This search did **not** find a peer-reviewed numerical table for the Standard-Model spin-independent weak contribution to H-like Pb/U `1s` `g` factors with separate finite proton/neutron Fermi weak densities and neutron-skin sensitivity.

Red-flag caveat: modern high-precision H-like `g` theory papers may include a small "weak" correction in total budgets. If such a term is traced to a detailed finite-density calculation, the novelty category would move from C to B. At present, the accessible search results did not show such a detailed Pb/U finite p/n density benchmark.

## 8. Recommended Paper Positioning

Safest title-level framing:

```text
Standard-Model weak contribution to the bound-electron g factor of heavy H-like ions:
finite proton-neutron density and neutron-skin sensitivity
```

Safer abstract-level claim:

```text
We provide an explicit finite-density benchmark for the spin-independent
Standard-Model weak contribution to the H-like 1s bound-electron g factor
in Pb and U, and quantify its sensitivity to neutron-skin variations.
```

Avoid:

- "first calculation" unless a deeper database/author-citation search confirms no prior finite-density weak `g` table;
- presenting the mass-derivative relation as new;
- presenting the computation as experimentally observable;
- implying that neutron-skin sensitivity in `g` factors is comparable to APV sensitivity without a scale comparison.

Best positioning:

```text
Improved finite-density benchmark / independent weak-source estimate with neutron-density sensitivity.
```

If a prior simplified weak `g` estimate is found later, the paper can still be positioned as:

```text
independent reproduction plus finite proton/neutron density and neutron-skin sensitivity.
```

## 9. Must-Cite Papers

- S. G. Karshenboim, R. N. Lee, and A. I. Milstein, "g factor of an electron or muon bound by an arbitrary central potential," *Phys. Rev. A* **72**, 042101 (2005). DOI: https://doi.org/10.1103/PhysRevA.72.042101 .
- V. M. Shabaev, D. A. Glazov, G. Plunien, and A. V. Volotka, "Theory of Bound-Electron g Factor in Highly Charged Ions," *J. Phys. Chem. Ref. Data* **44**, 031205 (2015). DOI: https://doi.org/10.1063/1.4921299 .
- D. A. Glazov and V. M. Shabaev, "Finite nuclear size correction to the bound-electron g factor in a hydrogenlike atom," *Phys. Lett. A* **297**, 408-411 (2002). DOI: https://doi.org/10.1016/S0375-9601(02)00021-X .
- I. A. Valuev, Z. Harman, C. H. Keitel, and N. S. Oreshkina, "Skyrme-type nuclear interaction as a tool for calculating the finite-nuclear-size correction to atomic energy levels and the bound-electron g factor," *Phys. Rev. A* **101**, 062502 (2020). DOI: https://doi.org/10.1103/PhysRevA.101.062502 .
- J. Zatorski, N. S. Oreshkina, C. H. Keitel, and Z. Harman, "Nuclear shape effect on the g factor of hydrogenlike ions," *Phys. Rev. Lett.* **108**, 063005 (2012). DOI: https://doi.org/10.1103/PhysRevLett.108.063005 .
- A. V. Nefiodov, G. Plunien, and G. Soff, "Nuclear-polarization correction to the bound-electron g factor in heavy hydrogenlike ions," *Phys. Rev. Lett.* **89**, 081802 (2002). DOI: https://doi.org/10.1103/PhysRevLett.89.081802 .
- J. Morgner et al., "Stringent test of QED with hydrogen-like tin," *Nature* **622**, 53-57 (2023). DOI: https://doi.org/10.1038/s41586-023-06453-2 .
- T. Sailer et al., "Measurement of the bound-electron g-factor difference in coupled ions," *Nature* **606**, 479-483 (2022). DOI: https://doi.org/10.1038/s41586-022-04807-w .
- D. S. Akulov et al., "g-Factor Isotopic Shifts: Theoretical Limits on New Physics Search," *Atoms* **13**, 52 (2025). DOI: https://doi.org/10.3390/atoms13060052 .
- A. V. Viatkina et al., "Dependence of atomic parity-violation effects on neutron skins and new physics," *Phys. Rev. C* **100**, 034318 (2019). DOI: https://doi.org/10.1103/PhysRevC.100.034318 .

## 10. Red Flags

- Search results and recent correction-budget discussions mention "weak" contributions in broad bound-electron `g` contexts. These may refer to the free-electron electroweak anomaly contribution or to a small leading bound-state correction. Before claiming "first," trace every "weak" item in recent H-like `g` tables, especially for Sn and projected heavy ions.
- The 2025 H-like tin precision theory context is moving fast. Any paper draft should update the bound-electron `g` correction hierarchy using the latest theory values.
- New-physics scalar-boson `g`-factor isotope-shift papers are close in observable class. The manuscript must sharply distinguish Standard-Model weak neutral current from BSM scalar exchange.
- Atomic parity violation literature already treats finite neutron densities and neutron skins carefully. The manuscript should borrow that language and avoid implying neutron-density sensitivity is new in atomic physics generally.
- The project's numerical values are far below major nuclear-size and QED correction scales. The paper should be framed as a Standard-Model hierarchy/benchmark calculation, not as an experimental proposal.

## 11. Bottom Line

No direct peer-reviewed prior calculation was found in this search for H-like Pb/U finite proton-neutron-density Standard-Model weak `g` corrections with neutron-skin sensitivity.

The most defensible novelty claim is:

```text
an explicit finite-density benchmark and neutron-skin sensitivity study
for the spin-independent Standard-Model weak contribution to the H-like
1s bound-electron g factor in Pb and U.
```

The safest paper type is a compact benchmark/positioning paper, not a broad discovery or observability paper.
