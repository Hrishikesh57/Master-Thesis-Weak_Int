# Milestone Roadmap

## M1 Reproduce One Existing H-like Weak g-Factor Result

- identify one documented local thesis or paper weak target for an H-like ion
- implement the calculation through reused DKB orbitals and weak-potential machinery
- compare the numerical result against the documented target and record the source

## M2 Compute Relativistic Weak Contribution For H-like `1s` Systematically

- evaluate Pb and U first
- scan numerical controls and document stability
- keep units, weak normalization, and sign conventions explicit

## M3 Compute Li-like `2s` Baseline In A Consistent Framework

- reuse the copied Li-like `f2py` suite
- collect baseline ingredients in the same project-facing data structure
- document any remaining ambiguity in mapping the copied lithium workflows into the final paper observable

## M4 Extract Finite Nuclear Size Response Numerically

- vary the RMS radius in a controlled finite-difference protocol
- compute `\partial g / \partial R` or the equivalent discrete response
- quantify step-size dependence

## M5 Determine `\xi` To Cancel Leading FNS Contribution

- use the same nucleus, same charge-model convention, and same differentiation protocol in both charge states
- construct the first `\xi` from Dirac+FNS only
- record cancellation checks explicitly

## M6 Test Whether The Weak Residual Survives In `\delta_\xi g`

- first study the case where only the H-like weak contribution is available
- document that as a partial residual estimate
- extend later once a Li-like weak contribution exists in a consistent framework

## M7 Generate Publication-Quality Tables And Plots

- absolute baseline ingredients
- finite nuclear size response
- `\xi`
- `\delta_\xi g`
- weak residual estimates and caveats
