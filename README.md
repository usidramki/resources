# From Thebes to Babylon: Resource distribution from noncombatants in complex conflicts

This repository accompanies the manuscript by S. Shelyag, J. C. McKerral, Y. Shi, M. Francis, T. McLennan-Smith and M. Zuparic. It provides the computational code used to produce the numerical results reported in the study.

## Model overview

The model couples:
- Attrition between two opposing forces (Blue and Red).
- Production, transfer and consumption of resources.
- Evolving support from a heterogeneous noncombatant population (Green).

The study examines victory, defeat and resource-driven conflict stall using the full model, a large-affinity asymptotic approximation and a further two-dimensional reduction. Historical applications are confined to qualitative illustrations based on selected pre-modern sieges from Alexander the Great's campaign against the Achaemenid Empire.

## Software requirements and dependencies

- **Reference runtime:** IDL 8.8.3.
- **External libraries or packages:** None.

The code is self-contained and requires no third-party IDL libraries, add-on packages or separately installed numerical solvers. The source supplied here is the implementation used for the reported computations.

## Numerical methods and reproduction

The differential equations are integrated using an in-house implementation of the embedded Runge–Kutta–Fehlberg 5(6) method.

Model equations, parameter choices and initial conditions are described in the manuscript and the relevant figure captions. Numerical implementation details are provided in the paper. Reproduction of a particular result requires the corresponding model variant, parameter values and initial conditions.

## Alternative runtime

[GNU Data Language (GDL)](https://github.com/gnudatalanguage/gdl) is a free, open-source IDL-compatible implementation. The reference environment used to generate the published results is IDL 8.8.3.

## Citation

If you use this code, please cite the accompanying manuscript:

S. Shelyag, J. C. McKerral, Y. Shi, M. Francis, T. McLennan-Smith and M. Zuparic, *From Thebes to Babylon: Resource distribution from noncombatants in complex conflicts*.

