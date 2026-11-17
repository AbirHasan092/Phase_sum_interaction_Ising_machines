# Four-stage data-driven phase-sum verification

Choose gain modulation when running MATLAB: `R=main([],1);` uses the first harmonic, and `R=main([],10);` uses all ten image-derived harmonics. `R=main;` defaults to 1. A custom data path can replace `[]`. `R.Kmax` records the chosen count. Each run replaces the standard result files and figures with that configuration. The independent Python evaluator currently remains configured for the one-harmonic case.

CURRENT CONFIGURATION: gain modulation uses only k=1, |G1|=0.85, so g(t)=1.7*cos(2*omega*t), with zero phase and DC. The CSV and NPZ results have been regenerated and all 28 independent temporal averaging checks passed. At M=8: RMS curl 15.0430384, maximum curl 30.1073335, normalized curl 0.08167526, cross-power fraction 0.0006084876. Run `R=main;` to regenerate MATLAB MAT/FIG/PNG outputs; previously generated MATLAB outputs remain from their last run. The original ten-harmonic assumptions and numerical table below are retained as historical context and do not describe this current configuration.

Run `R = main;` in MATLAB from this folder. Input defaults to `../Diffring_WV_PPV.matlab`. Original data and existing analysis are untouched. Results are written into `results/`. The supplied CSV results were independently computed with NumPy; MATLAB itself could not be run because MathWorks service startup failed in this session. Running main creates the MATLAB result structure and comparison figure.

## Explicit model assumptions

The export contains one complete cycle, endpoints included, with PPV time/value and waveform time/value columns. PPV is used as iPRC without rescaling. A four-stage differential ring is modeled as eight equally spaced taps, each an exact shifted copy of this single trace. Actual tap-specific behavior cannot be measured from this file. Frequencies are in inverse exported X units; absolute curl is in the resulting exported coupling units per radian, not calibrated physical units.

The image supplies amplitudes only. Approximate manually digitized positive complex coefficients |G_k|, k=1..10, are [0.85, 0.072, 0.039, 0.020, 0.010, 0.0055, 0.0030, 0.0016, 0.00090, 0.00048]. They are approximate, editable in main.m, and are not exact underlying source data. All phases are assumed zero; G0=0 because DC was not shown. The convention follows PhaseSumVerification: g(tau)=2 sum_k |G_k| cos(2 k tau). Thus plot harmonic k is interpreted as pump harmonic k at 2*k*omega, and its cosine amplitude is 2|G_k|. Different phases, DC, or fundamental-frequency convention change the result. Signed gain is allowed.

## Cutoff and method

All waveform, PPV, and gain coefficients receive a -100 dBc cutoff, referenced independently to the strongest non-DC coefficient in each signal. This is equivalent to power ratio >=1e-10 or amplitude ratio >=1e-5. DC is preserved only if it passes the same coefficient comparison; gain DC is zero. Interaction coefficients are not thresholded again before curl evaluation, avoiding hiding small derivative contributions.

Fourier coefficients are integrated over actual nonuniform times with trapezoidal quadrature. Both interval endpoints are retained. The candidate harmonic limit is floor((601-1)/4)=150, with 26 waveform and 87 PPV orders retained at the specified cutoff. The weak high-order estimates are limited by input sampling/quadrature accuracy; -100 dBc retention is a numerical selection rule, not a validated measurement noise floor.

For each M=2..8, nearest distinct taps among the eight available offsets are paired complementarily. An ideal comparator uses M perfectly uniformly spaced taps, which need not be realizable by the eight-tap ring. Exact subsets occur at M=2,4,8. Nearest selection is greedy, with deterministic column-major tie ordering matching MATLAB.

Temporal resonance requires n+m+2k=0. The forward coefficient is G_k times mean_r[Z_n Q_m exp(i*n*delta_p+i*m*delta_s)]. The reverse coefficient uses Q_n Z_m with the same tap factor. For ideal uniform paths, the tap factor is 1 when n-m is a multiple of M, otherwise 0. Two-dimensional inverse Fourier transforms give both phase-velocity surfaces; spectral differentiation gives curl = dF1/dtheta2 - dF2/dtheta1. Conductance self terms depend only on the receiving phase and contribute zero to this cross-derivative curl, so are omitted.

R_cross is off-diagonal coefficient mean-square divided by all non-DC interaction mean-square. Normalized curl is RMS(curl)/sqrt(mean((dF1/dtheta2)^2+(dF2/dtheta1)^2)). Nonzero curl rules out gradient dynamics under this model. Zero curl alone does not establish a globally periodic potential.

## Measured results for realizable eight-tap model

| M paths | Cross-power fraction | Curl RMS | Normalized curl |
|---|---:|---:|---:|
| 2 | 0.28910 | 124.309 | 0.57129 |
| 3 | 0.04373 | 43.802 | 0.23275 |
| 4 | 0.24073 | 123.593 | 0.56917 |
| 5 | 0.01658 | 28.911 | 0.15574 |
| 6 | 0.04373 | 43.802 | 0.23275 |
| 7 | 0.00882 | 23.185 | 0.12537 |
| 8 | 0.00061344 | 15.061 | 0.08177 |

At M=8 maximum absolute curl is approximately 31.470. Curl remains nonzero because finite tap sampling retains aliases with n-m a nonzero multiple of 8. Neither increasing M nor reducing cross power guarantees a monotonic normalized curl.

## Verification

`verify_and_measure.py` independently implements the Fourier calculation and checks direct periodic time averaging at two phase points for every M and both tap arrangements (28 comparisons). All checks passed. It creates CSV tables and `independent_results.npz`. It requires NumPy. This validates the algebra and implementation under the stated assumptions, not the unprovided modulation phases or PPV calibration. This is an averaged weak-coupling analysis, not a full coupled-circuit simulation.
