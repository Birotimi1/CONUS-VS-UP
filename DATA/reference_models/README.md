# ECM1 crust + AK135 reference mantle

The selected ECM1 cell is **19178**, **36.5 N, 82.5 W**, a continental
Appalachian column with **10 m of sediment** and a **43.7731 km Moho**.
It is a reproducible thin-sediment test column, not a fitted TA station.
The TSV files preserve both published source rows. The JSON records official
download URLs, archive/table SHA-256 hashes, units and rounding corrections.
Re-extract with `scripts/extract_ecm1_column.py` and the two official archives.

`ak135_mantle.csv` transcribes the continental velocity values at 35, 77.5,
120 km and the deeper values at 165, 210, 260 km from the
[ANU AK135 release](https://rses.anu.edu.au/seismology/ak135/ak135f.html).
At 120 km we retain the continental Vp=8.0500, not the spherical-average
8.0505. The reference is Kennett, Engdahl & Buland (1995).
We interpolate Vs linearly and sample it at the midpoint of five finite
mantle layers (Moho–77.5–120–165–210–250 km). The half-space uses Vs at250km.
This is a labelled piecewise-constant approximation to AK135, not a local
mantle measurement or an ECM1 depth profile. AK135 does not introduce an LVZ.

The matched inversion test preserves ECM1 thicknesses and Vs and uses the
inversion's THBI material closure for Vp/density. Its crustal Vp/Vs is the
ratio of published harmonic-mean crustal Vp and Vs. Sediment Vp/rho and
crustal rho use the imported reference relations; mantle Vp/rho use the
reference's depth-dependent relations. **These derived properties are not
the published ECM1/AK135 Vp and density.** Raw ECM1 values and AK135 velocity
nodes remain separate for inspection and future model-mismatch tests.

The `reference100_vti` synthetic additionally assigns crustal xi=1.06,
phi=0.94 and eta=1; sediments and mantle are isotropic. These are **assigned
synthetic anisotropy values**, not properties published in ECM1 or AK135.
The retained Vs/Vp parameters use THBI's velocity-average convention; see
`functions/vti_elastic.m` for their relation to Vsv/Vsh/Vpv/Vph.
