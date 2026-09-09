# Trajectory and weak freeze-out analysis

This directory contains post-processing tools for tracer-particle trajectories,
with a particular focus on the evolution and freeze-out of the electron fraction
$Y_e$.

The main Python workflow is

```text
analyze_weak_freezeout.py
        |
        | selected tracer IDs / particle tables
        v
compute_tracer_percentiles.py
        |
        | temperature-dependent percentile tables
        v
plot_weak_evolution.py
```

The Python scripts are kept in

```text
particle-tracing-in-3d/analysis_trajectory/py/
```

They can either be run directly from the repository or copied into an individual
model-analysis directory, as described below.

## 1. Setup

### 1.1 Running the scripts directly from the repository

The recommended setup is to keep a single copy of the analysis scripts in the
Git repository. Define the root directory of this repository, for example

```bash
export PTR_ROOT=/path/to/particle-tracing-in-3d
export PTR_PY=$PTR_ROOT/analysis_trajectory/py
```

Then move to the directory containing the tracer-analysis data and run the
scripts from the repository:

```bash
cd /path/to/model/analysis
python3 $PTR_PY/analyze_weak_freezeout.py ...
```

`tracer_analysis_utils.py` is located in the same directory as the main Python
scripts. When a script is executed as above, Python automatically adds that
script directory to its module search path, so there is no need to copy
`tracer_analysis_utils.py` into each analysis directory or to set `PYTHONPATH`.

### 1.2 Alternative: copying the scripts into the analysis directory

It is also possible to use the analysis tools without defining `PTR_ROOT` or
`PTR_PY`. In this case, copy the required Python files from

```text
particle-tracing-in-3d/analysis_trajectory/py/
```

to the model-analysis directory. For the workflow described in this README, the
required files are

```text
analyze_weak_freezeout.py
compute_tracer_percentiles.py
plot_weak_evolution.py
tracer_analysis_utils.py
```

All four files should be kept in the same directory because
`analyze_weak_freezeout.py` and `plot_weak_evolution.py` import
`tracer_analysis_utils.py`.

After copying them, the scripts can be executed directly, for example

```bash
cd /path/to/model/analysis
python3 analyze_weak_freezeout.py ...
python3 compute_tracer_percentiles.py ...
python3 plot_weak_evolution.py ...
```

This local-copy setup is convenient for a self-contained analysis directory,
whereas the `PTR_ROOT` / `PTR_PY` setup above is preferable when the same scripts
are used for many models because only one repository copy needs to be updated.

A typical analysis directory contains

```text
stat_traj.dat

weak_freezeout_dye.dat
weak_freezeout_timescale.dat
weak_freezeout_dng.dat
weak_freezeout_dnv.dat

traj_time
traj_Ye
traj_yecap
traj_texp
traj_caprate
traj_eta
traj_r
traj_vel
traj_avtexp
```

The `traj_*` temperature tables are generated from the tracer trajectories by
the trajectory-analysis code (in particular `src/make_table.F90`).

## 2. Weak freeze-out analysis

The main freeze-out analysis is performed with

```text
analyze_weak_freezeout.py
```

This script

- selects the ejecta population;
- compares several definitions of weak freeze-out;
- constructs freeze-out-temperature distributions;
- writes selected tracer IDs and particle tables; and
- produces diagnostic histograms and 5-GK scatter plots.

It reads

```text
stat_traj.dat
weak_freezeout_dye.dat
weak_freezeout_timescale.dat
weak_freezeout_dng.dat
weak_freezeout_dnv.dat
```

The four freeze-out files correspond to different definitions:

```text
weak_freezeout_dye.dat
    remaining net change in Ye

weak_freezeout_timescale.dat
    weak-interaction timescale versus expansion timescale

weak_freezeout_dng.dat
    remaining gross weak processing

weak_freezeout_dnv.dat
    remaining absolute variation in Ye
```

The current analysis starts from positive-mass tracers satisfying the Bernoulli
criterion used in the code (`bernoulli < 0`) and then applies the requested
additional cuts.

### 2.1 time cut

```bash
--time-rfl-min TMIN
--time-rfl-max TMAX
```

These options select tracers according to `time_rfl`, the time at which the
tracer crosses the extraction radius rfl, read from the
`t_last_rfl` column of `stat_traj.dat`.

For example,

```bash
--time-rfl-min 0.05 --time-rfl-max 0.10
```

selects

```text
0.05 s <= time_rfl <= 0.10 s
```

by default. The default lower bound is `0.0 s`; there is no default upper
bound.

If the `t_last_rfl` column cannot be identified from the header of
`stat_traj.dat`, its 1-based column number can be supplied explicitly:

```bash
--time-rfl-column N
```

### 2.2 Angular cut

```bash
--theta-min THETA_MIN
--theta-max THETA_MAX
```

The default definition is

```bash
--theta-mode axis
```

with

```text
theta = atan2(sqrt(x_fin^2 + y_fin^2), abs(z_fin))
```

so that `theta` ranges from 0 to 90 degrees. For example,

```bash
--theta-max 45
```

selects tracers within 45 degrees of either the positive or negative z axis.

Other available definitions are

```text
--theta-mode polar
--theta-mode latitude
```

but `axis` is normally used for a polar-ejecta selection.

### 2.3 Electron-fraction cut at 5 GK

```bash
--ye-5gk-min YE_MIN
--ye-5gk-max YE_MAX
```

These options are applied to

```text
Ye(5 GK)
```

from `stat_traj.dat`. For example,

```bash
--ye-5gk-min 0.20 --ye-5gk-max 0.35
```

selects

```text
0.20 <= Ye(5 GK) <= 0.35
```

### 2.4 Entropy cut at 5 GK

```bash
--s-5gk-min S_MIN
--s-5gk-max S_MAX
```

These options are applied to

```text
s(5 GK) / k_B
```

in units of `k_B` per nucleon. For example,

```bash
--s-5gk-min 10 --s-5gk-max 100
```

selects

```text
10 <= s(5 GK)/k_B <= 100
```

### 2.5 Final radial-velocity cut

```bash
--vr-min VR_MIN
--vr-max VR_MAX
```

These options are applied to the final radial velocity `v^r_fin` in cm/s.
For example,

```bash
--vr-min 1.0e9
```

selects tracers with

```text
v^r_fin >= 1.0e9 cm/s
```

### 2.6 Other useful options

Additional cuts include

```text
--time-fin-min / --time-fin-max
--r-fin-min    / --r-fin-max
```

By default, upper bounds are inclusive (`value <= maximum`). Use

```bash
--open-upper
```

to make all upper bounds strict (`value < maximum`).

A particular freeze-out-temperature range can be inspected with

```bash
--inspect-T-range TMIN TMAX
```

or the histogram bin containing a specified temperature can be selected with

```bash
--inspect-T T
```

The freeze-out definition used for the corresponding selected-particle table is
controlled by

```text
--inspect-source dye
--inspect-source time
--inspect-source both
```

with `dye` as the default.

### 2.7 Typical freeze-out selection

For example, to select relatively late reflected polar ejecta,

```bash
python3 $PTR_PY/analyze_weak_freezeout.py \
    --time-rfl-min 0.05 \
    --time-rfl-max 0.10 \
    --theta-max 45 \
    -o weakfo
```

A selection can additionally be restricted in 5-GK properties and final radial
velocity, for example

```bash
python3 $PTR_PY/analyze_weak_freezeout.py \
    --time-rfl-min 0.05 \
    --time-rfl-max 0.10 \
    --theta-max 45 \
    --ye-5gk-min 0.20 \
    --ye-5gk-max 0.35 \
    --s-5gk-min 10 \
    --s-5gk-max 100 \
    --vr-min 1.0e9 \
    -o weakfo
```

Important outputs include

```text
weakfo_common_selected_ids.txt
weakfo_selected_dye_particles.txt
weakfo_selected_dye_trajectories.txt
weakfo_TFO_histograms.txt
```

as well as freeze-out-temperature histograms and 5-GK diagnostic plots.

When no `--inspect-T` or `--inspect-T-range` is specified,
`weakfo_selected_dye_particles.txt` contains the tracers passing all active cuts
for which the `dye` freeze-out quantities are finite. This is normally used as
the PID source for the next step.

## 3. Computing temperature-dependent tracer percentiles

Use

```text
compute_tracer_percentiles.py
```

to compute mass-weighted statistics for a selected tracer population.

The input PID file can be either a simple one-column PID list or a wider
`*_particles.txt` table. Only the first non-comment column is interpreted as the
tracer ID.

For the population selected above, the standard use is

```bash
python3 $PTR_PY/compute_tracer_percentiles.py \
    --id-file weakfo_selected_dye_particles.txt \
    -o weakfo
```

The default mass-weighted percentiles are

```text
15.865 %
50.000 %
84.135 %
```

corresponding approximately to the median and the +/-1-sigma interval.

The script reads

```text
stat_traj.dat

traj_time
traj_texp
traj_caprate
traj_Ye
traj_yecap
traj_eta
traj_r
traj_vel
traj_avtexp
```

and produces

```text
weakfo_ye.txt
weakfo_yecap.txt
weakfo_texp.txt
weakfo_caprate.txt
weakfo_eta.txt
weakfo_r.txt
weakfo_vr.txt
weakfo_time.txt
weakfo_avtexp.txt
weakfo_ids.txt
```

These files are read directly by `plot_weak_evolution.py`.

### 3.1 Division of responsibility between the two selection stages

Both `analyze_weak_freezeout.py` and `compute_tracer_percentiles.py` can apply
cuts, but their intended roles are different.

The recommended convention is:

```text
Cuts defining the population used for the freeze-out analysis
    -> analyze_weak_freezeout.py

Optional further subdivision of an already selected PID population
    -> compute_tracer_percentiles.py
```

In particular, if the `T_FO` histograms and the temperature-dependent
percentile curves are intended to describe exactly the same tracer population,
apply all relevant cuts in `analyze_weak_freezeout.py` and do not apply an
additional cut in `compute_tracer_percentiles.py`.

Thus the normal second step is simply

```bash
python3 $PTR_PY/compute_tracer_percentiles.py \
    --id-file weakfo_selected_dye_particles.txt \
    -o weakfo
```

#### Optional final-state cuts in `compute_tracer_percentiles.py`

For exploratory subdivision, `compute_tracer_percentiles.py` supports

```text
--ye-min / --ye-max    -> Ye_fin
--s-min  / --s-max     -> s_fin
--vr-min / --vr-max    -> v^r_fin
```

#### Optional 5-GK cuts in `compute_tracer_percentiles.py`

Alternatively it supports

```text
--ye-5gk-min   / --ye-5gk-max   -> Ye(5 GK)
--s-5gk-min    / --s-5gk-max    -> s(5 GK)
--texp-5gk-min / --texp-5gk-max -> t_exp(5 GK)
```

The final-state and 5-GK cut families in `compute_tracer_percentiles.py` are
mutually exclusive.

For example, to examine only the `0.20 <= Ye(5 GK) <= 0.30` subset of an
already selected population,

```bash
python3 $PTR_PY/compute_tracer_percentiles.py \
    --id-file weakfo_selected_dye_particles.txt \
    --ye-5gk-min 0.20 \
    --ye-5gk-max 0.30 \
    -o weakfo_Ye020_030
```

In this case the resulting percentile curves describe a narrower population
than the original `weakfo_TFO_histograms.txt`. Therefore the original TFO
histogram should not be overplotted as if it represented exactly the same
population unless the same cut was already applied in
`analyze_weak_freezeout.py`.

## 4. Plotting the weak evolution

The final figures are produced by

```text
plot_weak_evolution.py
```

using the percentile tables generated above.

For example,

```bash
python3 $PTR_PY/plot_weak_evolution.py \
    --prefix weakfo \
    --tfo-hist-file weakfo_TFO_histograms.txt
```

The output includes figures showing

- `Ye(T)` and the capture-equilibrium electron fraction;
- `Ye(T) - Ye_eq(T)`;
- weak-interaction and expansion timescales;
- electron degeneracy;
- radius;
- radial velocity; and
- elapsed time.

The freeze-out-temperature distributions from `analyze_weak_freezeout.py` are
overlaid when `--tfo-hist-file` is supplied.

Important outputs include

```text
Ye-weakfo.pdf
Timescale-weakfo.pdf
Yediff-weakfo.pdf
eta-weakfo.pdf
r-weakfo.pdf
vr-weakfo.pdf
```

## 5. Standard workflow

A typical complete analysis using the scripts directly from the repository is

```bash
export PTR_ROOT=/path/to/particle-tracing-in-3d
export PTR_PY=$PTR_ROOT/analysis_trajectory/py

cd /path/to/model/analysis

python3 $PTR_PY/analyze_weak_freezeout.py \
    --time-rfl-min 0.05 \
    --time-rfl-max 0.10 \
    --theta-max 45 \
    -o weakfo

python3 $PTR_PY/compute_tracer_percentiles.py \
    --id-file weakfo_selected_dye_particles.txt \
    -o weakfo

python3 $PTR_PY/plot_weak_evolution.py \
    --prefix weakfo \
    --tfo-hist-file weakfo_TFO_histograms.txt
```

If the four Python files listed in Sec. 1.2 have instead been copied into the
analysis directory, the corresponding workflow is simply

```bash
python3 analyze_weak_freezeout.py \
    --time-rfl-min 0.05 \
    --time-rfl-max 0.10 \
    --theta-max 45 \
    -o weakfo

python3 compute_tracer_percentiles.py \
    --id-file weakfo_selected_dye_particles.txt \
    -o weakfo

python3 plot_weak_evolution.py \
    --prefix weakfo \
    --tfo-hist-file weakfo_TFO_histograms.txt
```

The repository-based setup keeps the analysis scripts under Git control in a
single location, while the local-copy setup provides a fully self-contained
analysis directory.
