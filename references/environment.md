# The environment: what is installed, what to do without it, and the output language

## Three roles a language can play

Sort every language question into one of three roles before answering it, because only one of them is fixed.

- **Host.** The language the analysis runs in. This is always Python: `bimgrl` is a Python package and so is gotranx, which supplies the symbolic Jacobian and generates every scheme. Neither can be ported to another language in a session, and neither is a fallback for the other.
- **Source.** The language the user's model is written in. It is needed exactly once, to produce the files `fidelity` checks the translation against (`translation.md`), or to serve the RHS adapter on the native path (`contracts.md`).
- **Target.** The language the user wants the stepper in. It is needed exactly once, at the end, to run the frozen stepper and write a trace for `bimgrl validate`.

Everything between those two ends (structure, invariants, reference, stiffness, blocks, sweeps, tuning) happens on the canonical `.ode` in the host. So a source or target language the agent cannot execute does not block the pipeline: it costs one file handoff per end. What does block it is a missing host, which is why the host is provisioned first.

## Preflight

```
sh <skill-dir>/scripts/env-check.sh <project-dir>
```

A POSIX shell script, so it runs before any Python exists. It prints which of uv, Python, a C compiler, Julia, MATLAB and Octave are present, whether the skill's environment is built, and which language markers the project directory carries. It only reports; the rules below decide. Run it once at the start and keep the two tables in view for the gate.

## The host: no uv, no Python

1. **uv is present.** Go to Setup in `SKILL.md`. uv downloads a Python of its own when the system has none or has one older than 3.11, so a missing `python` row is not a problem by itself.
2. **uv is missing.** Ask the user for permission to install it, saying that it goes into their home directory and needs no administrator rights: `curl -LsSf https://astral.sh/uv/install.sh | sh` (Windows: `powershell -ExecutionPolicy ByPass -c "irm https://astral.sh/uv/install.ps1 | iex"`). Then Setup.
3. **The user declines, the machine is offline, or the install fails.** Stop. Show the Tools table and what is missing. Do not translate the scripts into another language and do not estimate by reasoning what they would have printed: the skill's rule is that every number comes from a script, and without a host there are no scripts. The way forward is another machine for the analysis; the source and target ends still work by file handoff from this one.

The same stop applies when `bimgrl setup` cannot reach the gotranx fork (`config.toml`): report the message and the `GOTRANX_SRC` alternative.

## No C compiler: the generated-stepper path

Only `sweep` and `time` need a compiler, because they compile the generated C scheme and call it through ctypes. Every phase before the baseline (intake, structure, invariants, reference, stiffness, propose) is unaffected. Without a compiler, the same scheme is generated in Python instead, stepped by a short driver that writes the trace contract, and judged by `bimgrl validate`, which is the same metric code the sweep uses. On TT06 the Python and the compiled C stepper give the same biomarkers to 1e-11 relative, so nothing about the conclusions changes; what changes is the speed and the bookkeeping.

**Generate one stepper per scheme.** `bimgrl codegen` maps a label to a gotranx scheme and takes the base and the frozen boundary from the blocks file, so the labels of the sweep are reproduced like this (use absolute paths: `uv run --directory` changes the working directory to the skill):

| label | command |
|---|---|
| `RL` | `bimgrl codegen <model_dir> --backend python --scheme RL -o <out> -- --stiff-states g1 --stiff-states g2 ...` with the gates of `hh_gates` in `bimgrl/structure.json` (or `rl_gates` of `meta.toml` if set) |
| `GRL` | `bimgrl codegen <model_dir> --backend python --scheme GRL -o <out>` |
| `BIRL-M`, `BIMGRL` | `bimgrl codegen <model_dir> --backend python --scheme BIRL-M -o <out>` with the `blocks.toml` that `propose` wrote for that base |
| `BIRL`, `BIGRL` | the same with `--blocks` pointing at a copy of `blocks.toml` whose `[[substep_blocks]]` are removed |
| `BIRL-M-unfrozen` | the same with a copy that sets `freeze_boundary_fluxes = false` |

**Drive it.** Save this once as `<model_dir>/bimgrl/drive.py` and run it with the skill's Python (`<skill-dir>/.venv/bin/python`). It reads the stimulus from `meta.toml`, runs protocol B, and writes the last beat on the 0.1 ms grid, the state at every beat start, and the run meta:

```python
#!/usr/bin/env python3
"""usage: drive.py <generated.py> <meta.toml> <label> <dt> <beats> <out_stem>"""
import importlib.util, json, sys, time, tomllib
import numpy as np

gen_path, meta_path, label, dt, beats, stem = sys.argv[1], sys.argv[2], sys.argv[3], float(sys.argv[4]), int(sys.argv[5]), sys.argv[6]
spec = importlib.util.spec_from_file_location("gen", gen_path)
gen = importlib.util.module_from_spec(spec)
spec.loader.exec_module(gen)
step = next(getattr(gen, n) for n in ("block_rush_larsen", "hybrid_rush_larsen", "generalized_rush_larsen", "explicit_euler") if hasattr(gen, n))
meta = tomllib.load(open(meta_path, "rb"))
stim = dict(meta["stim"], period=meta["bcl"])
p = gen.init_parameter_values(**{f"stim_{k}": v for k, v in stim.items()})
y = gen.init_state_values()
names = sorted(gen.state, key=gen.state.get)        # the generated order, which the reference uses too

per_beat = round(meta["bcl"] / dt)
every = max(1, round(0.1 / dt)) if dt <= 0.1 else 1  # the last beat on the 0.1 ms grid
t_beats, y_beats, rows, status = [], [], [], "ok"
start = time.perf_counter()
for beat in range(beats):
    t_beats.append(beat * meta["bcl"]); y_beats.append(y.copy())
    last = beat == beats - 1
    if last:
        rows.append([beat * meta["bcl"], *y])
    for it in range(per_beat):
        y = np.asarray(step(y, it * dt, dt, p), dtype=float)   # t restarts every beat: the stimulus is periodic
        if last and (it + 1) % every == 0:
            rows.append([beat * meta["bcl"] + (it + 1) * dt, *y])
    if not np.all(np.isfinite(y)):
        status = "nan"
        break
wall = time.perf_counter() - start
t_beats.append(len(t_beats) * meta["bcl"]); y_beats.append(y.copy())

header = ",".join(["t", *names])
if rows:
    np.savetxt(stem + ".csv", np.array(rows), delimiter=",", header=header, comments="", fmt="%.17g")
np.savetxt(stem + "_beats.csv", np.column_stack([t_beats, y_beats]), delimiter=",", header=header, comments="", fmt="%.17g")
json.dump(dict(model=meta["id"], scheme=label, dt=dt, reference=False, solver={}, language="python", t_unit="ms", protocol="B",
               beats=beats, bcl=meta["bcl"], stimulus=stim, status=status, message="", n_steps=per_beat * beats,
               wall_s_total=wall, wall_s_per_beat=wall / beats, parameters={}, float_type="float64", blocks="blocks.toml",
               created=time.strftime("%Y-%m-%dT%H:%M:%S")), open(stem + ".json", "w"), indent=2)
print(f"{label} dt {dt:g}: {status}, {wall / beats:.3g} s per beat (python)")
```

```
<skill-dir>/.venv/bin/python <model_dir>/bimgrl/drive.py <generated.py> <model_dir>/meta.toml BIRL-M 0.1 10 <model_dir>/bimgrl/traces/BIRL-M_py_dt0.1
bimgrl validate <model_dir> --test <model_dir>/bimgrl/traces/BIRL-M_py_dt0.1.csv
```

**Walk the grid by hand, in the sweep's order.** Coarsest first over 0.5, 0.2, 0.1, 0.05, 0.02, 0.01, 0.005 (only steps that divide the stimulus duration), one `validate` per step, and stop after two consecutive tight passes. A run whose status is `nan`, or that `validate` classes `unphysical`, is a stability failure; a run that is `ok` and fails a threshold is an accuracy failure. Quote the tables `validate` prints; the judgment at the gates is unchanged. Nothing sizes these runs to the session clock, because the driver is not a `bimgrl` command: read the `clock:` line that ends each `validate` before starting the next step, estimate the step from the `wall_s_per_beat` of the previous one (halving the step doubles it), do not start a run that would cross the 15-minute mark, and apply the ladder of SKILL.md by hand (under 30 minutes, BIRL-M only).

**What this path does not give you, and what to say instead.**

- No stability bisection. Report the limit as a bracket: the coarsest grid step that ran and the next coarser one that did not. Compare the bracket with the phase-4 prediction.
- No `runs.csv`, so no `bimgrl report` and no `bimgrl time`. The hand-back is the `validation_*.md` files, one per scheme and step, and the speed-up is stated from them.
- Cost is `wall_s_per_beat` of the run meta, in Python. Ratios between schemes at their passing steps are meaningful (same interpreter, same machine); the absolute numbers are not comparable with the survey's compiled-C timings and must be labelled as Python. A step in Python is about a hundred times slower than in C: TT06 takes 0.04 to 0.15 ms per step depending on the scheme, so 10 beats at 0.005 ms is several minutes, and a 90-state model proportionally more. Start coarse.
- If Julia is installed and the target is Julia, the same driver can be written in Julia against the `--backend julia` output (`init_state_values!`, `init_parameter_values!`, `block_rush_larsen(states, t, dt, parameters, values)`, `parameter_index` for the `stim_*` entries) for a much faster run. The contract and the validation are the same.

## Gate P: the output language

Settle the target before phase 0, because it decides what `codegen` emits at the end and whether a hand translation is needed.

The target is **obvious**, and is used without asking, only when one of these holds:

- the user named it ("give me a Julia stepper", "I need this in CUDA");
- the Project markers table shows the markers and the source files of exactly one language, and the user's model file is in that same language (a Julia package with the model in `src/`, a Python project with the model as a module).

Anything else is a question to the user: markers of two languages, a bare model file in a scratch directory, a CellML or `.ode` input (which says nothing about where the stepper will run), a MATLAB model inside a Python repository. Ask once, with the languages gotranx generates (C, C++, CUDA, Python, Julia, CytoZoo) and MATLAB as options, the ones installed on this machine first. Then record it, once intake has written `meta.toml` (the log needs it):

```
bimgrl log <model_dir> --section Source --line "target language: <language> (named by the user | inferred from <markers> | chosen at Gate P)"
```

The target does not change the analysis. Whatever it is, the blocks are found, swept and frozen on the canonical `.ode`; the target only decides the last `codegen` and its check.

## MATLAB

MATLAB is common among the scientists who own these models and it is the one target gotranx does not generate, so it gets a path of its own. Two facts shape it. The agent can usually run MATLAB after all, and when it cannot, the user runs one line twice.

### The run ladder

Whenever a `.m` file has to be executed, take the first rung that exists (the Tools table says which):

1. **MATLAB itself**, headless: `matlab -batch "cd('<dir>'); <call>"` (R2019a and later). A licence error here means this rung is absent, not that the user must fix it now.
2. **Octave**: `octave --no-gui --quiet --eval "cd('<dir>'); <call>"`. Plain arithmetic `.m` code, which is all an RHS or a stepper is, runs identically. Octave older than 7 has no `jsonencode`, so write JSON with `fprintf`.
3. **Handoff.** Write the script so that it needs no arguments beyond what is in the call, takes its inputs from files beside it, and writes its outputs into the model directory. Give the user the one line to paste into the MATLAB prompt and the names of the files that will appear, then wait. Check the files when they arrive (header, row count, finite values) before using them.

Say which rung was used in `NOTES.md`. A result from Octave is a result; a file nobody ran is not.

### MATLAB as the source

Translate the model to `.ode` as `translation.md` describes. The two files `fidelity` needs from the original code come from an export script the agent writes around the user's RHS function, untouched: it evaluates the RHS at 50 or more states spread over a beat (from the user's own `ode15s` run), writes `rhs_samples.csv` with `%.17g`, and writes `params.json`. Run it by the ladder. That is the only time the source is executed.

### MATLAB as the target

After the blocks are frozen:

1. `bimgrl codegen <model_dir> --backend julia --scheme BIRL-M`. Julia is the closest generated language to MATLAB (1-based, `end`-delimited, element-wise operators).
2. Translate three functions to `.m`: the scheme (`block_rush_larsen`), and the two initialisers. This is a transliteration, line for line. The constructs that are not the same:
   - Julia mutates `values` in place; MATLAB returns it: `function values = birlm_step(states, t, dt, parameters)`, starting with `values = zeros(size(states));`.
   - Every assignment needs a trailing `;`, or MATLAB prints a hundred thousand values per beat.
   - Indexing `x[i]` becomes `x(i)`; comments `#` become `%`; `for k in a:b` becomes `for k = a:b`; `maximum` becomes `max`; `Int(...)` is dropped (the sub-step count is already whole after `ceil`).
   - MATLAB has no ternary. `c ? a : b` becomes a call to a three-line helper `birlm_tern(c, a, b)` that returns `a` or `b`; both branches are evaluated, which is harmless here (a division by zero in the branch not taken yields `Inf`, not an error).
   - Identifiers may not start with an underscore: the generated temporaries (`_b0_n`, `_g_3`) need a prefix, the same one everywhere.
   - One statement can span several lines: the trigger ratio of a `binary` sub-step block is a `max` of absolute values that the generator expands into nested ternaries, tens of kilobytes long. Carry it across whole; do not rewrite it from the blocks file.
3. Write `birlm_trace.m`, the driver: the same loop as `drive.py` above (protocol B, `t` restarting every beat, the last beat on the 0.1 ms grid, the state at every beat start), the state names in the generated order as the header, every number with `%.17g`, and the run meta with `"language": "matlab"`. Run it by the ladder at the step and scheme that were frozen.
4. `bimgrl validate <model_dir> --test <trace>.csv`.

**The proof.** The translated stepper is the same arithmetic as the generated one, so at the same scheme and step its biomarkers must be the same numbers, not merely pass the same thresholds. Compare `vpeak`, `apd90`, `cat` and `tact` in `bimgrl/validation_<name>.json` with the columns of the same names in the `runs.csv` row for that scheme and step (on the no-compiler path, with the validation of the Python stepper). They agree to 1e-9 relative or better; on TT06 the Octave run of a translated BIMGRL stepper agreed with compiled C to 1e-12. Compare the values, not the `err_*` columns: `validate` interpolates the stored reference trace and the sweep evaluates the dense reference, so the errors differ in the second digit while the values do not. A larger difference is a translation error (a dropped term, a misplaced parenthesis, an off-by-one index), never a numerics result, and passing the thresholds does not excuse it.

**Until that comparison exists, the stepper is unverified.** If the handoff is refused or the trace never arrives, hand the `.m` files back with that word in the first sentence and in `NOTES.md`, and with the one line that would verify them. Never describe a stepper that was not run as validated, and never report a speed-up for it.

Expect pure MATLAB to be slow (seconds per beat for a 19-state model in Octave). That is the cost of the language, not of the scheme. When the user needs speed inside MATLAB and has a compiler configured (`mex -setup`), the generated C stepper behind a MEX gateway is the alternative to offer; it involves no hand translation.
