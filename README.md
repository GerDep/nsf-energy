# NSF energy diagnostics

Python code for the retained connected graviton channels `2 -> 2` and
`2 -> 3`, their conditional mean outgoing energy, and the physical variance
of that energy. The units are `Ei = G = epsilon = 1`.

The calculation includes the exchange iteration `-T1^2/2` and the cubic
iteration `-i T1^3/6`. It does **not** include an independently evaluated
quintic contact `i T3`. Conditional cancellations of quartic `T2` sectors
are not a proof that all omitted physical amplitudes vanish. See
[docs/NUMERICS.md](docs/NUMERICS.md) before interpreting the results.

## Install

Use Python 3.10 or newer. Python 3.12 is used by the configured CI.
In Windows PowerShell, open this folder and run:

```powershell
py -3 -m venv .venv
.\.venv\Scripts\python.exe -m pip install -r requirements.txt
```

All commands below use `python`; on Windows you can replace it with
`.\.venv\Scripts\python.exe` without activating the environment.
On Linux/macOS, create a virtual environment with `python3 -m venv .venv`.
For the exact direct-dependency versions used in the local Python 3.12
checks, install `requirements-tested.txt` instead of `requirements.txt`.

## Quick example: tables without a new integration

```powershell
python analysis/render_tables.py
```

This uses the checked summary in `results/reference/variance.json` and
creates `results/generated/varianza_NSF_formato_nota.tex` and a plain-text
table. With a normal LaTeX installation, compile the generated `.tex` to PDF.
The standalone checked LaTeX is also in `docs/varianza_NSF_formato_nota.tex`.

To redraw the stored mean-energy scan without integrating:

```powershell
python graficos_nsf_local.py --solo-graficos
```

## Calculate a new point

```powershell
python graficos_nsf_local.py --ir 1e-8 --emax 2024 --workers 4
```

`--ir` means `kappa_IR / Ei`; `--emax` means the **total** final energy
bound `Emax / Ei`, not a bound per graviton. `--help` lists the controls.
Known points are reused. Add `--recalcular` to integrate them again.
The pure Python calculation can take hours. Samples are checkpointed and
the same command resumes an interrupted run. The program updates its
history, CSV tables and mean-energy figures in `resultados_nsf/`.

## Scan the upper-energy bound

```powershell
python barrido_energia_maxima.py --workers 4
```

This separate driver fixes `kappa_IR / Ei = 1e-8` and uses
`Emax / Ei = 2024, 4048, 16192, 18096`. The last two values follow the
original request literally: `4*4048` and `2*9048`.
For different choices use `--maximos 2024 4048 8096 16192`.
It preserves the original history and writes to
`resultados_barrido_emax_IR_1e-8/`.

The standalone drivers retain their original embedded references. Their
rounded prior value at 4048 is different from the newer paired-tail estimate
in the variance reference dataset. They are separate runs, with distinct
covariance information; they must not be silently merged.

## Reproduce the variance from original weighted samples

Download the separate `nsf_reference_data.zip` asset supplied with this
repository or attached to its GitHub release. It contains numerical data,
not manuscripts or compiled libraries. Keep it outside the source upload.

```powershell
python scripts/import_reference_data.py C:\ruta\nsf_reference_data.zip
python analysis/reorganize_variance.py
python analysis/render_tables.py --input results/generated/variance.json
```

The importer verifies the manifest and SHA-256 hashes before installing
the samples in `reference_data/`. The reduction rebuilds all moments and
covariances without new amplitude evaluations. The data contain four IR
cutoffs (`1e-5` through `1e-8`) at `Emax/Ei=2024`, the paired upper-bound
scan, and the angular-partition comparison. Full input precision is retained.

The reprocessor expects this archived layout. It is not an automatic
postprocessor for every future local-driver history.

## Checks

```powershell
python -m unittest discover -s tests -v
```

The tests check folded moments against explicit energy distributions,
nonzero BCH reconstruction, coherent interference and missing-input guards.
They do not establish the physical completeness of the retained amplitudes.
GitHub Actions is configured to run the tests on Linux and Windows.

## Files

| Path | Purpose |
|---|---|
| `graficos_nsf_local.py` | Original self-contained integrator and cumulative mean plots |
| `barrido_energia_maxima.py` | Separate upper-energy scan |
| `analysis/` | Archived-sample variance reduction and report generation |
| `algebra/` | BCH reconstruction and coherent-weight interfaces |
| `results/reference/` | Small, frozen numerical summary for quick reproduction |
| `docs/NUMERICS.md` | Measures, amplitudes, regulators and error definitions |
| `SUBIR_A_GITHUB.md` | Upload instructions in Spanish |

Generated samples, local histories, environments and binaries are excluded
by `.gitignore`. No software license has been selected by the authors in
this prepared package; choose the intended reuse terms before distributing
it as a licensed project.
