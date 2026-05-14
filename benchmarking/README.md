# BEAST Benchmarking

This folder contains the benchmarking infrastructure for measuring and validating BEAST performance.

## Folder layout

```
Benchmarking/
├── benchmark.py          # The benchmark runner script
├── README.md             # This file
├── Benchmarks/           # Benchmark XML files, organised by type
│   ├── Benchmark-big/
│   │   └── benchmark-big-1-min.xml
│   ├── Benchmark-small/
│   │   └── benchmark-small-1-min.xml
│   ├── Benchmark-tall/
│   ├── Benchmark-wide/
│   └── Benchmark-medium/
└── Results/              # Created automatically; one sub-folder per run
    └── {machine}_{type}_{number}_{model}_{timestamp}/
        ├── {prefix}benchmark-{type}-{number}-{model}.stdout
        ├── benchmark-{type}-{number}-{model}.log
        ├── benchmark-{type}-{number}-{model}.trees
        ├── {run_name}_report.md
        └── {run_name}_results.json
```

## Benchmark XML format

Each XML file begins with a structured header comment that records the dataset properties, expected MCMC settings, and the expected log-probability values for a canonical run at a fixed random seed.  The script parses this header to validate results.

```xml
<!-->
Benchmark: benchmark-small-1
Source data: YFV.fasta
Taxa: 71
Sites: 654
Site patterns: 280
Model: JC69, strict clock, coalescent constant population tree prior
Priors: rate ~ ApproximateReferencePrior, popSize ~ Gamma(0.001, 1000)
Operators:
  - SubtreeLeap on tree with weight 71.0
  ...
MCMC length: 10000000
Logging every: 10000
Seed: 1234
Initial Joint density: -29529.1683
Final Joint density: -6981.7346
<-->
```

The **Initial Joint density** is the log-probability of the initial state drawn at the given seed. The **Final Joint density** is the log-probability at the end of the chain. For a given seed, both values are fully deterministic and must match exactly across platforms if the implementation is correct. The benchmark runner reports `PASS` or `FAIL` for each.

### XML file naming convention

```
benchmark-{type}-{number}-{model}.xml
```

| Field    | Description |
| --- | --- |
| `type`   | Dataset category: `big`, `small`, `tall`, `wide`, `medium` |
| `number` | Dataset index within a type (1, 2, …) |
| `model`  | Model label — e.g. `min` for the minimum/baseline model configuration |

## Running a benchmark

```bash
python3 benchmark.py <type> <number> <model> [options]
```

**Example — run the small-1 minimum-model benchmark:**

```bash
python3 Benchmarking/benchmark.py small 1 min
```

**Example — tag the run with a machine name:**

```bash
python3 Benchmarking/benchmark.py small 1 min --machine bivalve
```

**Full usage:**

```
usage: benchmark.py [-h] [--seed SEED] [--beast BEAST] [--states STATES]
                    [--prefix PREFIX] [--benchmarks DIR] [--results DIR]
                    [--machine NAME] [--no-overwrite]
                    type number model

positional arguments:
  type              Benchmark type: big, small, tall, wide, medium
  number            Benchmark number (e.g. 1)
  model             Model label used in XML filename (e.g. min)

optional arguments:
  --seed SEED       Random seed (default: value from XML header, or 1234)
  --beast BEAST     Path to beast executable (default: 'beast' on PATH)
  --states STATES   Override MCMC chain length (for quick test runs)
  --prefix PREFIX   Override the run directory name
                    (default: {machine_}{type}_{number}_{model}_{timestamp})
  --benchmarks DIR  Path to Benchmarks folder
                    (default: Benchmarking/Benchmarks)
  --results DIR     Path to results parent folder
                    (default: Benchmarking/Results)
  --machine NAME    Machine identifier — prefixed to the run folder and
                    output file names
  --no-overwrite    Do not pass -overwrite to BEAST
```

## Output files

Each run creates a sub-folder inside `Results/`.  The folder name is:

```
{machine_}{type}_{number}_{model}_{timestamp}
```

where `{machine_}` is omitted if `--machine` is not supplied.

| File | Description |
| --- | --- |
| `{prefix}benchmark-{type}-{number}-{model}.stdout` | Full captured stdout from BEAST |
| `benchmark-{type}-{number}-{model}.log` | BEAST MCMC log (written by BEAST) |
| `benchmark-{type}-{number}-{model}.trees` | Tree log (written by BEAST) |
| `{run_name}_report.md` | Human-readable markdown report (see below) |
| `{run_name}_results.json` | Machine-readable JSON summary (see below) |

### Markdown report sections

- **Header table** — date, run directory, full command line
- **Dataset** — benchmark name, type, number, model, taxa, sites, site patterns
- **MCMC Settings** — chain length, logging frequency, seed
- **System** — macOS version, platform string, CPU, RAM, core counts, Python version
- **Software** — BEAST version, BEAGLE version, BEAGLE resource, BEAGLE thread count
- **Results / Joint Density** — table with observed vs expected initial and final joint log-density, deviation, and `✓ PASS` / `✗ FAIL` status
- **Results / Performance** — BEAST-reported runtime, wall-clock time, sample rate (min / million states)
- **Operator Analysis** — verbatim operator acceptance rate block from BEAST stdout

### JSON results file

```json
{
  "benchmark":        "benchmark-small-1-min",
  "type":             "small",
  "number":           1,
  "model":            "min",
  "taxa":             71,
  "sites":            654,
  "patterns":         280,
  "chain_length":     10000000,
  "sample_frequency": 10000,
  "seed":             1234,
  "timestamp":        "2026-05-14 12:00:00",
  "platform":         "macOS-26.3.1-arm64-arm-64bit",
  "cpu":              "Apple M3 Max",
  "ram_gb":           36,
  "cores_physical":   14,
  "cores_logical":    14,
  "beast_version":    "v10.5.0",
  "beagle_version":   "v4.0.1 (PRE-RELEASE)",
  "beagle_resource":  "CPU (arm64)",
  "beagle_threads":   14,
  "initial_joint":    -29529.1683,
  "final_joint":      -6981.7346,
  "runtime":          3.55,
  "wall_clock_time":  213.0,
  "sample_rate":      0.355
}
```

`runtime` is the BEAST-reported execution time in minutes. `wall_clock_time` is the wall-clock elapsed time in seconds. `sample_rate` is minutes per million MCMC states.

## Validation logic

Both the initial and final joint log-densities are **fully deterministic** for a given XML and seed. The runner checks both against the expected values recorded in the XML header:

- Deviation < 0.01 → **PASS**
- Deviation ≥ 0.01 → **FAIL** (exit code 2)

A failure indicates that either the BEAST implementation has changed in a way that alters the likelihood or MCMC trajectory, or that the benchmark XML header values need updating.

## Adding a new benchmark

1. Create a folder `Benchmarks/Benchmark-{type}/`.
2. Place the XML file as `benchmark-{type}-{number}-{model}.xml`.
3. Add the structured header comment (see format above).
4. Run once with the canonical seed to obtain the initial and final joint densities:
   ```bash
   python3 Benchmarking/benchmark.py {type} {number} {model} --seed 1234
   ```
5. Copy the reported joint densities into the XML header.
6. Re-run to confirm both values show `✓ PASS`.
