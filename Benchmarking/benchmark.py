#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
benchmark.py — Run a BEAST benchmark XML and report initial/final Joint density.

Usage:
    python benchmark.py <type> <number> [--seed SEED] [--beast PATH] [--states STATES]

Arguments:
    type     Benchmark type: big, small, tall, wide, medium
    number   Benchmark number (e.g. 1 or 2)

Options:
    --seed SEED       Random seed passed to BEAST (default: value from XML header, or 1234)
    --beast PATH      Path to beast executable or jar (default: 'beast' on PATH)
    --states STATES   Override number of MCMC states from XML header (for shorter test runs)
    --prefix PREFIX   Override the output file prefix (default: auto timestamp)
    --no-overwrite    Do not pass -overwrite to BEAST

The XML header comment (see benchmark-big-1.xml) is parsed for expected values:
    MCMC length, Logging every, Seed, Initial Joint density, Final Joint density
"""

import argparse
import csv
import json
import os
import platform
import re
import shutil
import subprocess
import sys
from datetime import datetime
from pathlib import Path
from typing import Optional


# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------

BENCHMARKS_DIR = Path("Benchmarks")
RESULTS_DIR    = Path("Results")
CSV_PATH       = Path("results.csv")


# ---------------------------------------------------------------------------
# XML header parsing
# ---------------------------------------------------------------------------

# Pattern for the structured header block: <!> ... <-->
# e.g.  <!>
#       Benchmark: benchmark-big-1
#       ...
#       <-->
HEADER_BLOCK_RE = re.compile(r'<!-->\s*(.*?)\s*<-->', re.DOTALL)

# Individual field patterns within the header block
FIELD_PATTERNS = {
    "benchmark":     re.compile(r'^Benchmark:\s*(.+)$',                  re.MULTILINE),
    "taxa":          re.compile(r'^Taxa:\s*(\d+)$',                       re.MULTILINE),
    "sites":         re.compile(r'^Sites:\s*(\d+)$',                      re.MULTILINE),
    "patterns":      re.compile(r'^Site patterns:\s*(\d+)$',              re.MULTILINE),
    "mcmc_length":   re.compile(r'^MCMC length:\s*([\d,]+)',              re.MULTILINE),
    "log_every":     re.compile(r'^Logging every:\s*([\d,]+)',            re.MULTILINE),
    "seed":          re.compile(r'^Seed:\s*(\d+)$',                       re.MULTILINE),
    "initial_joint": re.compile(r'^Initial Joint density:\s*([-\d.]+)$', re.MULTILINE),
    "final_joint":   re.compile(r'^Final Joint density:\s*([-\d.]+)$',   re.MULTILINE),
}


def parse_xml_header(xml_path: Path) -> dict:
    """Return a dict of metadata fields parsed from the structured XML header comment."""
    try:
        with open(xml_path, encoding="utf-8") as fh:
            # Only need the first ~2000 chars
            head = fh.read(2000)
    except OSError as exc:
        sys.exit(f"Error reading XML: {exc}")

    block_match = HEADER_BLOCK_RE.search(head)
    if not block_match:
        return {}

    block = block_match.group(1)
    result = {}
    for key, pat in FIELD_PATTERNS.items():
        m = pat.search(block)
        if m:
            val = m.group(1).replace(",", "")
            result[key] = int(val) if val.lstrip("-").isdigit() else (
                float(val) if re.fullmatch(r'-?\d+\.\d+', val) else val
            )
    return result


# ---------------------------------------------------------------------------
# stdout parsing
# ---------------------------------------------------------------------------

# Column header line looks like:
#   state\tJoint       \tPrior       \t...
STATE_HEADER_RE = re.compile(r'^state\t')

# Data rows: first token is an integer state number
STATE_ROW_RE = re.compile(r'^(\d+)\t')

# Timing line at the very end of BEAST output: "X.Y minutes"
BEAST_MINUTES_RE = re.compile(r'^([\d.]+)\s+minutes\s*$')

# "X minutes/million states" embedded in state rows (last field)
RATE_RE = re.compile(r'([\d.]+)\s+minutes/million states')

# BEAST / BEAGLE version strings
BEAST_VERSION_RE   = re.compile(r'BEAST (v[\d.]+)')
BEAGLE_LIB_RE      = re.compile(r'Using BEAGLE library (v\S+(?:\s+\([^)]+\))?)\s+for')
BEAGLE_RESOURCE_RE = re.compile(r'Using BEAGLE resource \d+:\s*(.+)')
BEAGLE_THREADS_RE  = re.compile(r'Using (\d+) threads for')


def parse_versions(text: str) -> dict:
    """Extract BEAST and BEAGLE version/resource info from BEAST stdout."""
    v: dict = {}
    m = BEAST_VERSION_RE.search(text)
    if m:
        v["beast"] = m.group(1)
    m = BEAGLE_LIB_RE.search(text)
    if m:
        v["beagle"] = m.group(1).strip()
    m = BEAGLE_RESOURCE_RE.search(text)
    if m:
        v["beagle_resource"] = m.group(1).strip()
    m = BEAGLE_THREADS_RE.search(text)
    if m:
        v["beagle_threads"] = int(m.group(1))
    return v


def parse_stdout(text: str) -> dict:
    """
    Parse BEAST stdout.  Returns:
        columns        list of column names (from header)
        joint_col      index of the Joint column
        initial_row    dict of column->value for state 0
        final_row      dict of column->value for the last state
        beast_minutes  total runtime in minutes (float) as reported by BEAST
        final_rate     minutes/million states at the final state (float or None)
    """
    lines = text.splitlines()

    columns = []
    joint_col = None
    rows = []
    beast_minutes = None

    for line in lines:
        # Locate column header
        if not columns and STATE_HEADER_RE.match(line):
            raw_cols = line.rstrip().split('\t')
            columns = [c.strip() for c in raw_cols]
            try:
                joint_col = columns.index('Joint')
            except ValueError:
                joint_col = 1  # fallback: second column
            continue

        # Data rows (only after header is found)
        if columns and STATE_ROW_RE.match(line):
            parts = line.rstrip().split('\t')
            row = {}
            for i, col in enumerate(columns):
                if i < len(parts):
                    row[col] = parts[i].strip()
            rows.append(row)
            continue

        # "X.Y minutes" total runtime line
        m = BEAST_MINUTES_RE.match(line.strip())
        if m:
            beast_minutes = float(m.group(1))

    initial_row = rows[0] if rows else {}
    final_row = rows[-1] if rows else {}

    # Extract rate from the last data row's trailing field
    final_rate = None
    if rows:
        last_line_content = '\t'.join(rows[-1].get(c, '') for c in columns)
        rm = RATE_RE.search(last_line_content)
        if rm:
            final_rate = float(rm.group(1))
        # also check the raw lines
        for line in reversed(lines):
            rm = RATE_RE.search(line)
            if rm:
                final_rate = float(rm.group(1))
                break

    return {
        "columns":       columns,
        "joint_col":     joint_col,
        "initial_row":   initial_row,
        "final_row":     final_row,
        "beast_minutes": beast_minutes,
        "final_rate":    final_rate,
    }


# ---------------------------------------------------------------------------
# System information
# ---------------------------------------------------------------------------

def get_system_info() -> dict:
    """Collect hardware and OS information (cross-platform)."""
    import os as _os
    info: dict = {
        "os":      platform.platform(),
        "machine": platform.machine(),
        "python":  platform.python_version(),
        "cpu":     platform.processor() or platform.machine() or "unknown",
    }

    sys_name = platform.system()  # 'Darwin', 'Linux', 'Windows'

    # ---- macOS via sysctl ------------------------------------------------
    if sys_name == "Darwin":
        def _sysctl(key: str) -> Optional[str]:
            try:
                return subprocess.check_output(
                    ["sysctl", "-n", key], stderr=subprocess.DEVNULL
                ).decode().strip()
            except Exception:
                return None

        cpu = _sysctl("machdep.cpu.brand_string")
        if cpu:
            info["cpu"] = cpu

        mem = _sysctl("hw.memsize")
        if mem:
            info["ram_gb"] = int(mem) // (1024 ** 3)

        phys = _sysctl("hw.physicalcpu")
        if phys:
            info["cores_physical"] = int(phys)

        logical = _sysctl("hw.logicalcpu")
        if logical:
            info["cores_logical"] = int(logical)

        try:
            sw = subprocess.check_output(
                ["sw_vers"], stderr=subprocess.DEVNULL
            ).decode().strip()
            parts: dict = {}
            for ln in sw.splitlines():
                if ":" in ln:
                    k, _, v = ln.partition(":")
                    parts[k.strip()] = v.strip()
            name  = parts.get("ProductName", "")
            ver   = parts.get("ProductVersion", "")
            build = parts.get("BuildVersion", "")
            info["macos"] = f"{name} {ver} ({build})".strip()
        except Exception:
            pass

    # ---- Linux via /proc -------------------------------------------------
    elif sys_name == "Linux":
        # CPU brand
        try:
            with open("/proc/cpuinfo", encoding="utf-8") as fh:
                for line in fh:
                    if line.startswith("model name"):
                        info["cpu"] = line.split(":", 1)[1].strip()
                        break
        except Exception:
            pass

        # RAM
        try:
            with open("/proc/meminfo", encoding="utf-8") as fh:
                for line in fh:
                    if line.startswith("MemTotal"):
                        kb = int(line.split()[1])
                        info["ram_gb"] = round(kb / (1024 ** 2), 1)
                        break
        except Exception:
            pass

        # Cores — logical from os.cpu_count(), physical via /proc/cpuinfo
        info["cores_logical"] = _os.cpu_count() or 1
        try:
            cores = set()
            with open("/proc/cpuinfo", encoding="utf-8") as fh:
                for line in fh:
                    if line.startswith("core id"):
                        cores.add(line.split(":", 1)[1].strip())
            if cores:
                info["cores_physical"] = len(cores)
        except Exception:
            info["cores_physical"] = info["cores_logical"]

    # ---- Windows via wmi / winreg ----------------------------------------
    elif sys_name == "Windows":
        import ctypes

        # CPU brand from registry
        try:
            import winreg
            key = winreg.OpenKey(
                winreg.HKEY_LOCAL_MACHINE,
                r"HARDWARE\DESCRIPTION\System\CentralProcessor\0",
            )
            info["cpu"] = winreg.QueryValueEx(key, "ProcessorNameString")[0].strip()
            winreg.CloseKey(key)
        except Exception:
            pass

        # RAM via GlobalMemoryStatusEx
        try:
            class _MemStatus(ctypes.Structure):
                _fields_ = [
                    ("dwLength",                ctypes.c_ulong),
                    ("dwMemoryLoad",            ctypes.c_ulong),
                    ("ullTotalPhys",            ctypes.c_ulonglong),
                    ("ullAvailPhys",            ctypes.c_ulonglong),
                    ("ullTotalPageFile",        ctypes.c_ulonglong),
                    ("ullAvailPageFile",        ctypes.c_ulonglong),
                    ("ullTotalVirtual",         ctypes.c_ulonglong),
                    ("ullAvailVirtual",         ctypes.c_ulonglong),
                    ("ullAvailExtendedVirtual", ctypes.c_ulonglong),
                ]
            ms = _MemStatus()
            ms.dwLength = ctypes.sizeof(ms)
            ctypes.windll.kernel32.GlobalMemoryStatusEx(ctypes.byref(ms))
            info["ram_gb"] = round(ms.ullTotalPhys / (1024 ** 3), 1)
        except Exception:
            pass

        # Cores via WMI (optional but preferred)
        try:
            import subprocess as _sp
            out = _sp.check_output(
                ["wmic", "cpu", "get", "NumberOfCores,NumberOfLogicalProcessors", "/format:csv"],
                stderr=_sp.DEVNULL,
            ).decode()
            rows = [r.strip() for r in out.splitlines() if r.strip() and r.strip().lower() != "node,numberofcores,numberoflogicalprocessors"]
            if rows:
                _, phys_s, logical_s = rows[0].split(",")
                info["cores_physical"] = int(phys_s)
                info["cores_logical"]  = int(logical_s)
        except Exception:
            info["cores_logical"] = _os.cpu_count() or 1

    # ---- psutil fallback for any missing fields --------------------------
    try:
        import psutil
        if "ram_gb" not in info:
            info["ram_gb"] = round(psutil.virtual_memory().total / (1024 ** 3), 1)
        if "cores_logical" not in info:
            info["cores_logical"] = psutil.cpu_count(logical=True) or 1
        if "cores_physical" not in info:
            info["cores_physical"] = psutil.cpu_count(logical=False) or info["cores_logical"]
    except ImportError:
        pass

    # ---- GPU detection ---------------------------------------------------
    import json as _json
    gpu_names: list = []
    gpu_cores_list: list = []

    if sys_name == "Darwin":
        try:
            sp_raw = subprocess.check_output(
                ["system_profiler", "SPDisplaysDataType", "-json"],
                stderr=subprocess.DEVNULL,
            ).decode()
            for entry in _json.loads(sp_raw).get("SPDisplaysDataType", []):
                name = entry.get("sppci_model") or entry.get("_name", "")
                if name:
                    gpu_names.append(name)
                cores_str = entry.get("sppci_cores")
                if cores_str:
                    try:
                        gpu_cores_list.append(int(str(cores_str).replace(",", "")))
                    except ValueError:
                        pass
        except Exception:
            pass

    elif sys_name == "Linux":
        # NVIDIA via nvidia-smi
        try:
            nv_raw = subprocess.check_output(
                ["nvidia-smi", "--query-gpu=name", "--format=csv,noheader"],
                stderr=subprocess.DEVNULL,
            ).decode().strip()
            if nv_raw:
                gpu_names.extend(n.strip() for n in nv_raw.splitlines() if n.strip())
        except Exception:
            pass
        # lspci fallback
        if not gpu_names:
            try:
                lspci_raw = subprocess.check_output(
                    ["lspci"], stderr=subprocess.DEVNULL
                ).decode()
                for line in lspci_raw.splitlines():
                    if any(t in line.lower() for t in ("vga", "3d controller", "display controller")):
                        gpu_names.append(line.split(":", 2)[-1].strip())
            except Exception:
                pass

    elif sys_name == "Windows":
        try:
            wmic_raw = subprocess.check_output(
                ["wmic", "path", "win32_VideoController", "get", "Name", "/format:csv"],
                stderr=subprocess.DEVNULL,
            ).decode()
            for row in wmic_raw.splitlines():
                row = row.strip()
                if row and "," in row:
                    name = row.split(",", 1)[1].strip()
                    if name and name.lower() != "name":
                        gpu_names.append(name)
        except Exception:
            pass
        # nvidia-smi provides cleaner names than wmic
        try:
            nv_raw = subprocess.check_output(
                ["nvidia-smi", "--query-gpu=name", "--format=csv,noheader"],
                stderr=subprocess.DEVNULL,
            ).decode().strip()
            if nv_raw:
                gpu_names = [n.strip() for n in nv_raw.splitlines() if n.strip()]
        except Exception:
            pass

    # pynvml: CUDA core count for NVIDIA GPUs (any platform, optional)
    if not gpu_cores_list:
        try:
            import pynvml
            pynvml.nvmlInit()
            for i in range(pynvml.nvmlDeviceGetCount()):
                handle = pynvml.nvmlDeviceGetHandleByIndex(i)
                gpu_cores_list.append(pynvml.nvmlDeviceGetNumGpuCores(handle))
            pynvml.nvmlShutdown()
        except Exception:
            pass

    if gpu_names:
        info["gpu"] = "; ".join(gpu_names)
    if gpu_cores_list:
        info["gpu_cores"] = sum(gpu_cores_list)

    return info


# ---------------------------------------------------------------------------
# Progress bar
# ---------------------------------------------------------------------------

def print_progress(state: int, total_states: Optional[int],
                   rate_str: Optional[str], elapsed_s: float) -> None:
    """Overwrite the current terminal line with a progress bar."""
    term_width = shutil.get_terminal_size((80, 20)).columns
    bar_width = 30
    elapsed_str = f"{elapsed_s / 60:.1f} min"
    rate_part = f"  {rate_str}" if rate_str else ""

    if total_states:
        frac = min(state / total_states, 1.0)
        filled = int(bar_width * frac)
        bar = '\u2588' * filled + '\u2591' * (bar_width - filled)
        pct = f"{100 * frac:.1f}%"
        state_str = f"{state:,}/{total_states:,}"
    else:
        bar = '\u2591' * bar_width
        pct = "   ?%"
        state_str = f"{state:,}"

    line = f"\r  [{bar}] {pct}  state {state_str}  elapsed {elapsed_str}{rate_part}"
    # Pad to clear any previous longer content
    padded = line + ' ' * max(0, term_width - len(line) - 1)
    print(padded, end='', flush=True)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def format_deviation(actual: float, expected: float) -> str:
    diff = actual - expected
    sign = "+" if diff >= 0 else ""
    return f"{sign}{diff:.4f}  ({sign}{100 * diff / abs(expected):.2f}%)"


def get_joint(row: dict) -> Optional[float]:
    val = row.get("Joint")
    if val is None:
        return None
    try:
        return float(val)
    except ValueError:
        return None


# ---------------------------------------------------------------------------
# Report generation
# ---------------------------------------------------------------------------

def _md_table(headers: list, rows: list) -> list:
    """Return a list of markdown table lines."""
    sep = "|" + "|".join(" --- " for _ in headers) + "|"
    out = ["| " + " | ".join(str(h) for h in headers) + " |", sep]
    for row in rows:
        out.append("| " + " | ".join(str(c) for c in row) + " |")
    return out


def write_report(report_path: Path, run_name: str, btype: str, bnum: int,
                 model: str, seed: int, timestamp: str, meta: dict, sys_info: dict,
                 versions: dict, parsed: dict, elapsed: float,
                 stdout_text: str, cmd: list) -> None:
    """Write a markdown benchmark report."""
    initial_joint = get_joint(parsed["initial_row"])
    final_joint   = get_joint(parsed["final_row"])
    final_state   = parsed["final_row"].get("state", "?")
    try:
        final_state_fmt = f"{int(final_state):,}"
    except (ValueError, TypeError):
        final_state_fmt = str(final_state)

    out = [
        f"# Benchmark Report: benchmark-{btype}-{bnum}-{model}",
        "",
        "| | |",
        "| --- | --- |",
        f"| **Date** | {timestamp} |",
        f"| **Run directory** | `{run_name}` |",
        f"| **Command** | `{' '.join(cmd)}` |",
        "",
    ]

    # Dataset
    out += ["## Dataset", ""]
    ds_rows: list = [("Benchmark", f"benchmark-{btype}-{bnum}-{model}"),
                     ("Type", btype), ("Number", str(bnum)), ("Model", model)]
    if "taxa"     in meta: ds_rows.append(("Taxa",          f"{meta['taxa']:,}"))
    if "sites"    in meta: ds_rows.append(("Sites",         f"{meta['sites']:,}"))
    if "patterns" in meta: ds_rows.append(("Site patterns", f"{meta['patterns']:,}"))
    out += _md_table(["Property", "Value"], ds_rows)

    # MCMC settings
    out += ["", "## MCMC Settings", ""]
    mcmc_rows: list = []
    if "mcmc_length" in meta: mcmc_rows.append(("Chain length", f"{meta['mcmc_length']:,}"))
    if "log_every"   in meta: mcmc_rows.append(("Log every",    f"{meta['log_every']:,}"))
    mcmc_rows.append(("Seed", str(seed)))
    out += _md_table(["Property", "Value"], mcmc_rows)

    # System
    out += ["", "## System", ""]
    sys_rows: list = []
    if "macos" in sys_info:
        sys_rows.append(("macOS", sys_info["macos"]))
    sys_rows.append(("Platform",       sys_info.get("os", "unknown")))
    sys_rows.append(("CPU",            sys_info.get("cpu", "unknown")))
    if "ram_gb"         in sys_info: sys_rows.append(("RAM",            f"{sys_info['ram_gb']} GB"))
    if "cores_physical" in sys_info: sys_rows.append(("Physical cores", str(sys_info["cores_physical"])))
    if "cores_logical"  in sys_info: sys_rows.append(("Logical cores",  str(sys_info["cores_logical"])))
    if "gpu"            in sys_info: sys_rows.append(("GPU",             sys_info["gpu"]))
    if "gpu_cores"      in sys_info: sys_rows.append(("GPU cores",       str(sys_info["gpu_cores"])))
    sys_rows.append(("Python", sys_info.get("python", "unknown")))
    out += _md_table(["Property", "Value"], sys_rows)

    # Software versions
    out += ["", "## Software", ""]
    sw_rows = [
        ("BEAST",           versions.get("beast",           "unknown")),
        ("BEAGLE",          versions.get("beagle",          "unknown")),
        ("BEAGLE resource", versions.get("beagle_resource", "unknown")),
        ("BEAGLE threads",  str(versions.get("beagle_threads", "?"))),
    ]
    out += _md_table(["Software", "Version / Details"], sw_rows)

    # Results — Joint densities
    out += ["", "## Results", "", "### Joint Density", ""]
    joint_rows: list = []
    if initial_joint is not None:
        exp = meta.get("initial_joint")
        if exp is not None:
            diff = initial_joint - exp
            dev  = f"`{diff:+.4f}` ({100 * diff / abs(exp):+.2f}%)"
            status = "✓ PASS" if abs(diff) < 0.01 else "✗ FAIL"
            joint_rows.append(("0 (initial)", f"`{initial_joint:.4f}`",
                                f"`{exp:.4f}`", dev, status))
        else:
            joint_rows.append(("0 (initial)", f"`{initial_joint:.4f}`", "—", "—", "—"))
    if final_joint is not None:
        exp = meta.get("final_joint")
        if exp is not None:
            diff   = final_joint - exp
            dev    = f"`{diff:+.4f}` ({100 * diff / abs(exp):+.2f}%)"
            status = "✓ PASS" if abs(diff) < 0.01 else "✗ FAIL"
            joint_rows.append((f"{final_state_fmt} (final)", f"`{final_joint:.4f}`",
                                f"`{exp:.4f}`", dev, status))
        else:
            joint_rows.append((f"{final_state_fmt} (final)",
                                f"`{final_joint:.4f}`", "—", "—", "—"))
    out += _md_table(["State", "Joint", "Expected", "Deviation", "Status"], joint_rows)

    # Performance
    out += ["", "### Performance", ""]
    perf_rows: list = []
    beast_mins = parsed.get("beast_minutes")
    if beast_mins is not None:
        perf_rows.append(("BEAST runtime",   f"{beast_mins:.2f} min ({beast_mins * 60:.0f} s)"))
    perf_rows.append(("Wall-clock time", f"{elapsed / 60:.2f} min ({elapsed:.0f} s)"))
    if parsed.get("final_rate") is not None:
        perf_rows.append(("Rate", f"{parsed['final_rate']:.3f} min / million states"))
    out += _md_table(["Metric", "Value"], perf_rows)

    # Operator analysis
    op_idx = stdout_text.find("Operator analysis")
    if op_idx != -1:
        op_block = stdout_text[op_idx:].split("\n\n")[0].rstrip()
        out += ["", "## Operator Analysis", "", "```", op_block, "```"]

    report_path.write_text('\n'.join(out) + '\n', encoding='utf-8')


def write_json(json_path: Path, btype: str, bnum: int, model: str, seed: int,
              timestamp: str, meta: dict, sys_info: dict, versions: dict,
              parsed: dict, elapsed: float) -> None:
    """Write a machine-readable JSON results file."""
    initial_joint = get_joint(parsed["initial_row"])
    final_joint   = get_joint(parsed["final_row"])
    data = {
        "benchmark":        f"benchmark-{btype}-{bnum}-{model}",
        "type":             btype,
        "number":           bnum,
        "model":            model,
        "taxa":             meta.get("taxa"),
        "sites":            meta.get("sites"),
        "patterns":         meta.get("patterns"),
        "chain_length":     meta.get("mcmc_length"),
        "sample_frequency": meta.get("log_every"),
        "seed":             seed,
        "timestamp":        timestamp,
        "platform":         sys_info.get("os"),
        "cpu":              sys_info.get("cpu"),
        "ram_gb":           sys_info.get("ram_gb"),
        "cores_physical":   sys_info.get("cores_physical"),
        "cores_logical":    sys_info.get("cores_logical"),
        "gpu":              sys_info.get("gpu"),
        "gpu_cores":        sys_info.get("gpu_cores"),
        "beast_version":    versions.get("beast"),
        "beagle_version":   versions.get("beagle"),
        "beagle_resource":  versions.get("beagle_resource"),
        "beagle_threads":   versions.get("beagle_threads"),
        "initial_joint":      initial_joint,
        "initial_deviation":  round(initial_joint - meta["initial_joint"], 6) if initial_joint is not None and meta.get("initial_joint") is not None else None,
        "final_joint":        final_joint,
        "final_deviation":    round(final_joint   - meta["final_joint"],   6) if final_joint   is not None and meta.get("final_joint")   is not None else None,
        "runtime":          parsed.get("beast_minutes"),
        "wall_clock_time":  round(elapsed, 3),
        "sample_rate":      parsed.get("final_rate"),
    }
    json_path.write_text(json.dumps(data, indent=2), encoding='utf-8')
    return data


def append_csv(csv_path: Path, data: dict) -> None:
    """Append a results row to the CSV, creating it with a header row if absent."""
    headers = list(data.keys())
    values  = [data[h] for h in headers]

    if csv_path.exists():
        with open(csv_path, encoding='utf-8', newline='') as fh:
            existing_headers = next(csv.reader(fh), None)
        if existing_headers != headers:
            print(f"\n  WARNING: Cannot append to {csv_path.name} — "
                  f"column headings differ from expected.")
            print(f"  Expected: {headers}")
            print(f"  Found:    {existing_headers}")
            return
        with open(csv_path, 'a', encoding='utf-8', newline='') as fh:
            csv.writer(fh).writerow(values)
    else:
        with open(csv_path, 'w', encoding='utf-8', newline='') as fh:
            writer = csv.writer(fh)
            writer.writerow(headers)
            writer.writerow(values)

    print(f"  Appended to:  {csv_path.name}")

def main():
    parser = argparse.ArgumentParser(
        description="Run a BEAST benchmark and report Joint density results.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("type",   help="Benchmark type: big, small, tall, wide, medium")
    parser.add_argument("number", type=int, help="Benchmark number (e.g. 1)")
    parser.add_argument("model",  help="Model label used in XML filename (e.g. min)")
    parser.add_argument("--seed",   type=int, default=None,
                        help="Random seed (default: from XML header or 1234)")
    parser.add_argument("--beast",  default="beast",
                        help="Path to beast executable (default: 'beast' on PATH)")
    parser.add_argument("--states", type=int, default=None,
                        help="Override MCMC length (for quick test runs)")
    parser.add_argument("--prefix", default=None,
                        help="Override the run directory name (default: {type}_{number}_{model}_{timestamp})")
    parser.add_argument("--benchmarks", default=None, metavar="DIR",
                        help="Path to Benchmarks folder (default: Benchmarking/Benchmarks)")
    parser.add_argument("--results", default=None, metavar="DIR",
                        help="Path to results parent folder (default: Benchmarking/Results)")
    parser.add_argument("--machine", default=None, metavar="NAME",
                        help="Machine identifier — prefixed to the run folder and output file names")
    parser.add_argument("--no-overwrite", action="store_true",
                        help="Do not pass -overwrite to BEAST")
    parser.add_argument("--append", action="store_true",
                        help="Append results to results.csv in the current directory")
    args = parser.parse_args()

    btype = args.type.lower()
    bnum  = args.number
    model = args.model
    machine      = args.machine or ""
    file_prefix  = f"{machine}_" if machine else ""
    valid_types = {"big", "small", "tall", "wide", "medium"}
    if btype not in valid_types:
        sys.exit(f"Unknown type '{btype}'. Must be one of: {', '.join(sorted(valid_types))}")

    # Resolve configurable paths
    benchmarks_dir = Path(args.benchmarks) if args.benchmarks else BENCHMARKS_DIR
    results_dir    = Path(args.results)    if args.results    else RESULTS_DIR

    # Locate XML
    xml_path = benchmarks_dir / f"Benchmark-{btype}" / f"benchmark-{btype}-{bnum}-{model}.xml"
    if not xml_path.exists():
        sys.exit(f"XML not found: {xml_path}")

    # Parse header metadata
    meta = parse_xml_header(xml_path)
    if meta:
        print("  Benchmark metadata from XML header:")
        if "taxa"        in meta: print(f"    Taxa:            {meta['taxa']:,}")
        if "sites"       in meta: print(f"    Sites:           {meta['sites']:,}")
        if "patterns"    in meta: print(f"    Site patterns:   {meta['patterns']:,}")
        if "mcmc_length" in meta: print(f"    MCMC length:     {meta['mcmc_length']:,}")
        if "log_every"   in meta: print(f"    Logging every:   {meta['log_every']:,}")
    else:
        print("  No structured header comment found in XML — will report observed values only.")

    # Resolve seed
    seed = args.seed if args.seed is not None else meta.get("seed", 1234)
    print(f"  Seed:              {seed}")

    # Run directory: Results/{type}_{number}_{model}_{timestamp}/
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    run_name  = args.prefix if args.prefix else f"{file_prefix}{btype}_{bnum}_{model}_{timestamp}"
    results_dir.mkdir(parents=True, exist_ok=True)
    run_dir   = results_dir / run_name
    run_dir.mkdir(parents=True, exist_ok=True)

    # Collect system info before the run
    sys_info = get_system_info()

    # Build BEAST command (absolute XML path; cwd will be run_dir)
    cmd = [
        args.beast,
        "-seed", str(seed),
        "-citations_off",
        "-tests", "0",
    ]
    if not args.no_overwrite:
        cmd.append("-overwrite")
    cmd.append(str(xml_path))

    print(f"\n  Run directory: {run_dir.name}")
    print(f"  Command: {' '.join(cmd)}")
    print("-" * 70)
    print("  Running BEAST...\n")

    total_states = args.states if args.states else meta.get("mcmc_length")

    # Stream BEAST output line-by-line to drive the progress bar
    start_time   = datetime.now()
    stdout_lines: list = []
    header_found = False
    current_rate = None

    try:
        proc = subprocess.Popen(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            bufsize=1,
            cwd=run_dir,
        )
    except FileNotFoundError:
        sys.exit(f"BEAST executable not found: '{args.beast}'\n"
                 "Use --beast to specify the path.")

    for line in proc.stdout:
        stdout_lines.append(line)
        stripped = line.rstrip()

        if STATE_HEADER_RE.match(stripped):
            header_found = True
            continue

        if header_found and STATE_ROW_RE.match(stripped):
            m = STATE_ROW_RE.match(stripped)
            state_num = int(m.group(1))
            rm = RATE_RE.search(stripped)
            if rm:
                current_rate = rm.group(1) + " min/Mstates"
            elapsed_s = (datetime.now() - start_time).total_seconds()
            print_progress(state_num, total_states, current_rate, elapsed_s)

    proc.wait()
    print()  # newline after final progress bar line

    elapsed     = (datetime.now() - start_time).total_seconds()
    stdout_text = ''.join(stdout_lines)

    # Save stdout inside the run directory
    stdout_file = run_dir / f"{file_prefix}benchmark-{btype}-{bnum}-{model}.stdout"
    stdout_file.write_text(stdout_text, encoding="utf-8")
    print(f"  Stdout saved to: {run_dir.name}/{stdout_file.name}")

    if proc.returncode != 0:
        print(f"\n  WARNING: BEAST exited with code {proc.returncode}")
        tail = stdout_text.strip().splitlines()[-20:]
        print("  --- Last output ---")
        for line in tail:
            print(f"  {line}")
        sys.exit(1)

    # Parse output
    parsed   = parse_stdout(stdout_text)
    versions = parse_versions(stdout_text)

    print()
    print("=" * 70)
    print(f"  RESULTS: benchmark-{btype}-{bnum}-{model}  [seed={seed}]")
    print("=" * 70)

    if parsed["columns"]:
        print(f"  Log columns: {', '.join(parsed['columns'])}")

    initial_joint = get_joint(parsed["initial_row"])
    final_joint   = get_joint(parsed["final_row"])
    final_state   = parsed["final_row"].get("state", "?")

    print()
    print("  Initial state (state 0):")
    if initial_joint is not None:
        print(f"    Joint = {initial_joint:.4f}", end="")
        exp = meta.get("initial_joint")
        if exp is not None:
            dev    = format_deviation(initial_joint, exp)
            status = "OK" if abs(initial_joint - exp) < 0.01 else "MISMATCH"
            print(f"    expected: {exp:.4f}   deviation: {dev}  [{status}]", end="")
        print()
    else:
        print("    (no Joint value found)")

    print()
    print(f"  Final state (state {final_state}):")
    if final_joint is not None:
        print(f"    Joint = {final_joint:.4f}", end="")
        exp = meta.get("final_joint")
        if exp is not None:
            dev    = format_deviation(final_joint, exp)
            status = "OK" if abs(final_joint - exp) < 0.01 else "MISMATCH"
            print(f"    expected: {exp:.4f}   deviation: {dev}  [{status}]", end="")
        print()
    else:
        print("    (no Joint value found)")

    print()
    beast_mins = parsed.get("beast_minutes")
    if beast_mins is not None:
        print(f"  BEAST runtime:  {beast_mins:.2f} minutes ({beast_mins * 60:.0f} s)")
    else:
        print(f"  Wall-clock time: {elapsed:.1f} s ({elapsed / 60:.2f} minutes)")
    if parsed["final_rate"] is not None:
        print(f"  Rate:           {parsed['final_rate']:.3f} minutes/million states")

    op_start = stdout_text.find("Operator analysis")
    if op_start != -1:
        op_section = stdout_text[op_start:].split("\n\n")[0]
        print()
        print("  " + op_section.replace("\n", "\n  "))

    print("=" * 70)

    # Write markdown report
    report_path   = run_dir / f"{run_name}_report.md"
    run_ts_pretty = datetime.strptime(timestamp, "%Y%m%d_%H%M%S").strftime("%Y-%m-%d %H:%M:%S")
    write_report(
        report_path=report_path,
        run_name=run_name,
        btype=btype,
        bnum=bnum,
        model=model,
        seed=seed,
        timestamp=run_ts_pretty,
        meta=meta,
        sys_info=sys_info,
        versions=versions,
        parsed=parsed,
        elapsed=elapsed,
        stdout_text=stdout_text,
        cmd=cmd,
    )
    print(f"  Report:       {run_dir.name}/{report_path.name}")

    json_path = run_dir / f"{run_name}_results.json"
    json_data = write_json(
        json_path=json_path,
        btype=btype, bnum=bnum, model=model, seed=seed,
        timestamp=run_ts_pretty,
        meta=meta, sys_info=sys_info, versions=versions,
        parsed=parsed, elapsed=elapsed,
    )
    print(f"  Results JSON: {run_dir.name}/{json_path.name}")

    if args.append:
        append_csv(CSV_PATH, json_data)

    # Non-zero exit if either Joint density doesn't match expected (both are deterministic for a given seed)
    failed = False
    if initial_joint is not None and meta.get("initial_joint") is not None:
        if abs(initial_joint - meta["initial_joint"]) >= 0.01:
            print("\n  FAIL: Initial Joint density does not match expected value.")
            failed = True
    if final_joint is not None and meta.get("final_joint") is not None:
        if abs(final_joint - meta["final_joint"]) >= 0.01:
            print("\n  FAIL: Final Joint density does not match expected value.")
            failed = True
    if failed:
        sys.exit(2)


if __name__ == "__main__":
    main()
