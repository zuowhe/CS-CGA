import csv
import re
from collections import OrderedDict
from pathlib import Path
from statistics import mean, stdev


ROOT = Path(__file__).resolve().parents[2]
PAPER_ROOT = ROOT / "my paper" / "Result"
PAPER_TEX = ROOT / "my paper" / "CSCGA-elsarticle-num.tex"
RESULTS_ROOT = ROOT / "results"
OUTPUT_DIR = RESULTS_ROOT / "comparison_tables"

MAIN_NETWORK_ORDER = [
    "Asia",
    "INS",
    "Water",
    "Alarm",
    "Hailfinder",
    "HEPAR",
    "Win95pts",
    "AND",
]
MAIN_SIZES = [500, 1000, 3000]

PAPER_METHOD_DIRS = OrderedDict(
    [
        ("Paper_CS-CGA", "CS-CGA"),
        ("Paper_MIGA", "MIGA"),
        ("Paper_hybrid-SLA", "SLA"),
        ("Paper_EKGA", "EKGA"),
        ("Paper_AESL-GA", "AESL"),
        ("Paper_PSX", "PSX"),
        ("Paper_MAGABN", "MAGABN"),
    ]
)

PAPER_TEX_LABELS = {
    "F1": "tab:comparison_algorithms_percent",
    "SHD": "tab:shd_comparison",
    "time": "tab:runtime",
}

PAPER_TEX_COLUMNS = [
    "Paper_CS-CGA",
    "Paper_MIGA",
    "Paper_hybrid-SLA",
    "Paper_EKGA",
    "Paper_AESL-GA",
    "Paper_PSX",
    "Paper_MAGABN",
]

LATEX_DATASET_NAMES = {
    "AS": "Asia",
    "IN": "INS",
    "WA": "Water",
    "AL": "Alarm",
    "HA": "Hailfinder",
    "HE": "HEPAR",
    "WI": "Win95pts",
    "AN": "AND",
}

NEW_METHODS = [
    "CS-CGA",
    "GES",
    "IGES_RCI",
    "MMHC",
    "GTT",
    "BayesianSearch",
    "GSP",
    "SaiyanH",
    "GOBNILP",
    "GOBNILP_CSCIP",
]

NEW_METHOD_COLUMNS = OrderedDict((f"New_{m}" if m == "CS-CGA" else m, m) for m in NEW_METHODS)
SKIP_NAME_PARTS = ["_CE_", "_native", "_core", "_adjusted"]


def metric_file_name(metric):
    if metric == "F1":
        return "F1_result_summary.csv"
    if metric == "SHD":
        return "SHD_result_summary.csv"
    if metric == "time":
        return "Time_result_summary.csv"
    raise ValueError(metric)


def parse_paper_summary(path):
    values = {}
    if not path.exists():
        return values
    with path.open("r", encoding="utf-8-sig", newline="") as f:
        reader = csv.reader(f)
        next(reader, None)
        for row in reader:
            if len(row) < 3:
                continue
            network = row[0].strip()
            try:
                size = int(float(row[1]))
            except ValueError:
                continue
            if network in MAIN_NETWORK_ORDER and size in MAIN_SIZES:
                values[(network, size)] = row[2].strip()
    return values


def metric_value_index(metric):
    if metric == "F1":
        return 1
    if metric == "SHD":
        return 5
    if metric == "time":
        return 7
    raise ValueError(metric)


def format_metric_value(avg, sd, metric):
    if avg is None:
        return "--"
    if sd is None:
        sd = 0.0
    if metric == "F1":
        return f"{avg:.4f}({sd:.4f})"
    if metric == "SHD":
        return f"{avg:.2f}({sd:.2f})"
    return f"{avg:.1f}({sd:.1f})"


def strip_latex_cell(cell):
    text = cell.strip()
    text = text.replace(r"\textbf{", "")
    text = text.replace(r"\%", "%")
    text = text.replace("%", "")
    text = text.replace("{", "")
    text = text.replace("}", "")
    text = text.replace("$", "")
    return text.strip()


def format_latex_metric_cell(cell, metric):
    text = strip_latex_cell(cell)
    match = re.search(r"([0-9]+(?:\.[0-9]+)?)\s*\(\s*([0-9]+(?:\.[0-9]+)?)\s*\)", text)
    if not match:
        return text if text else "--"
    avg = float(match.group(1))
    sd = float(match.group(2))
    if metric == "F1":
        return format_metric_value(avg / 100.0, sd / 100.0, metric)
    return format_metric_value(avg, sd, metric)


def parse_paper_tex_table(metric):
    values = {column: {} for column in PAPER_TEX_COLUMNS}
    if not PAPER_TEX.exists():
        return values

    label = PAPER_TEX_LABELS[metric]
    tex = PAPER_TEX.read_text(encoding="utf-8", errors="ignore")
    label_pos = tex.find(rf"\label{{{label}}}")
    if label_pos < 0:
        return values
    tabular_start = tex.find(r"\begin{tabular}", label_pos)
    tabular_end = tex.find(r"\end{tabular}", tabular_start)
    if tabular_start < 0 or tabular_end < 0:
        return values
    table_text = tex[tabular_start:tabular_end]

    row_pattern = re.compile(r"^\s*(AS|IN|WA|AL|HA|HE|WI|AN)-(\d+)\s*&(.+?)\\\\", re.MULTILINE)
    for match in row_pattern.finditer(table_text):
        network = LATEX_DATASET_NAMES[match.group(1)]
        size = int(match.group(2))
        cells = [part.strip() for part in match.group(3).split("&")]
        if len(cells) < len(PAPER_TEX_COLUMNS):
            continue
        for column, cell in zip(PAPER_TEX_COLUMNS, cells):
            values[column][(network, size)] = format_latex_metric_cell(cell, metric)
    return values


def parse_number(token):
    token = token.strip()
    if not token:
        return None
    try:
        return float(token)
    except ValueError:
        return None


def split_result_line(line):
    return [x.strip() for x in line.strip().split(",")]


def summarize_raw_result(path, metric):
    trial_records = OrderedDict()
    if not path.exists():
        return "--"

    avg_row = None
    std_row = None
    value_index = metric_value_index(metric)

    with path.open("r", encoding="utf-8-sig", errors="ignore") as f:
        for line in f:
            if not line.strip():
                continue
            parts = split_result_line(line)
            if not parts:
                continue
            row_type = parts[0].lower()
            if row_type == "avg":
                avg_row = [parse_number(p) for p in parts]
                continue
            if row_type == "std":
                std_row = [parse_number(p) for p in parts]
                continue
            try:
                trial = int(float(parts[0]))
            except ValueError:
                continue

            if len(parts) > 1 and parts[1].strip().lower() == "timeout":
                trial_records[trial] = {"status": "timeout", "values": []}
                continue

            numeric = [parse_number(p) for p in parts]
            if len(numeric) < 8 or numeric[1] is None:
                trial_records[trial] = {"status": "invalid", "values": numeric}
            else:
                trial_records[trial] = {"status": "valid", "values": numeric}

    if avg_row is not None and len(avg_row) > value_index and avg_row[value_index] is not None:
        sd = None
        if std_row is not None and len(std_row) > value_index:
            sd = std_row[value_index]
        return format_metric_value(avg_row[value_index], sd, metric)

    if not trial_records:
        return "--"

    total_trials = len(trial_records)
    valid_values = []
    for rec in trial_records.values():
        if rec["status"] != "valid":
            continue
        nums = rec["values"]
        value = nums[value_index] if len(nums) > value_index else None
        if value is not None:
            valid_values.append(value)

    if not valid_values:
        return "Timeout" if any(r["status"] == "timeout" for r in trial_records.values()) else "--"

    avg = mean(valid_values)
    sd = stdev(valid_values) if len(valid_values) > 1 else 0.0
    text = format_metric_value(avg, sd, metric)

    if len(valid_values) != total_trials:
        text = f"{text} [n={len(valid_values)}/{total_trials}]"
    return text


def parse_paper_raw_filename(path):
    stem = path.stem
    for network in sorted(MAIN_NETWORK_ORDER, key=len, reverse=True):
        if not stem.startswith(network):
            continue
        rest = stem[len(network) :]
        match = re.match(r"^(\d+)_100_200_", rest)
        if not match:
            continue
        size = int(match.group(1))
        if size in MAIN_SIZES:
            return network, size
    return None


def parse_result_filename(path):
    name = path.name
    if any(part in name for part in SKIP_NAME_PARTS):
        return None
    stem = path.stem
    for method in sorted(NEW_METHODS, key=len, reverse=True):
        suffix = f"_100_200_{method}"
        if not stem.endswith(suffix):
            continue
        prefix = stem[: -len(suffix)]
        for network in sorted(MAIN_NETWORK_ORDER, key=len, reverse=True):
            if not prefix.startswith(network):
                continue
            size_text = prefix[len(network) :]
            if not size_text.isdigit():
                continue
            size = int(size_text)
            if size not in MAIN_SIZES:
                continue
            return network, size, method
    return None


def collect_latest_new_result_files():
    latest = {}
    for date_dir in RESULTS_ROOT.iterdir():
        if not date_dir.is_dir() or not re.fullmatch(r"\d{8}", date_dir.name):
            continue
        for path in date_dir.glob("*_100_200_*.csv"):
            parsed = parse_result_filename(path)
            if parsed is None:
                continue
            network, size, method = parsed
            key = (network, size, method)
            old = latest.get(key)
            if old is None or path.stat().st_mtime > old.stat().st_mtime:
                latest[key] = path
    return latest


def collect_paper_raw_results(method_dir, metric):
    values = {}
    source_dir = PAPER_ROOT / method_dir
    if not source_dir.exists():
        return values
    for path in source_dir.glob("*_100_200_*.csv"):
        parsed = parse_paper_raw_filename(path)
        if parsed is None:
            continue
        values[parsed] = summarize_raw_result(path, metric)
    return values


def build_metric_table(metric, latest_new_files):
    columns = ["Network", "DataSize"] + list(PAPER_METHOD_DIRS.keys()) + list(NEW_METHOD_COLUMNS.keys())
    paper_values_by_method = parse_paper_tex_table(metric)
    for method_column, method_dir in PAPER_METHOD_DIRS.items():
        if not paper_values_by_method.get(method_column):
            paper_values_by_method[method_column] = collect_paper_raw_results(method_dir, metric)

    rows = []
    for network in MAIN_NETWORK_ORDER:
        for size in MAIN_SIZES:
            row = OrderedDict((col, "--") for col in columns)
            row["Network"] = network
            row["DataSize"] = size

            for method_column, values in paper_values_by_method.items():
                row[method_column] = values.get((network, size), "--")

            for column, method in NEW_METHOD_COLUMNS.items():
                path = latest_new_files.get((network, size, method))
                if path is not None:
                    row[column] = summarize_raw_result(path, metric)

            rows.append(row)
    return columns, rows


def write_csv(path, columns, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8-sig", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=columns)
        writer.writeheader()
        writer.writerows(rows)


def main():
    latest_new_files = collect_latest_new_result_files()
    outputs = {}
    for metric, file_name in [
        ("F1", "comparison_F1.csv"),
        ("SHD", "comparison_SHD.csv"),
        ("time", "comparison_time.csv"),
    ]:
        columns, rows = build_metric_table(metric, latest_new_files)
        output_path = OUTPUT_DIR / file_name
        write_csv(output_path, columns, rows)
        outputs[metric] = output_path

    print("Generated comparison tables:")
    for metric, output_path in outputs.items():
        print(f"{metric}: {output_path}")
    print(f"New result files used: {len(latest_new_files)}")


if __name__ == "__main__":
    main()
