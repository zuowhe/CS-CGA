import csv
import math
from pathlib import Path

from PIL import Image, ImageDraw, ImageFont
from reportlab.pdfgen import canvas


ROOT = Path(__file__).resolve().parent
OUTPUT_DIR = ROOT / "output"
DATA_DIR = ROOT / "data"
OUTPUT_DIR.mkdir(exist_ok=True)

CONFIGURATIONS = [
    "HSM-10", "HSM-20", "HSM-50", "HSM-100", "HSM-200",
    "HSM-1n", "HSM-2n", "HSM-5n", "HSM-10n", "HSM-50n", "HSM-100n", "HSM-200n",
]
DATASETS = ["Asia", "INS", "Water", "Alarm", "Hailfinder", "HEPAR", "Win95pts", "AND"]
DATASET_LABELS = ["Asia", "Insurance", "Water", "Alarm", "Hailfinder", "Hepar II", "Win95pts", "Andes"]


def load_rows():
    with (DATA_DIR / "hsm_capacity_latest_valid.csv").open(newline="", encoding="utf-8-sig") as source:
        return list(csv.DictReader(source))


def interpolate(start, end, position):
    position = max(0.0, min(1.0, position))
    red = round(start[0] + (end[0] - start[0]) * position)
    green = round(start[1] + (end[1] - start[1]) * position)
    blue = round(start[2] + (end[2] - start[2]) * position)
    return red, green, blue


def color_hex(red, green, blue):
    return f"#{red:02x}{green:02x}{blue:02x}"


def color_for(value, low, high, metric):
    position = (value - low) / (high - low) if high != low else 0.5
    if metric == "local_hit":
        return interpolate((166, 192, 254), (246, 128, 132), position)
    if metric == "entry_retention":
        return interpolate((172, 182, 229), (134, 253, 232), position)
    return interpolate((253, 219, 146), (209, 253, 255), position)


def text_color(rgb):
    red, green, blue = rgb
    return "#ffffff" if 0.299 * red + 0.587 * green + 0.114 * blue < 145 else "#111111"


def font(size, bold=False):
    name = "arialbd.ttf" if bold else "arial.ttf"
    return ImageFont.truetype(Path("C:/Windows/Fonts") / name, size)


def data_matrix(rows, method, metric):
    selected = [row for row in rows if row["method"] == method and row["sampleSize"] == "1000"]
    lookup = {(row["variant"], row["dataset"]): row for row in selected}
    matrix = []
    for configuration in CONFIGURATIONS:
        values = []
        for dataset in DATASETS:
            current = lookup[(configuration, dataset)]
            baseline = lookup[("HSM-large", dataset)]
            if metric == "runtime":
                value = float(current["runtimeSec"]) / float(baseline["runtimeSec"])
            elif metric == "local_hit":
                value = 100 * (float(current["localHitRate"]) - float(baseline["localHitRate"]))
            else:
                value = min(100.0, 100 * float(current["localEntries"]) / float(baseline["localEntries"]))
            values.append(value)
        matrix.append(values)
    return matrix


def value_label(value, metric):
    if metric == "runtime":
        return f"{value:.2f}x"
    if metric == "entry_retention":
        return f"{value:.1f}%"
    return f"{value:.2f}"


def chart_range(matrix, metric):
    values = [value for row in matrix for value in row]
    if metric == "runtime":
        return 0.8, max(1.5, math.ceil(max(values) * 10) / 10)
    if metric == "local_hit":
        return min(-20.0, math.floor(min(values) / 5) * 5), 0.0
    return 0.0, max(100.0, math.ceil(max(values) / 10) * 10)


def draw_heatmap(rows, method, metric, title, colorbar_label, filename):
    matrix = data_matrix(rows, method, metric)
    low, high = chart_range(matrix, metric)
    width, height = 1750, 1080
    left, top, cell_width, cell_height = 270, 105, 165, 68
    chart_width = cell_width * len(DATASETS)
    chart_height = cell_height * len(CONFIGURATIONS)
    chunks = [
        '<svg xmlns="http://www.w3.org/2000/svg" width="1750" height="1200" viewBox="0 0 1750 1200">',
        '<rect width="100%" height="100%" fill="white"/>',
        '<style>text { font-family: Arial, Helvetica, sans-serif; } .axis { font-size: 32px; } .network { font-size: 30px; } .label { font-size: 26px; } .value { font-size: 32px; font-weight: 500; } .tick { font-size: 32px; }</style>',
    ]

    for index, label in enumerate(DATASET_LABELS):
        x = left + index * cell_width + cell_width / 2
        chunks.append(f'<text x="{x}" y="84" text-anchor="middle" class="network">{label}</text>')

    for row_index, configuration in enumerate(CONFIGURATIONS):
        y = top + row_index * cell_height
        chunks.append(f'<text x="{left - 18}" y="{y + 44}" text-anchor="end" class="label">{configuration}</text>')
        for column_index, value in enumerate(matrix[row_index]):
            x = left + column_index * cell_width
            rgb = color_for(value, low, high, metric)
            chunks.append(f'<rect x="{x}" y="{y}" width="{cell_width}" height="{cell_height}" fill="{color_hex(*rgb)}" stroke="#ffffff" stroke-width="1"/>')
            chunks.append(f'<text x="{x + cell_width / 2}" y="{y + 44}" text-anchor="middle" class="value" fill="{text_color(rgb)}">{value_label(value, metric)}</text>')

    bar_x, bar_y, bar_width, bar_height = 1600, top + 30, 60, chart_height - 60
    steps = 120
    for index in range(steps):
        fraction = index / (steps - 1)
        value = high - fraction * (high - low)
        rgb = color_for(value, low, high, metric)
        y = bar_y + index * bar_height / steps
        chunks.append(f'<rect x="{bar_x}" y="{y}" width="{bar_width}" height="{bar_height / steps + 1}" fill="{color_hex(*rgb)}"/>')
    chunks.append(f'<rect x="{bar_x}" y="{bar_y}" width="{bar_width}" height="{bar_height}" fill="none" stroke="#333333"/>')
    chunks.append(f'<text x="{bar_x + bar_width / 2}" y="{bar_y - 12}" text-anchor="middle" class="tick">{high:.2f}</text>')
    chunks.append(f'<text x="{bar_x + bar_width / 2}" y="{bar_y + bar_height + 32}" text-anchor="middle" class="tick">{low:.2f}</text>')
    chunks.append(f'<text x="{bar_x + bar_width / 2 - 10}" y="{bar_y + bar_height / 2}" text-anchor="middle" transform="rotate(-90 {bar_x + bar_width / 2 - 10} {bar_y + bar_height / 2})" class="axis">{colorbar_label}</text>')
    chunks.append('</svg>')
    (OUTPUT_DIR / filename).write_text("\n".join(chunks), encoding="utf-8")


def draw_heatmap_png(rows, method, metric, title, colorbar_label, filename):
    matrix = data_matrix(rows, method, metric)
    low, high = chart_range(matrix, metric)
    width, height = 1750, 1080
    left, top, cell_width, cell_height = 270, 105, 165, 68
    chart_width = cell_width * len(DATASETS)
    chart_height = cell_height * len(CONFIGURATIONS)
    image = Image.new("RGB", (width, height), "white")
    draw = ImageDraw.Draw(image)
    network_font = font(30)
    label_font = font(26)
    value_font = font(32, bold=True)
    tick_font = font(32)

    for index, label in enumerate(DATASET_LABELS):
        x = left + index * cell_width + cell_width / 2
        draw.text((x, 84), label, fill="#111111", font=network_font, anchor="mm")

    for row_index, configuration in enumerate(CONFIGURATIONS):
        y = top + row_index * cell_height
        draw.text((left - 18, y + 34), configuration, fill="#111111", font=label_font, anchor="rm")
        for column_index, value in enumerate(matrix[row_index]):
            x = left + column_index * cell_width
            rgb = color_for(value, low, high, metric)
            draw.rectangle((x, y, x + cell_width, y + cell_height), fill=rgb, outline="white")
            draw.text(
                (x + cell_width / 2, y + 34),
                value_label(value, metric),
                fill=text_color(rgb),
                font=value_font,
                anchor="mm",
            )

    bar_x, bar_y, bar_width, bar_height = 1600, top + 30, 60, chart_height - 60
    steps = 120
    for index in range(steps):
        fraction = index / (steps - 1)
        value = high - fraction * (high - low)
        rgb = color_for(value, low, high, metric)
        y = bar_y + index * bar_height / steps
        draw.rectangle((bar_x, y, bar_x + bar_width, y + bar_height / steps + 1), fill=rgb)
    draw.rectangle((bar_x, bar_y, bar_x + bar_width, bar_y + bar_height), outline="#333333", width=1)
    draw.text((bar_x + bar_width / 2, bar_y - 12), f"{high:.2f}", fill="#111111", font=tick_font, anchor="ms")
    draw.text((bar_x + bar_width / 2, bar_y + bar_height + 32), f"{low:.2f}", fill="#111111", font=tick_font, anchor="ma")
    image.save(OUTPUT_DIR / filename)


def draw_heatmap_pdf(rows, method, metric, title, colorbar_label, filename):
    matrix = data_matrix(rows, method, metric)
    low, high = chart_range(matrix, metric)
    width, height = 1750, 1080
    crop_left, crop_bottom, crop_right, crop_top = 100, 135, 1750, 1025
    left, top, cell_width, cell_height = 270, 105, 165, 68
    chart_width = cell_width * len(DATASETS)
    chart_height = cell_height * len(CONFIGURATIONS)
    pdf = canvas.Canvas(
        str(OUTPUT_DIR / filename),
        pagesize=(crop_right - crop_left, crop_top - crop_bottom),
    )
    pdf.setTitle(title)
    pdf.translate(-crop_left, -crop_bottom)

    def set_color(rgb):
        pdf.setFillColorRGB(rgb[0] / 255, rgb[1] / 255, rgb[2] / 255)

    pdf.setFillColorRGB(1, 1, 1)
    pdf.rect(0, 0, width, height, stroke=0, fill=1)
    pdf.setFillColorRGB(17 / 255, 17 / 255, 17 / 255)
    pdf.setFont("Helvetica", 32)
    for index, label in enumerate(DATASET_LABELS):
        x = left + index * cell_width + cell_width / 2
        pdf.drawCentredString(x, height - 92, label)

    for row_index, configuration in enumerate(CONFIGURATIONS):
        y_top = top + row_index * cell_height
        pdf.setFillColorRGB(17 / 255, 17 / 255, 17 / 255)
        pdf.setFont("Helvetica", 28)
        pdf.drawRightString(left - 18, height - y_top - 43, configuration)
        for column_index, value in enumerate(matrix[row_index]):
            x = left + column_index * cell_width
            rgb = color_for(value, low, high, metric)
            set_color(rgb)
            pdf.setStrokeColorRGB(1, 1, 1)
            pdf.rect(x, height - y_top - cell_height, cell_width, cell_height, stroke=1, fill=1)
            pdf.setFont("Helvetica-Bold", 36)
            pdf.setFillColor(text_color(rgb))
            pdf.drawCentredString(x + cell_width / 2, height - y_top - 43, value_label(value, metric))

    bar_x, bar_y, bar_width, bar_height = 1600, top + 30, 60, chart_height - 60
    steps = 120
    for index in range(steps):
        fraction = index / (steps - 1)
        value = high - fraction * (high - low)
        rgb = color_for(value, low, high, metric)
        set_color(rgb)
        y_top = bar_y + index * bar_height / steps
        pdf.rect(bar_x, height - y_top - bar_height / steps - 1, bar_width, bar_height / steps + 1, stroke=0, fill=1)
    pdf.setStrokeColorRGB(51 / 255, 51 / 255, 51 / 255)
    pdf.rect(bar_x, height - bar_y - bar_height, bar_width, bar_height, stroke=1, fill=0)
    pdf.setFillColorRGB(17 / 255, 17 / 255, 17 / 255)
    pdf.setFont("Helvetica", 32)
    pdf.drawCentredString(bar_x + bar_width / 2, height - bar_y + 8, f"{high:.2f}")
    pdf.drawCentredString(bar_x + bar_width / 2, height - bar_y - bar_height - 32, f"{low:.2f}")
    pdf.setFont("Helvetica", 32)
    pdf.saveState()
    pdf.translate(bar_x + bar_width / 2 - 10, height - bar_y - bar_height / 2)
    pdf.rotate(90)
    pdf.drawCentredString(0, -15, colorbar_label)
    pdf.restoreState()
    pdf.save()


def main():
    rows = load_rows()
    figures = [
        ("CS-CGA", "runtime", "CS-CGA: Runtime Ratio to HSM-large", "Runtime ratio", "cs-cga-runtime-ratio.svg"),
        ("CS-CGA", "local_hit", "CS-CGA: Local Cache Hit-Rate Difference from HSM-large", "Difference (percentage points)", "cs-cga-local-hit-gap.svg"),
        ("CS-CGA", "entry_retention", "CS-CGA: Relative Local Parent-Set Entry Count", "Entry count relative to HSM-large (%)", "cs-cga-local-entry-retention.svg"),
        ("AESL-GA", "runtime", "AESL-GA: Runtime Ratio to HSM-large", "Runtime ratio", "aesl-ga-runtime-ratio.svg"),
        ("AESL-GA", "local_hit", "AESL-GA: Local Cache Hit-Rate Difference from HSM-large", "Difference (percentage points)", "aesl-ga-local-hit-gap.svg"),
        ("AESL-GA", "entry_retention", "AESL-GA: Relative Local Parent-Set Entry Count", "Entry count relative to HSM-large (%)", "aesl-ga-local-entry-retention.svg"),
    ]
    for method, metric, title, colorbar_label, filename in figures:
        draw_heatmap(rows, method, metric, title, colorbar_label, filename)
        draw_heatmap_png(rows, method, metric, title, colorbar_label, filename.replace(".svg", ".png"))
        draw_heatmap_pdf(rows, method, metric, title, colorbar_label, filename.replace(".svg", ".pdf"))


if __name__ == "__main__":
    main()
