"""Generate every paper result figure with centrally configured font sizes."""

import argparse
from copy import deepcopy

import plot_adaptive_threshold_miga
import plot_adaptive_threshold_sla
import plot_bic_convergence
import plot_normalized_radar_charts
import plot_time_ablation_runtime


FONT_CONFIG = {
    "bic_convergence": {
        "base": 22,
        "title": 28,
        "axis_label": 26,
        "tick": 22,
        "legend": 20,
        "offset": 22,
    },
    "adaptive_threshold": {
        "base": 20,
        "title": 26,
        "axis_label": 26,
        "tick": 23,
        "legend": 19,
    },
    "radar_charts": {
        "base": 22,
        "tick": 22,
        "legend": 20,
    },
    "time_ablation": {
        "base": 22,
        "axes_label": 29,
        "x_tick": 30,
        "y_tick": 32,
        "legend": 29,
        "y_label": 34,
    },
}

# Change this single value to increase or decrease every configured font size.
FONT_DELTA = 0


def build_font_config(font_config=None, font_delta=0):
    config = deepcopy(FONT_CONFIG)
    for figure_group, values in (font_config or {}).items():
        config.setdefault(figure_group, {}).update(values)

    if font_delta:
        for values in config.values():
            for key, value in values.items():
                values[key] = value + font_delta
    return config


def main(font_config=None, font_delta=FONT_DELTA):
    config = build_font_config(font_config, font_delta)
    plot_bic_convergence.main(config["bic_convergence"])
    plot_adaptive_threshold_sla.main(config["adaptive_threshold"])
    plot_adaptive_threshold_miga.main(config["adaptive_threshold"])
    plot_normalized_radar_charts.main(config["radar_charts"])
    plot_time_ablation_runtime.main(config["time_ablation"])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--font-delta",
        type=int,
        default=FONT_DELTA,
        help="Add this value to every configured font size.",
    )
    args = parser.parse_args()
    main(font_delta=args.font_delta)
