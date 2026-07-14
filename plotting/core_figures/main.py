from plot_all_figures import main as plot_all_figures


def main():
    font_config = {
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

    plot_all_figures(
        font_config=font_config,
        font_delta=0,
    )


if __name__ == "__main__":
    main()