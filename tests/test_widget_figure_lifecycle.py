"""Keep earlier ipympl canvases alive while another widget redraws its output."""
import json
from pathlib import Path
from textwrap import dedent
from types import SimpleNamespace

import matplotlib.pyplot as plt
from matplotlib.figure import Figure
import pytest

from daliscope.widgets import factory

ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture
def plot_runtime(monkeypatch):
    previous_backend = plt.get_backend()
    previous_interactive = plt.isinteractive()
    previous_numbers = set(plt.get_fignums())
    plt.switch_backend("module://ipympl.backend_nbagg")
    plt.ioff()
    displayed = []
    widget_roots = []

    def display(item):
        displayed.append((item, item.canvas.comm is not None if isinstance(item, Figure) else None))

    monkeypatch.setattr(factory, "display", display)
    try:
        yield SimpleNamespace(displayed=displayed, widget_roots=widget_roots)
    finally:
        for number in set(plt.get_fignums()) - previous_numbers:
            plt.close(number)
        for root in widget_roots:
            for child in root.children:
                child.close()
            root.close()
        plt.switch_backend(previous_backend)
        plt.interactive(previous_interactive)


def test_widget_redraw_keeps_another_canvas_alive_and_closes_its_own_figures(plot_runtime):
    earlier = plt.figure()
    live_comm = earlier.canvas.comm
    rendered = []

    def plot(value):
        fig, ax = plt.subplots()
        ax.plot([0, 1], [0, value])
        rendered.append(fig)
        return fig

    root, get_state, controls = factory.widget_view(
        plot, params={"value": {"type": "intslider", "min": 0, "max": 2, "value": 0}},
        display_result=False, return_state=True, show_default_buttons=False,
    )
    plot_runtime.widget_roots.append(root)
    assert earlier.canvas.comm is live_comm
    assert rendered[0].canvas.comm is None

    controls["value"].value = 1
    assert get_state()["value"] == 1
    assert len(rendered) == 2
    assert earlier.canvas.comm is live_comm
    assert all(fig.canvas.comm is None for fig in rendered)
    assert [(fig, alive) for fig, alive in plot_runtime.displayed] == [(fig, True) for fig in rendered]


@pytest.mark.parametrize("interactive", [False, True])
def test_launcher_retains_prior_canvas_and_interactive_mode(plot_runtime, interactive):
    earlier = plt.figure()
    live_comm = earlier.canvas.comm
    rendered = []
    plt.interactive(interactive)

    def plot():
        fig = plt.figure()
        rendered.append(fig)
        return fig

    factory.launch_interactive_view(plot, title="Another plot", params={})

    assert earlier.canvas.comm is live_comm
    assert rendered[0].canvas.comm is None
    assert plt.isinteractive() is interactive
    plot_runtime.widget_roots.append(plot_runtime.displayed[-1][0])


def test_failed_render_cleans_only_new_figures_and_restores_interactive_mode(plot_runtime):
    earlier = plt.figure()
    live_comm = earlier.canvas.comm
    rendered = []
    plt.ion()

    def plot():
        rendered.append(plt.figure())
        raise RuntimeError("Rendering failed")

    with pytest.raises(RuntimeError, match="Rendering failed"):
        factory.widget_view(plot, params={}, display_result=False, show_default_buttons=False)

    assert earlier.canvas.comm is live_comm
    assert rendered[0].canvas.comm is None
    assert plt.isinteractive()


def test_profile_cell_preserves_prior_canvas_and_cleans_its_new_figures(plot_runtime):
    notebook = json.loads((ROOT / "notebooks/RunDaliScope-1.ipynb").read_text(encoding="utf-8"))
    source = next("".join(cell["source"]) for cell in notebook["cells"]
                  if "derive_signature_profile(" in "".join(cell.get("source", [])))
    start = source.index("    nr_threshold =")
    stop = source.index('    show_heading("Validation in 3D:')
    earlier = plt.figure()
    live_comm = earlier.canvas.comm
    rendered = []

    def derive(*args, **kwargs):
        rendered.append(plt.figure())
        return (None,) * 6

    namespace = {
        "plt": plt, "project": object(), "SEED_VIEW": "full",
        "daliscope": SimpleNamespace(analysis=SimpleNamespace(signature=SimpleNamespace(derive_signature_profile=derive))),
    }
    exec(dedent(source[start:stop]), namespace)

    assert earlier.canvas.comm is live_comm
    assert rendered[0].canvas.comm is None
