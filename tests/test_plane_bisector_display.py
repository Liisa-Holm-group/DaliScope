"""Display the live bisector canvas without an ipympl PNG snapshot."""
import importlib
from types import SimpleNamespace

import ipywidgets as widgets
from ipympl.backend_nbagg import Canvas
import matplotlib.pyplot as plt
import pandas as pd
import pytest

bisector_module = importlib.import_module("daliscope.widgets.plane_bisector")


@pytest.mark.parametrize("update_during_snapshot", [False, True])
def test_bisector_display_keeps_live_canvas_without_png_snapshot(monkeypatch, update_during_snapshot):
    previous_backend = plt.get_backend()
    previous_interactive = plt.isinteractive()
    previous_figures = set(plt.get_fignums())
    plt.switch_backend("module://ipympl.backend_nbagg")
    plt.ioff()
    displayed = []
    snapshots = []
    angle_slider = None
    original_repr = Canvas._repr_mimebundle_
    containers = []
    selector = None

    def snapshot(canvas, **kwargs):
        snapshots.append(canvas)
        if update_during_snapshot:
            # Schedule a real slider callback after savefig has captured an
            # artist for drawing, matching the detached-collection traceback.
            collection = canvas.figure.axes[0].collections[0]
            original_draw = collection.draw
            updated = False

            def draw(renderer):
                nonlocal updated
                if not updated:
                    updated = True
                    angle_slider.value = 0
                return original_draw(renderer)

            monkeypatch.setattr(collection, "draw", draw)
        return original_repr(canvas, **kwargs)

    def display(item):
        nonlocal angle_slider
        displayed.append(item)
        if isinstance(item, widgets.FloatSlider) and item.description == "Angle":
            angle_slider = item
        if isinstance(item, widgets.Box):
            containers.append(item)
        # Exercise the actual rich-display method rather than a display spy.
        item._repr_mimebundle_()

    monkeypatch.setattr(Canvas, "_repr_mimebundle_", snapshot)
    monkeypatch.setattr(bisector_module, "display", display)
    frame = pd.DataFrame({
        "target_id": list("abcdef"),
        "z_score": [1., 2., 3., 4., 5., 6.],
        "query_coverage": [.1, .2, .3, .4, .5, .6],
        "pfam": ["p", "q", "p", "q", "p", "q"],
    })
    try:
        selector = bisector_module.plane_bisector(SimpleNamespace(views={"full": frame}))
        assert snapshots == []
        canvas = selector.fig.canvas
        assert any(isinstance(item, widgets.Box) and canvas in item.children for item in displayed)
        live_comm = canvas.comm
        assert live_comm is not None
        assert selector.get_selected()[0]["target_id"].tolist() == list("def")

        # The child is still the original, functioning ipympl canvas.
        canvas._handle_message(None, {"type": "initialized"}, [])
        canvas._handle_message(None, {"type": "draw"}, [])
        selector.controls["angle"].value = 90
        canvas._handle_message(None, {"type": "draw"}, [])
        assert selector.get_selected()[0]["target_id"].tolist() == list("abc")
        assert canvas.comm is live_comm
        assert snapshots == []
    finally:
        for number in set(plt.get_fignums()) - previous_figures:
            plt.close(number)
        if selector is not None:
            for control in selector.controls.values():
                control.close()
        for item in containers:
            item.close()
        plt.switch_backend(previous_backend)
        plt.interactive(previous_interactive)
