"""Exercise the plane bisector's selection API with a real widget and Project."""
import importlib

import matplotlib.pyplot as plt
import pandas as pd
import pytest

from daliscope.core.project import Project

bisector_module = importlib.import_module("daliscope.widgets.plane_bisector")


@pytest.fixture
def selector(monkeypatch):
    previous_backend = plt.get_backend()
    previous_interactive = plt.isinteractive()
    plt.switch_backend("Agg")
    monkeypatch.setattr(bisector_module, "display", lambda *_: None)
    project = Project()
    frame = pd.DataFrame({
        "target_id": list("abcdef"),
        "z_score": [1., 2., 3., 4., 5., 6.],
        "query_coverage": [.1, .2, .3, .4, .5, .6],
        "alignment_length": [10, 20, 30, 40, 50, 60],
        "sequence_identity": [10.] * 6,
        "description": ["Test hit"] * 6,
    })
    project.add_view("full", frame)
    widget = bisector_module.plane_bisector(project)
    try:
        yield widget
    finally:
        plt.close(widget.fig)
        for control in widget.controls.values():
            control.close()
        plt.switch_backend(previous_backend)
        plt.interactive(previous_interactive)


def test_get_selected_tracks_empty_and_restored_selection(selector):
    selected, params = selector.get_selected()
    assert selected["target_id"].tolist() == list("def")
    assert params["view_name"] == "full"

    selector.controls["angle"].value = 0
    selector.controls["y_sect"].value = selector.controls["y_sect"].max
    empty, empty_params = selector.get_selected()
    assert empty.empty
    assert empty_params["angle"] == 0

    selector.controls["angle"].value = -90
    selector.controls["y_sect"].value = .35
    restored, restored_params = selector.get_selected()
    pd.testing.assert_frame_equal(restored, selected)
    assert restored_params["angle"] == -90
