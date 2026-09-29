"""Exercise the notebook's bisector harvest with a real widget and Project."""
import importlib
import json
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import pytest

import daliscope
from daliscope.core.project import Project

bisector_module = importlib.import_module("daliscope.widgets.plane_bisector")
ROOT = Path(__file__).resolve().parents[1]


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
        yield project, widget
    finally:
        plt.close(widget.fig)
        for control in widget.controls.values():
            control.close()
        plt.switch_backend(previous_backend)
        plt.interactive(previous_interactive)


def harvest(namespace):
    notebook = json.loads((ROOT / "notebooks/RunDaliScope.ipynb").read_text(encoding="utf-8"))
    source = next("".join(cell.get("source", [])) for cell in notebook["cells"]
                  if 'show_heading("Locking selection")' in "".join(cell.get("source", [])))
    namespace.update(daliscope=daliscope, show_heading=lambda *_: None)
    namespace.setdefault("display", lambda *_: None)
    exec(compile(source, "RunDaliScope.ipynb:Locking selection", "exec"), namespace)


@pytest.mark.parametrize("prior_nonempty", [False, True])
def test_empty_bisector_selection_does_not_register_a_view(selector, capsys, prior_nonempty):
    project, widget = selector
    displayed = []
    namespace = {"project": project, "bisector": widget, "display": displayed.append}
    if prior_nonempty:
        harvest(namespace)
        displayed.clear()
    widget.controls["angle"].value = 0
    widget.controls["y_sect"].value = widget.controls["y_sect"].max
    assert widget.get_selected()[0].empty
    before_views = dict(project.views)
    before_provenance = dict(project.provenance)

    harvest(namespace)

    assert project.views.keys() == before_views.keys()
    assert project.provenance == before_provenance
    assert namespace["view_name"] is None
    assert displayed == []
    assert "No hits selected" in capsys.readouterr().out

    widget.controls["angle"].value = -90
    widget.controls["y_sect"].value = .35
    harvest(namespace)
    expected_name = "full_filtered_1" if prior_nonempty else "full_filtered_0"
    assert namespace["view_name"] == expected_name
    assert project.views[expected_name]["target_id"].tolist() == list("def")


def test_nonempty_selection_registers_rows_and_provenance_on_rerun(selector):
    project, widget = selector
    expected, params = widget.get_selected()
    assert expected["target_id"].tolist() == list("def")
    namespace = {"project": project, "bisector": widget}

    harvest(namespace)
    harvest(namespace)

    for name in ("full_filtered_0", "full_filtered_1"):
        pd.testing.assert_frame_equal(project.views[name], expected.reset_index(drop=True))
        assert project.provenance[name].function == "plane_bisector"
        assert project.provenance[name].parameters == params
        assert project.provenance[name].parent_objects == ["full"]
