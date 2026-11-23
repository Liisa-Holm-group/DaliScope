from importlib import import_module

import ipywidgets as widgets
import pandas as pd
import plotly.graph_objects as go

from daliscope.core.project import Project

explorer = import_module("daliscope.widgets.plotly_explorer")


def test_lasso_and_save_keep_the_scatter_controls_visible(monkeypatch):
    # A callback output that Colab doesn't capture in an Output widget reaches
    # the cell. Clearing that cell removes its entire interactive panel.
    cell = []
    monkeypatch.setattr(explorer, "display", cell.append)

    project = Project()
    project.add_view("full", pd.DataFrame({
        "target_id": ["one", "two", "three", "four"],
        "z_score": [1., 2., 3., 4.], "query_coverage": [.1, .2, .3, .4],
        "clan": ["a", "a", "b", "b"], "pfam": ["f"] * 4,
    }))
    explorer.launch_interactive_selectable_scatter(project, "full")

    def descendants(widget):
        yield widget
        for child in getattr(widget, "children", ()):
            yield from descendants(child)

    def visible_widgets():
        return [item for root in cell for item in descendants(root)]

    def button(label):
        return next(w for w in visible_widgets()
                    if isinstance(w, widgets.Button) and w.description == label)

    assert any(isinstance(w, widgets.Button) and w.description == "Save Subset"
               for w in visible_widgets()), "initial rendering cleared the control panel"
    save = button("Save Subset")
    assert save.disabled
    figure = next(w for w in visible_widgets() if isinstance(w, go.FigureWidget))
    monkeypatch.setattr(widgets.widgets.widget_output, "clear_output",
                        lambda *args, **kwargs: cell.clear())
    figure._js2py_pointsCallback = {
        "event_type": "plotly_selected",
        "selector": {"type": "lasso", "selector_state": {"xs": [0, 4, 4], "ys": [0, 0, 1]}},
        "points": {"xs": [1., 4.], "ys": [.1, .4], "point_indexes": [0, 1], "trace_indexes": [0, 1]},
    }
    assert figure in visible_widgets()
    assert button("Save Subset") is save and not save.disabled
    assert any("Selected 2" in str(getattr(w, "value", "")) for w in visible_widgets())
    name = next(w for w in visible_widgets() if isinstance(w, widgets.Text))
    assert name.placeholder == "BOX_full_2_HITS"
    name.value = "full"
    save.click()
    assert not save.disabled and len(project.views["full"]) == 4
    name.value = "my_lasso"
    save.click()
    assert project.views["my_lasso"]["target_id"].tolist() == ["one", "four"]
    assert project.provenance["my_lasso"].parent_objects == ["full"]
    assert button("Save Subset") is save and save.disabled
    assert any(isinstance(w, go.FigureWidget) for w in visible_widgets())
    assert any("my_lasso" in str(getattr(w, "value", "")) for w in visible_widgets())
