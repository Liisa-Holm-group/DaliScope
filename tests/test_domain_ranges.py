"""Regression checks for domain boundaries shared by UI and analysis APIs."""
import numpy as np
import pandas as pd
import pytest


def test_viewer_converts_displayed_domains_without_changing_ui_values():
    from daliscope.optics.domnet_viewer import DomNetViewer

    viewer = DomNetViewer(domain_string="1-129_389-404, 247-388")
    assert viewer.clipping_domain_string == "0-129_388-404, 246-388"
    assert viewer.domain_string == "1-129_389-404, 247-388"
    assert viewer.canvas.domain_string == viewer.domain_string
    viewer.domain_string = "23-123"
    assert viewer.clipping_domain_string == "22-123"


def test_domain_string_roundtrip_keeps_discontinuous_and_single_residue_domains():
    from daliscope.mechanics.domain_ranges import (
        clipping_domain_string_to_pdb,
        pdb_domain_string_to_clipping,
    )

    displayed = "1-1_100-150, 151-200"
    clipped = pdb_domain_string_to_clipping(displayed)
    assert clipped == "0-1_99-150, 150-200"
    assert clipping_domain_string_to_pdb(clipped) == displayed
    assert pdb_domain_string_to_clipping("") == ""


def test_clipping_and_community_select_the_same_gold_query_residues():
    from daliscope.analysis.structural_fingerprints import build_query_mask
    from daliscope.mechanics.domain_ranges import (
        clipping_ranges_to_pdb,
        pdb_ranges_to_clipping,
    )
    from daliscope.mechanics.metrics import clip_segments_to_domain_ranges

    clipping_ranges = pdb_ranges_to_clipping([(23, 50), (60, 123)])
    assert clipping_ranges == [(22, 50), (59, 123)]
    segments = pd.DataFrame({"alignment_id": [1], "q_start": [0], "s_start": [0], "length": [123]})
    clipped = clip_segments_to_domain_ranges(segments, clipping_ranges)
    clipped_mask = np.zeros(123, dtype=bool)
    for segment in clipped.itertuples():
        clipped_mask[segment.q_start:segment.q_start + segment.length] = True
    community_mask = build_query_mask(123, clipping_ranges_to_pdb(clipping_ranges))
    np.testing.assert_array_equal(clipped_mask, community_mask)
    assert np.flatnonzero(clipped_mask)[0] + 1 == 23
    assert np.flatnonzero(clipped_mask)[-1] + 1 == 123


@pytest.mark.parametrize("text", ["0-50", "5-4", "1-two", "1-4,,8-12", "1-2_", "1-2_0-4"])
def test_invalid_displayed_domain_ranges_fail_before_registration(text):
    from daliscope.mechanics.domain_ranges import pdb_domain_string_to_clipping

    with pytest.raises(ValueError):
        pdb_domain_string_to_clipping(text)


@pytest.mark.parametrize("ranges", [[(-1, 10)], [(3, 3)], [(1.5, 3)], [(True, 3)]])
def test_invalid_clipping_ranges_fail(ranges):
    from daliscope.mechanics.domain_ranges import clipping_ranges_to_pdb

    with pytest.raises(ValueError):
        clipping_ranges_to_pdb(ranges)
