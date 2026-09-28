"""
Portable standalone DomNet viewer.

This class has no Jupyter or ipywidgets dependency.
"""

from pathlib import Path

from .domnet_viewer_core import (
    DEFAULT_COLORS,
    build_standalone_html,
)


class DomNetViewerHTML:
    """
    Generate a standalone interactive DomNet HTML viewer.

    Example
    -------

        viewer = DomNetViewerHTML(
            pdb_text=pdb_text,
            domain_string="1-402_725-761,403-724_762-823",
        )

        viewer.save("viewer.html")
    """

    def __init__(
        self,
        pdb_text,
        domain_string,
        width=900,
        height=650,
        colors=None,
    ):

        self.pdb_text = str(
            pdb_text
        )

        self.domain_string = str(
            domain_string
        )

        self.width = int(width)

        self.height = int(height)

        self.colors = (
            colors
            if colors is not None
            else DEFAULT_COLORS
        )


    def html(self):
        """
        Return the complete standalone HTML document.
        """

        return build_standalone_html(
            pdb_text=self.pdb_text,
            domain_string=self.domain_string,
            width=self.width,
            height=self.height,
            colors=self.colors,
        )


    def save(self, filename="viewer.html"):
        """
        Write the standalone viewer to disk.
        """

        path = Path(filename)

        path.write_text(
            self.html(),
            encoding="utf-8",
        )

        return path
