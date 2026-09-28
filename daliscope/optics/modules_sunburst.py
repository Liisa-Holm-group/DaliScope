import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Wedge
from matplotlib.colors import to_rgb

import daliscope.analysis.communities


def _lighten_color(color, amount=0.55):
    """
    Mix a color with white.

    amount=0   -> original color
    amount=1   -> white
    """
    rgb = np.asarray(to_rgb(color))
    return tuple(rgb + (1.0 - rgb) * amount)


def plot_module_sunbursts(
    modules_df,
    *,
    modules=None,
    ncols=5,
    figsize=(12, 12),
    start_angle=90,
    inner_radius=0.48,
    outer_radius=0.90,
    ring_gap=0.025,
    clan_colors=None,
    clan_cmap="tab20",
    clan_label_threshold=0.08,
    pfam_label_threshold=0.10,
    pfam_lighten=0.55,
    edgecolor="white",
    linewidth=0.45,
    center_label=True,
    show_clan_counts=False,
    show_pfam_counts=False,
    fontsize_clan=8, #6.5,
    fontsize_pfam=6.5, #4.5,
    module_fontsize=8,
):
    """
    Plot one two-level sunburst per module in a grid.

    Input
    -----
    modules_df : DataFrame
        One row per target/domain. Required columns:

            module
            target_id
            clan
            pfam

        `pfam` is a subclass of `clan`.

    modules : sequence, optional
        Modules to plot. If None, all modules are plotted in sorted order.

    Hierarchy
    ---------
        clan
          └── pfam

    Inner ring
    ----------
        Pfam clan, opaque color.

    Outer ring
    ----------
        Pfam family, lightened version of its parent clan color.

    Sector size
    -----------
        Number of targets in the module.

    Returns
    -------
    fig, axes, clan_colors
    """

    required = {"module", "target_id", "clan", "pfam"}
    missing = required - set(modules_df.columns)

    if missing:
        raise ValueError(
            f"modules_df is missing required columns: {sorted(missing)}"
        )

    df = modules_df.copy()

    # ------------------------------------------------------------
    # Normalize annotations
    # ------------------------------------------------------------

    df["clan"] = df["clan"].fillna("Unassigned")
    df["pfam"] = df["pfam"].fillna("Unassigned")

    # ------------------------------------------------------------
    # Determine modules
    # ------------------------------------------------------------

    if modules is None:
        modules = sorted(df["module"].dropna().unique())

    modules = list(modules)

    if len(modules) == 0:
        raise ValueError("No modules found.")

    # ------------------------------------------------------------
    # Global clan color map
    #
    # Important: the same clan gets the same color in every module.
    # ------------------------------------------------------------

    all_clans = (
        df.loc[df["module"].isin(modules), "clan"]
        .drop_duplicates()
        .tolist()
    )

    if clan_colors is None:
        cmap = plt.get_cmap(clan_cmap)

        # Indices of tab20 colors to use.
        # Skip index 7, which is gray.
        palette_indices = [
            i for i in range(cmap.N)
            if i not in {14,15}
        ]


        assigned_clans = [
            clan for clan in all_clans
            if clan != "Unassigned"
        ]

        clan_colors = {
            clan: cmap(palette_indices[i % len(palette_indices)])
            for i, clan in enumerate(assigned_clans)
        }
        # patch
        clan_colors["CL0422"] = (0.8392156862745098, 0.15294117647058825, 0.1568627450980392, 1.0) 

        # Fixed neutral color for missing/unassigned annotation.
        if "Unassigned" in all_clans:
            clan_colors["Unassigned"] = "0.75"

    else:
        missing_clans = set(all_clans) - set(clan_colors)

        if missing_clans:
            raise ValueError(
                "clan_colors is missing colors for: "
                f"{sorted(missing_clans)}"
            )


    # ------------------------------------------------------------
    # Grid
    # ------------------------------------------------------------

    nrows = int(np.ceil(len(modules) / ncols))

    fig, axes = plt.subplots(
        nrows,
        ncols,
        figsize=figsize,
        squeeze=False,
    )

    axes = axes.ravel()

    # ------------------------------------------------------------
    # Plot each module
    # ------------------------------------------------------------

    for ax, module in zip(axes, modules):

        module_df = df.loc[
            df["module"] == module,
            ["target_id", "clan", "pfam"],
        ].copy()

        if module_df.empty:
            ax.axis("off")
            continue

        total = len(module_df)

        # --------------------------------------------------------
        # Count Pfams within clans
        # --------------------------------------------------------

        family_counts = (
            module_df
            .groupby(["clan", "pfam"], sort=False)
            .size()
            .reset_index(name="count")
        )

        clan_counts = (
            family_counts
            .groupby("clan", sort=False)["count"]
            .sum()
            .reset_index()
        )

        clans = clan_counts["clan"].tolist()

        # --------------------------------------------------------
        # Calculate clan angular ranges
        # --------------------------------------------------------

        clan_angles = {}

        angle = start_angle

        for _, row in clan_counts.iterrows():

            clan = row["clan"]
            count = row["count"]

            width = 360.0 * count / total

            clan_angles[clan] = (
                angle,
                angle + width,
            )

            angle += width

        # --------------------------------------------------------
        # Inner ring: CLANS
        # --------------------------------------------------------

        for _, row in clan_counts.iterrows():

            clan = row["clan"]
            count = row["count"]

            theta1, theta2 = clan_angles[clan]

            wedge = Wedge(
                center=(0, 0),
                r=inner_radius,
                theta1=theta1,
                theta2=theta2,
                width=inner_radius,
                facecolor=clan_colors[clan],
                edgecolor=edgecolor,
                linewidth=linewidth,
            )

            ax.add_patch(wedge)

            # Clan label
            fraction = count / total

            if fraction >= clan_label_threshold:

                theta_mid = np.deg2rad(
                    (theta1 + theta2) / 2
                )

                r_label = inner_radius * 0.50

                x = r_label * np.cos(theta_mid)
                y = r_label * np.sin(theta_mid)

                label = str(clan)

                if show_clan_counts:
                    label += f"\n{count}"

                ax.text(
                    x,
                    y,
                    label,
                    ha="center",
                    va="center",
                    fontsize=fontsize_clan,
                    clip_on=True,
                )

        # --------------------------------------------------------
        # Outer ring: PFAM FAMILIES
        # --------------------------------------------------------

        family_inner_radius = inner_radius + ring_gap
        family_width = outer_radius - family_inner_radius

        for clan in clans:

            theta1_clan, theta2_clan = clan_angles[clan]

            clan_families = family_counts.loc[
                family_counts["clan"] == clan
            ]

            clan_total = clan_families["count"].sum()

            angle = theta1_clan

            for _, row in clan_families.iterrows():

                pfam = row["pfam"]
                count = row["count"]

                width = (
                    (theta2_clan - theta1_clan)
                    * count
                    / clan_total
                )

                theta1 = angle
                theta2 = angle + width

                # Same hue as the parent clan, but lighter.
                pfam_color = _lighten_color(
                    clan_colors[clan],
                    amount=pfam_lighten,
                )

                wedge = Wedge(
                    center=(0, 0),
                    r=outer_radius,
                    theta1=theta1,
                    theta2=theta2,
                    width=family_width,
                    facecolor=pfam_color,
                    edgecolor=edgecolor,
                    linewidth=linewidth,
                )

                ax.add_patch(wedge)

                # Pfam label
                fraction = count / total

                if fraction >= pfam_label_threshold:

                    theta_mid = np.deg2rad(
                        (theta1 + theta2) / 2
                    )

                    r_label = (
                        family_inner_radius
                        + family_width * 0.53
                    )

                    x = r_label * np.cos(theta_mid)
                    y = r_label * np.sin(theta_mid)

                    label = str(pfam)

                    if show_pfam_counts:
                        label += f"\n{count}"

                    # Rotate labels tangentially but keep them
                    # approximately upright.
                    rotation = np.rad2deg(theta_mid) - 90

                    if 90 < rotation % 360 < 270:
                        rotation += 180

                    ax.text(
                        x,
                        y,
                        label,
                        ha="center",
                        va="center",
                        fontsize=fontsize_pfam,
                        rotation=rotation,
                        rotation_mode="anchor",
                        clip_on=True,
                    )

                angle = theta2

        # --------------------------------------------------------
        # Center
        # --------------------------------------------------------

        if center_label:
            ax.text(
                0,
                0,
                f"{total:,}",
                ha="center",
                va="center",
                fontsize=fontsize_clan,
                fontweight="bold",
            )

        # --------------------------------------------------------
        # Module label
        # --------------------------------------------------------

        ax.text(
            0.02,
            0.98,
            f"Module {module}",
            transform=ax.transAxes,
            ha="left",
            va="top",
            fontsize=module_fontsize,
            fontweight="bold",
        )

        # --------------------------------------------------------
        # Axes formatting
        # --------------------------------------------------------

        ax.set_aspect("equal")
        ax.axis("off")

        limit = outer_radius * 1.06
        ax.set_xlim(-limit, limit)
        ax.set_ylim(-limit, limit)

    # ------------------------------------------------------------
    # Turn off unused panels
    # ------------------------------------------------------------

    for ax in axes[len(modules):]:
        ax.axis("off")

    fig.subplots_adjust(
        left=0.01,
        right=0.99,
        bottom=0.01,
        top=0.99,
        wspace=0.02,
        hspace=0.04,
    )

    return fig, axes, clan_colors

#----------------

def module_recursion(project, i):
    view_name = f"dom_FOLD_module_{i}"

    modules_df_temp, n_modules_temp = (
        daliscope.analysis.communities.run_community_detection(
            view_name=view_name,
            domain_range=(23, 123),
            project=project,
            k_neighbors=200,
            max_workers=1,
            verbose=True,
        )
    )

    fig, axes, clan_colors = (
        plot_module_sunbursts(
            modules_df_temp,
            modules=range(1, n_modules_temp + 1),
            figsize=(8, 2),
        )
    )
