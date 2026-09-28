import ipywidgets as widgets
from IPython.display import display, clear_output

def msa_widget(viz, query_length: int, mode='logo'):
    """
    Interactive logo/heatmap viewer with left/right range sliders.
    
    Parameters
    ----------
    viz          : msa.logo or msa.heatmap function with (left, right) arguments
    query_length : total length of query sequence — sets slider bounds
    """

    style  = {"description_width": "50px"}
    layout = widgets.Layout(width="400px")

    left_slider = widgets.IntSlider(
        value=0,
        min=0,
        max=query_length - 1,
        step=1,
        description="Left:",
        style=style,
        layout=layout,
      )
    right_slider = widgets.IntSlider(
        value=query_length,
        min=1,
        max=query_length,
        step=1,
        description="Right:",
        style=style,
        layout=layout,
    )

    plot_out = widgets.Output()

    def update(_):
        left  = left_slider.value
        right = right_slider.value

        # enforce left < right
        if left >= right:
            with plot_out:
                clear_output(wait=True)
                print(f"⚠️  Left ({left}) must be less than Right ({right})")
            return

        with plot_out:
            clear_output(wait=True)
            viz(left=left, right=right) # display function

    # link sliders — update on any change
    left_slider.observe(update,  names="value")
    right_slider.observe(update, names="value")

    # reset button
    reset_btn = widgets.Button(
        description="Reset",
        button_style="",
        icon="refresh",
        layout=widgets.Layout(width="90px"),
    )
    def on_reset(_):
        left_slider.value  = 0
        right_slider.value = query_length
    reset_btn.on_click(on_reset)

    ui = widgets.VBox([
        widgets.HTML(f"<b>Sequence {mode} range</b>"),
        left_slider,
        right_slider,
        reset_btn,
        plot_out,
    ])

    display(ui)
    update(None)   # draw immediately with default values
