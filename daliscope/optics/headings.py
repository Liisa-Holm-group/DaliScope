from IPython.display import display, HTML

def heading_plain(title):        # style_0: just bold text, no decoration
    display(HTML(f"<h2 style='font-size:24px; font-weight:700; margin-bottom:8px;'>{title}</h2>"))

def heading_banner(title):       # style_1: colored background block
    display(HTML(f"""
    <div style='background:#2c5f8a; color:white; padding:10px 16px; 
                border-radius:6px; margin-bottom:12px;'>
        <h2 style='margin:0; font-size:22px;'>{title}</h2>
    </div>
    """))    

def heading_underlined(title):   # style_2: text with ruled line beneath
        display(HTML(f"""
    <h2 style='font-size:26px; font-weight:700; margin-bottom:4px;'>{title}</h2>
    <hr style='border: 2px solid steelblue; margin-top:0; margin-bottom:12px;'>
    """))