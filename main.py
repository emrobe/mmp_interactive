import os
from os.path import dirname, join

import numpy as np
import pandas.io.sql as psql
import pandas
import sqlite3 as sql
from collections import Counter

from bokeh.plotting import figure
from bokeh.layouts import layout, column
from bokeh.models import ColumnDataSource, Div, OpenURL, TapTool, NumeralTickFormatter
from bokeh.models.widgets import Slider, Select, TextInput
from bokeh.io import curdoc
from bokeh.palettes import Oranges, Greens, Blues, Purples

# Paths to data
parsed_path = join(os.getcwd(), 'mmp_interactive/inputdata/input.tsv')

# Parse data
parsed = pandas.read_csv(parsed_path, sep='\t')
total = len(parsed)

axis_map = {
	#"tRNA Count": "tRNAs",
	"Total tRNAs":"Total_tRNAs",
	"Unique tRNAs":"Unique_tRNAs",
	"5S rRNA count": "rRNA5S",
	"16S rRNA count": "rRNA16S",
	"23S rRNA count": "rRNA23S",
        "Sample Depth": "depth",
        "Salinity": "env_salinity",
        "Temperature": "env_temp",
        "Replicons": "num_replicons",
        "Gene Count": "genes",
        "CDS Count": "cds",
        "CDS Count (pseudo)": "pseudo_genes",
        "CDS Count (frameshifted)": "frameshifted_genes",
        "Temperature (optimal)": "optimal_temperature",
        "Plasmid Count": "plasmids",
        "Sequencing Depth": "sequencing_depth",
	"GC Content": "GC",
	"Contigs Count": "contigs",
	"Genome Length": "Assembly_length",
	#"N50 Assembly statistics": "N50",
	#"CheckM genomes": "genomes",
	#"CheckM markers": "markers",
	#"CheckM Marker sets": "marker_sets",
	"CheckM Completeness": "Completeness",
	"CheckM Contamination": "Contamination",
	"CheckM Strain Heterogeneity": "Strain_heterogeneity",
	"Quality Score": "QS"

}


desc = Div(text=open(join(dirname(__file__), "description.html")).read(), sizing_mode="stretch_width")
qual_desc = Div(text=open(join(dirname(__file__), "quality_explaination.html")).read(), sizing_mode="stretch_width")

# Create Input controls
database = Select(title="Database", value="All", options=['All', 'MarRef', 'MarDB'])
analysis_type = Select(title="Type", value="All", options=['All', 'Whole genome sequencing (WGS)', 'Metagenome assembled genome (MAG)', 'Single amplified genome (SAG)'])
draft = Select(title="Draft Quality *", value="All", options=['All', 'Finished', 'High Quality Draft', 'Near Complete', 'Medium Quality Draft', 'Low Quality Draft', 'Very Low Quality Draft'])
minqsscore = Slider(title="Minimum QS score", value=min(parsed['QS']), start=min(parsed['QS']), end=max(parsed['QS']), step=5)
maxqsscore = Slider(title="Maximum QS score", value=max(parsed['QS']), start=-50, end=max(parsed['QS']), step=5)
minlength = Slider(title="Minimum Genome Length", value=1, start=1, end=15000000, step=100000)
maxlength = Slider(title="Maximum Genome Length", value=max(parsed['Assembly_length']), start=0, end=max(parsed['Assembly_length']), step=500000)
mincontigs = Slider(title="Minimum Number of Contigs", value=1, start=1, end=2000, step=5)
maxcontigs = Slider(title="Maximum Number of Contigs", value=max(parsed['contigs']), start=1, end=max(parsed['contigs']), step=50)
x_axis = Select(title="X Axis", options=sorted(axis_map.keys()), value="Genome Length")
y_axis = Select(title="Y Axis", options=sorted(axis_map.keys()), value="CheckM Completeness")

# Create Column Data Source that will be used by the plot
source = ColumnDataSource(data=dict(x=[], y=[], color=[], db=[], mmpid=[], apt=[], fsn=[], alpha=[], label=[], qual=[], comp=[], cont=[]))

TOOLTIPS=[
    ("Database", "@db"),
    ("Scientific name", "@fsn"),
    ("ID", "@mmpid"),
    ("Type", "@apt"),
    ("Quality", "@qual"),
    ("Completeness", "@comp"),
    ("Contamination", "@cont")
]

p = figure(plot_height=800, plot_width=1500, title="", toolbar_location=None, tooltips=TOOLTIPS, sizing_mode="fixed", tools="tap")
p.circle(x="x", y="y", source=source, size=7, color="color", line_color=None, fill_alpha="alpha", legend='label')


def select_entries():
    # Get widget values
    db_val = database.value
    type_val = analysis_type.value
    qual_val = draft.value

    # Filter dataset based on widget values
    selected = parsed[
        (parsed.QS >= minqsscore.value) &
        (parsed.QS <= maxqsscore.value) &
        (parsed.contigs >= mincontigs.value) &
        (parsed.contigs <= maxcontigs.value) &
        (parsed.Assembly_length >= minlength.value) &
        (parsed.Assembly_length <= maxlength.value)
    ]
    
    # Filter based on categories
    if (qual_val != "All"):
        selected = selected[selected.quality.str.contains(qual_val)==True]

    if (db_val != "All"):
        selected = selected[selected.DB.str.contains(db_val)==True]
    if (type_val != "All"):
        if (type_val == "Whole genome sequencing (WGS)"):
            selected = selected[selected.analysis_project_type.str.contains("WGS")==True]
        if (type_val == "Metagenome assembled genome (MAG)"):
            selected = selected[selected.analysis_project_type.str.contains("MAG")==True]
        if (type_val == "Single amplified genome (SAG)"):
            selected = selected[selected.analysis_project_type.str.contains("SAG")==True]
    return selected


def update():
    df = select_entries()
    x_name = axis_map[x_axis.value]
    y_name = axis_map[y_axis.value]

    # Calculate selection statistics
    missing=0
    for pos, values in enumerate(zip(df[x_name], df[y_name]),0):
       if np.isnan(values[0]) or np.isnan(values[1]):
               missing+=1

    # Axes text from widgets
    p.xaxis.axis_label = x_axis.value
    p.yaxis.axis_label = y_axis.value
    # 10, 1k, 124m formatting for numbers on x/y axes
    p.xaxis.formatter=NumeralTickFormatter(format="0a")
    p.yaxis.formatter=NumeralTickFormatter(format="0a")
    # Text fontsize on axes
    p.xaxis.axis_label_text_font_size = "14pt"
    p.yaxis.axis_label_text_font_size = "14pt"
    # Number (major ticks) fontsize on axes
    p.xaxis.major_label_text_font_size = "10pt"
    p.yaxis.major_label_text_font_size = "10pt"
    # Title
    p.title.text = "%d entries shown" % len(df)
    p.title.text = "Showing {} entries out of {}. ({} entries are filtered using widgets or have missing data for either X or Y)".format(total-((total-len(df))+missing), total, (total-len(df))+missing)
    # Taptool specifics, links etc.
    url = "https://mmp.sfb.uit.no/databases/@db/#/records/@mmpid"
    taptool = p.select(type=TapTool)
    taptool.callback = OpenURL(url=url)

    source.data = dict(
        x=df[x_name],
        y=df[y_name],
        color=df["color"],
        db=df["db"],
        mmpid=df["mmp_ID"],
        fsn=df["full_scientific_name"],
        apt=df["analysis_project_type"],
        alpha=df["alpha"],
        label=df['label'],
        qual=df['quality'],
        comp=df['Completeness'],
        cont=df['Contamination']
    )

controls = [database, analysis_type, draft, minqsscore, maxqsscore, minlength, maxlength, mincontigs, maxcontigs, x_axis, y_axis]
for control in controls:
    control.on_change('value', lambda attr, old, new: update())
controls.append(qual_desc)

inputs = column(*controls, width=320, height=1000)
inputs.sizing_mode = "fixed"
l = layout([
    [desc],
    [inputs, p],
], sizing_mode="scale_both", height=600, width=800)

# initial load of the data
update()  

curdoc().add_root(l)
curdoc().title = "MMP INTERACTIVE EXPLORER OF NUMERIC METADATA"
