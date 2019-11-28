import os
from os.path import dirname, join

import numpy as np
import pandas.io.sql as psql
import pandas
import sqlite3 as sql
from collections import Counter

from bokeh.plotting import figure
from bokeh.layouts import layout, column
from bokeh.models import ColumnDataSource, Div, OpenURL, TapTool
from bokeh.models.widgets import Slider, Select, TextInput
from bokeh.io import curdoc
from bokeh.palettes import Oranges, Greens, Blues, Purples

# Paths to data
mardb_path = join(os.getcwd(), 'mmp_interactive/MarDB_merged_complete.tsv')
marref_path = join(os.getcwd(), 'mmp_interactive/MarRef_merged_done.tsv')

# Parse data
parsed_marref = pandas.read_csv(marref_path, sep='\t')
parsed_mardb = pandas.read_csv(mardb_path, sep='\t')

# Concatonate
parsed = pandas.concat([parsed_marref, parsed_mardb], ignore_index=True)

# Calculate Draft quality based on data from dataframe
quality = []
for index, row in parsed.iterrows():
    if row['DB'] == 'MarRef':
        quality.append('Finished')
    #elif row['rRNA5S'] > 0 and row['rRNA16S'] > 0 and row['rRNA23S'] > 0 and row['tRNAs'] > 18 and row['Completeness'] > 90 and row['Contamination'] < 5:
    elif row['Completeness'] > 90 and row['Contamination'] < 5:
        quality.append('High Quality Draft')
    elif row['Completeness'] >= 50 and row['Contamination'] < 10:
        quality.append('Medium Quality Draft')
    elif row['Completeness'] < 50 and row['Contamination'] < 10:
        quality.append('Low Quality Draft')
    elif row['Contamination'] > 10:
        quality.append('Very Low Quality Draft')
    else:
        quality.append('NA')
        #print (row)
#print (Counter(quality))
#exit()

# Set indivival colors for discrete databases
colorarray = []
for DB, apt in zip(parsed['DB'], parsed['analysis_project_type']):
	if DB == "MarRef" and apt == "missing":
		colorarray.append("grey")
	elif DB == "MarDB" and apt == "missing":
		colorarray.append("grey")
	elif DB == "MarRef":
		colorarray.append("blue")
	elif DB == "MarDB":
		colorarray.append("green")
	else:
		colorarray.append("grey")

# Tune alphas for analysis_project_type
alpha = []
for analysis_type in parsed['analysis_project_type']:
        if analysis_type == "Whole genome sequencing (WGS)":
                alpha.append(1)
        elif analysis_type == "Metagenome assembled genome (MAG)":
                alpha.append(0.6)
        elif analysis_type == "Single amplified genome (SAG)":
                alpha.append(0.4)
        else:
                alpha.append(1)

# Make labels for legend
label= []
for DB, apt in zip(parsed['DB'], parsed['analysis_project_type']):
	if DB == "MarRef" and apt == "missing":
		label.append("MarRef (Unknown type)")
	elif DB == "MarDB" and apt == "missing":
		label.append("MarDB (Unknown type)")
	elif DB == "MarRef" and apt == "Whole genome sequencing (WGS)":
		label.append("MarRef (Whole genome sequencing (WGS))")
	elif DB == "MarDB" and apt == "Whole genome sequencing (WGS)":
		label.append("MarDB (Whole genome sequencing (WGS))")
	elif DB == "MarRef" and apt == "Metagenome assembled genome (MAG)":
		label.append("MarRef (Metagenome assembled genome (MAG))")
	elif DB == "MarDB" and apt == "Metagenome assembled genome (MAG)":
		label.append("MarDB (Metagenome assembled genome (MAG))")
	elif DB == "MarRef" and apt == "Single amplified genome (SAG)":
		label.append("MarRef (Single amplified genome (SAG))")
	elif DB == "MarDB" and apt == "Single amplified genome (SAG)":
		label.append("MarDB (Single amplified genome (SAG))")
	else:
		label.append("Unknown")

# Include Color, alpha, labels (legend) and quality score (Completeness - 5xContamination)
parsed['color'] = colorarray
parsed['alpha'] = alpha
parsed['label'] = label
parsed['QS'] = parsed['Completeness'] - 5*parsed['Contamination']
parsed['quality'] = quality

# Quickfix: Add a column with lowercase db str we can use to neatly manage tooltip links to fit with mmp-address. Fix this in the future, it is unneccessary. 
parsed['db'] = parsed['DB'].str.lower()

# Hack to remove missing data (Remove #N/A)
missing=0
for pos, values in enumerate(zip(parsed['rRNA16S'], parsed['Completeness']),0):
	if np.isnan(values[0]) or np.isnan(values[1]):
		missing+=1
		#colorarray[pos] = 'red'
		#label[pos] = 'Missing data'

axis_map = {
	"tRNA count": "tRNAs",
	"5SrRNA count": "rRNA5S",
	"16SrRNA count": "rRNA16S",
	"23S rRNA count": "rRNA23S",
	"GC content": "GC",
	"Number of contigs": "contigs",
	"Assembly Length": "Assembly_length",
	"N50 Assembly statistics": "N50",
	"CheckM genomes": "genomes",
	"CheckM markers": "markers",
	"CheckM Marker sets": "marker_sets",
	"Checkm Completeness": "Completeness",
	"CheckM Contamination": "Contamination",
	"CheckM Strain Heterogeneity": "Strain_heterogeneity",
	"Quality Score": "QS"

}


desc = Div(text=open(join(dirname(__file__), "description.html")).read(), sizing_mode="stretch_width")
qual_desc = Div(text=open(join(dirname(__file__), "quality_explaination.html")).read(), sizing_mode="stretch_width")

# Create Input controls
database = Select(title="Database", value="All", options=['All', 'MarRef', 'MarDB'])
analysis_type = Select(title="Type", value="All", options=['All', 'Whole genome sequencing (WGS)', 'Metagenome assembled genome (MAG)', 'Single amplified genome (SAG)'])
draft = Select(title="Draft Quality *", value="All", options=['All', 'Finished', 'High Quality Draft', 'Medium Quality Draft', 'Low Quality Draft', 'Very Low Quality Draft'])
minqsscore = Slider(title="Minimum QS score", value=min(parsed['QS']), start=min(parsed['QS']), end=max(parsed['QS']), step=5)
maxqsscore = Slider(title="Maximum QS score", value=max(parsed['QS']), start=-50, end=max(parsed['QS']), step=5)
minlength = Slider(title="Minimum Assembly Length", value=1, start=1, end=15000000, step=100000)
maxlength = Slider(title="Maximum Assembly Length", value=max(parsed['Assembly_length']), start=0, end=max(parsed['Assembly_length']), step=500000)
mincontigs = Slider(title="Minimum Number of Contigs", value=1, start=1, end=2000, step=5)
maxcontigs = Slider(title="Maximum Number of Contigs", value=max(parsed['contigs']), start=1, end=max(parsed['contigs']), step=50)
x_axis = Select(title="X Axis", options=sorted(axis_map.keys()), value="Assembly Length")
y_axis = Select(title="Y Axis", options=sorted(axis_map.keys()), value="Checkm Completeness")

# Create Column Data Source that will be used by the plot
source = ColumnDataSource(data=dict(x=[], y=[], color=[], db=[], mmpid=[], apt=[], alpha=[], label=[], qual=[], comp=[], cont=[]))

TOOLTIPS=[
    ("Database", "@db"),
    ("ID", "@mmpid"),
    ("Type", "@apt"),
    ("Quality", "@qual"),
    ("Completeness", "@comp"),
    ("Contamination", "@cont")
]

p = figure(plot_height=800, plot_width=1500, title="", toolbar_location=None, tooltips=TOOLTIPS, sizing_mode="fixed", tools="tap")
p.circle(x="x", y="y", source=source, size=7, color="color", line_color=None, fill_alpha="alpha", legend='label')


def select_entries():
    db_val = database.value
    type_val = analysis_type.value
    qual_val = draft.value

    selected = parsed[
        (parsed.QS >= minqsscore.value) &
        (parsed.QS <= maxqsscore.value) &
        (parsed.contigs >= mincontigs.value) &
        (parsed.contigs <= maxcontigs.value) &
        (parsed.Assembly_length >= minlength.value) &
        (parsed.Assembly_length <= maxlength.value)
    ]
    
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

    p.xaxis.axis_label = x_axis.value
    p.yaxis.axis_label = y_axis.value
    p.title.text = "%d entries shown" % len(df)
    p.title.text = "{} entries shown, {} entries have missing data (not shown)".format(len(df), missing)
    url = "https://mmp.sfb.uit.no/databases/@db/#/records/@mmpid"
    taptool = p.select(type=TapTool)
    taptool.callback = OpenURL(url=url)
    source.data = dict(
        x=df[x_name],
        y=df[y_name],
        color=df["color"],
        db=df["db"],
        mmpid=df["mmp_ID"],
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

update()  # initial load of the data

curdoc().add_root(l)
curdoc().title = "MMP"
