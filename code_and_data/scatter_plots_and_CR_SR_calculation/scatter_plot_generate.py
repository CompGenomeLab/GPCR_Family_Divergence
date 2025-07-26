# -*- coding: utf-8 -*-
"""
Created on Wed Dec 20 14:47:08 2023

@author: selcuk.1
"""

import pandas as pd
import plotly.express as px
import numpy as np


def distribute_overlapping_points(df, cons_col='conservation', ent_col='entropy', annot_col='annotation', cons_adjust=2.5, ent_adjust=0.06):
    # Remove duplicate (conservation, entropy, annotation) rows
    df = df.drop_duplicates(subset=[cons_col, ent_col, annot_col])

    # Identify duplicated (conservation, entropy) pairs
    duplicated = df[df.duplicated(subset=[cons_col, ent_col], keep=False)]

    # Find unique (conservation, entropy) pairs within the duplicated points
    unique_points = duplicated.groupby([cons_col, ent_col]).size().reset_index(name='count')

    # Create a list to store adjusted points
    adjusted_points = []

    for _, row in unique_points.iterrows():
        cons_val = row[cons_col]
        ent_val = row[ent_col]
        count = row['count']

        # Get the corresponding annotations for the duplicated points
        annotations = df[(df[cons_col] == cons_val) & (df[ent_col] == ent_val)][annot_col].values

        if count == 1:
            adjusted_points.append((cons_val, ent_val, annotations[0]))
        else:
            # Distribute points around the original point in a circular manner
            angle_increment = 2 * np.pi / count
            for i in range(int(count)):
                angle = i * angle_increment
                new_cons_val = cons_val + cons_adjust * np.cos(angle)
                new_ent_val = ent_val + ent_adjust * np.sin(angle)
                adjusted_points.append((new_cons_val, new_ent_val, annotations[i]))

    # Create a new dataframe with adjusted points for duplicated pairs
    adjusted_df = pd.DataFrame(adjusted_points, columns=[cons_col, ent_col, annot_col])

    # Replace the duplicated rows in the original dataframe with the adjusted points
    df.loc[duplicated.index, [cons_col, ent_col, annot_col]] = adjusted_df.values

    return df

import os

# Get the absolute path of the directory where this script resides
script_dir = os.path.dirname(os.path.abspath(__file__))

# Change the current working directory to that directory
os.chdir(script_dir)

#%%
scatter_data="classA_labels.txt"
cons_thr=75.27
ent_thr1=1.78
ent_thr2=2.63
ordered_df = pd.read_csv(scatter_data,sep="\t")
out_figure="classA_scatter.svg"

unique_annotations = ordered_df['annotation'].unique()
print(unique_annotations)
# List of values to change to "Other"
values_to_change = ["CWxP","NPxxY","PIF","DRY","Na+ Pocket" ]

# Modify the 'annotation' column
ordered_df['annotation'] = ordered_df['annotation'].apply(lambda x: 'Known Motifs' if x in values_to_change else x)

values_to_change = ["Common Activation Mechanism"]
ordered_df['annotation'] = ordered_df['annotation'].apply(lambda x: 'Other' if x in values_to_change else x)
unique_annotations = ordered_df['annotation'].unique()
print(unique_annotations)
#%%
scatter_data="classOlf_labels.txt"
cons_thr=70
ent_thr1=1
ent_thr2=2.62
ordered_df = pd.read_csv(scatter_data,sep="\t")
out_figure="classOlf_scatter.svg"
#%%
scatter_data="classT_labels.txt"
cons_thr=31.25
ent_thr1=1
ent_thr2=2.4
ordered_df = pd.read_csv(scatter_data,sep="\t")
ordered_df=distribute_overlapping_points(ordered_df)
out_figure="classT_scatter.svg"
#%%
scatter_data="classB1_labels.txt"
cons_thr=73.33
ent_thr1=0.92
ent_thr2=2
ordered_df=pd.read_csv(scatter_data,sep="\t")
ordered_df=distribute_overlapping_points(ordered_df)
out_figure="classB1_scatter.svg"
#%%
scatter_data="classB2_labels.txt"
cons_thr=71.88
ent_thr1=0.93
ent_thr2=2.5
ordered_df=pd.read_csv(scatter_data,sep="\t")
ordered_df=distribute_overlapping_points(ordered_df)
out_figure="classB2_scatter.svg"
#%%
scatter_data="classF_labels.txt"
cons_thr=72.73
ent_thr1=1.1
ent_thr2=1.79
ordered_df=pd.read_csv(scatter_data,sep="\t")
ordered_df=distribute_overlapping_points(ordered_df)
out_figure="classF_scatter.svg"
#%%
scatter_data="classC_labels.txt"
cons_thr=73.33
ent_thr1=1.05
ent_thr2=2.13
ordered_df=pd.read_csv(scatter_data,sep="\t")
ordered_df=distribute_overlapping_points(ordered_df,cons_adjust=1)
out_figure="classC_scatter.svg"
#%%
# Map annotations to colors
color_map = {}
color_map["Transducer"] = "#fdb415"
color_map['Disulfide Bridge'] = "#df6131"
color_map['Cholesterol'] = "#7baec4"
color_map['Known Motifs'] = "#131110"
color_map['ICL2'] = "#839d28"
color_map['Tethered Agonist'] = "#839d28"
color_map['VFT Ligand'] = "#839d28"
color_map['Ligand'] = "#7d2a86"
color_map['Allosteric Modulator'] = "#7d2a86"
color_map["WNT Binding"] = "#7d2a86"
color_map['Other'] = "#9DA39A"

fig = px.scatter(ordered_df, y="conservation", x="entropy",color="annotation"
,template="simple_white",color_discrete_map=color_map, render_mode = 'svg')

fig.update_traces(marker_size=5)
fig.update_traces(textposition='top center')
fig.update_layout(legend=dict(title='Annotation'))
axis_font_dict = dict(family="Arial", size=12, color="black")
fig.update_xaxes(title_font=axis_font_dict, tickfont=axis_font_dict)
fig.update_yaxes(title_font=axis_font_dict, tickfont=axis_font_dict)

fig.update_layout( 
    showlegend=False,
    xaxis_title=dict(
        text="Family-Wide Entropy",
        font=dict(family="Arial", size=12, color="black", weight="bold"),standoff=5
    ),
    yaxis=dict(nticks=6),
    yaxis_title=dict(
        text="Family-Wide Conservation",
        font=dict(family="Arial", size=12, color="black", weight="bold"),standoff=5
    )
    
)
fig.update_traces(selector=dict(marker_color= "#A5AA99"))

fig.add_shape(
    type="line",
    x0=ent_thr2, x1=ent_thr2, y0=ordered_df['conservation'].min(), y1=ordered_df['conservation'].max()+10,
    xref="x", yref="y",
    line=dict(color="Black", width=2, dash="dash")
)

# Add a vertical dashed line at x=2
fig.add_shape(
    type="line",
    x0=ent_thr1 , x1=ent_thr1 , y0=ordered_df['conservation'].min(), y1=ordered_df['conservation'].max()+10,
    xref="x", yref="y",
    line=dict(color="Black", width=2, dash="dash")
)

# Add a horizontal dashed line at y=0.5
fig.add_shape(
    type="line",
    x0=ordered_df['entropy'].min(), x1=ordered_df['entropy'].max(), y0=cons_thr, y1=cons_thr,
    xref="x", yref="y",
    line=dict(color="Black", width=2, dash="dash")
)
fig.update_layout(
    legend=dict( title=None,
        x=0,          # Horizontal position (0 to 1, where 0 is left and 1 is right)
        y=0.65,          # Vertical position (0 to 1, where 0 is bottom and 1 is top)
        xanchor='left', # Horizontal anchor
        yanchor='top',  # Vertical anchor
        traceorder='normal' # Default order (top to bottom)
    )
)
# legend_order = ['Known Motifs', 'Ligand', 'G protein', 'Cholesterol','Disulfide Bridge', 'Other']
# for trace in fig.data:
#     if trace.name in legend_order:
#         trace.legendgroup = trace.name
#         trace.showlegend = True
#     else:
#         trace.showlegend = False

# # Update traces to reorder legends
# fig.data = sorted(fig.data, key=lambda trace: legend_order.index(trace.name) if trace.name in legend_order else -1)

fig.update_layout(
    width=400,
    height=300,
    legend=dict(
            font=dict(size=10),
            bgcolor='rgba(0,0,0,0)',
            itemwidth=30),
    yaxis=dict(nticks=6,tickvals=[0,20,40,60,80,100]),
)

fig.update_traces(selector=dict(marker_color="#9DA39A"), marker_opacity=0.5)
fig.write_image(out_figure, format="svg")