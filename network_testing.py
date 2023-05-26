import base64
import io
import pandas as pd
import flask
import dash
import dash_dangerously_set_inner_html as dhtml
from dash import dcc
from dash import html
from dash.dependencies import Input, Output, State
import plotly.express as px
import plotly.graph_objects as go
from rdkit import Chem
from rdkit.Chem.Draw import rdDepictor
import networkx as nx
from PIL import Image
import tempfile
import os
import atexit

def smi2svg(smi):
    mol = Chem.MolFromSmiles(smi)
    rdDepictor.Compute2DCoords(mol)
    mc = Chem.Mol(mol.ToBinary())
    Chem.Kekulize(mc)
    drawer = Chem.Draw.MolDraw2DSVG(200, 200)
    drawer.DrawMolecule(mc)
    drawer.FinishDrawing()
    svg = drawer.GetDrawingText().replace('svg:', '')
    return svg

upload_style = {
    "width": "20%",
    "height": "120px",
    "lineHeight": "60px",
    "borderWidth": "1px",
    "borderStyle": "dashed",
    "borderRadius": "5px",
    "textAlign": "center",
    "margin": "10px",
    "margin": "3% auto",
}

external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']
app = dash.Dash(__name__, external_stylesheets=external_stylesheets)
app.config.suppress_callback_exceptions = True

vals = {'energy': 'energy'}

app.layout = html.Div(children=[
    html.H1(children='CTYviewer'),
    html.Div([
        html.Div([
            html.Div(id="molimg"),
            dcc.Graph(id='mol_graph', style={'width': '100%', 'height': '800px'})
        ], className="six columns"),
        html.Div([
            dcc.Upload(
                id='csv',
                children=html.Div(['Upload CSV']),
                style=upload_style,
            ),
            html.Div(id='molimg-details')
        ], className="six columns")
    ], className="row"),
    html.Div([
        dcc.Dropdown(
            id='y-column',
            value='energy',
            options=[{'label': key, 'value': key} for key in vals.keys()],
            style={'width': '20%'}
        ),
    ], className="row"),
])



def parse_content(contents, filename):
    content_type, content_string = contents.split(",")
    decoded = base64.b64decode(content_string)
    bio = io.BytesIO(decoded)
    bio.seek(0)
    try:
        df = pd.read_csv(bio)
    except Exception as e:
        print(e)
        return html.Div([f"{filename} error occurred during file reading"])
    df.columns = df.columns.str.strip()
    return df


@app.callback(
    [Output('mol_graph', 'figure'), Output('molimg', 'children')],
    [Input('csv', 'contents'), Input('y-column', 'value')],
    [State('csv', 'filename')],
    prevent_initial_call=True
)
def update_graph_and_img(contents, y_column, filename):
    ctx = dash.callback_context
    trigger_id = ctx.triggered[0]['prop_id'].split('.')[0]

    if trigger_id == 'csv':
        df = parse_content(contents, filename)
        df['index'] = range(len(df))

        G = nx.from_pandas_edgelist(df, source='index', target='connectivity', create_using=nx.Graph())

        pos = nx.spring_layout(G, k=0.3)  # Layout algorithm to determine node positions

        node_trace = go.Scatter(
            x=[pos[node][0] for node in G.nodes()],
            y=[pos[node][1] for node in G.nodes()],
            mode='markers',
            marker=dict(
                size=15,
                opacity=0.5
            ),
            text=[f"Index: {index}<br>Connectivity: {row['connectivity']}<br>Energy: {row[y_column]}" for index, row in df.iterrows()],
            hoverinfo='text'
        )

        # Calculate the line width based on the absolute energy value
        line_width = [max(1, abs(row['energy'])) for _, row in df.iterrows()]

        min_width = min(line_width)
        max_width = max(line_width)
        normalized_width = [int((width - min_width) / (max_width - min_width) * 10) for width in line_width]



        # Create separate traces for each line segment
        edge_traces = []
        for edge, width in zip(G.edges(), normalized_width):
            edge_trace = go.Scatter(
                x=[pos[edge[0]][0], pos[edge[1]][0], None],
                y=[pos[edge[0]][1], pos[edge[1]][1], None],
                mode='lines',
                line=dict(width=width, color='rgb(0, 0, 0)'),
                hoverinfo='none'
            )
            edge_traces.append(edge_trace)
            
        # Create a single trace for all line segments
        all_edges_trace = go.Scatter(
            x=[pos[edge[0]][0] for edge in G.edges()] + [pos[edge[1]][0] for edge in G.edges()] + [None],
            y=[pos[edge[0]][1] for edge in G.edges()] + [pos[edge[1]][1] for edge in G.edges()] + [None],
            mode='lines',
            line=dict(width=1, color='rgb(0, 0, 0)'),  # Set a default width for the remaining lines
            hoverinfo='none'
        )





        fig = go.Figure(data=[node_trace] + edge_traces + [all_edges_trace],
                        layout=go.Layout(
                            showlegend=False,
                            hovermode='closest',
                            xaxis={'title': 'X-axis'},
                            yaxis={'title': 'Y-axis'},
                            title='Network Visualization'
                        ))

        for index, row in df.iterrows():
            svg = smi2svg(row['smiles'])
            with open(f"temp_{index}.svg", "w") as f:
                f.write(svg)
            cmd = f'cairosvg -f png -o temp_{index}.png temp_{index}.svg'
            os.system(cmd)

        return fig, dash.no_update

    return dash.no_update, dash.no_update


@app.callback(
    Output('molimg-details', 'children'),
    [Input('mol_graph', 'hoverData')],
)
def update_img(hoverData):
    if hoverData is None:
        svg = '2D Image'
    else:
        try:
            index = hoverData['points'][0]['text'].split("<br>")[0].split(":")[1].strip()
            energy = hoverData['points'][0]['text'].split("<br>")[2].split(":")[1].strip()
            connectivity = hoverData['points'][0]['text'].split("<br>")[1].split(":")[1].strip()

            image = Image.open(f"temp_{index}.png")
            with tempfile.NamedTemporaryFile(suffix=".png", delete=False) as temp_file:
                temp_filename = temp_file.name
                image.save(temp_filename)

            with open(temp_filename, 'rb') as image_file:
                encoded_image = base64.b64encode(image_file.read()).decode('ascii')

            svg = f'<div style="display: flex; align-items: center;">' \
                  f'<img src="data:image/png;base64,{encoded_image}" style="width:200px;height:200px;">' \
                  f'<div style="margin-left: 10px;">' \
                  f'<p><strong>Index:</strong> {index}</p>' \
                  f'<p><strong>Energy:</strong> {energy}</p>' \
                  f'<p><strong>Connectivity:</strong> {connectivity}</p>' \
                  f'</div>' \
                  f'</div>'
        except:
            svg = '2D Image'

    return dhtml.DangerouslySetInnerHTML(svg)


if __name__ == "__main__":
    app.run_server(debug=True)
