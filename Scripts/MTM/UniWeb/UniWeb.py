# -*- coding: utf-8 -*-
import dash
import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Input, Output, State
import base64, io
import numpy as np

external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']

app = dash.Dash(__name__, external_stylesheets=external_stylesheets)

app.layout = html.Div(children=[
    html.H1(children='unidec Web'),

    html.Div(children='''
        unidec Web App Testing
    '''),
    dcc.Upload(
        id='upload-data',
        children=html.Div([
            'Drag and Drop or ',
            html.A('Select Files')
        ]),
        style={
            'width': '100%',
            'height': '60px',
            'lineHeight': '60px',
            'borderWidth': '1px',
            'borderStyle': 'dashed',
            'borderRadius': '5px',
            'textAlign': 'center',
            'margin': '10px'
        },
        # Allow multiple files to be uploaded
        multiple=True
    ),
    html.Div(id='output-data-upload'),

])


def parse_contents(contents, filename, date):
    print(filename)
    content_type, content_string = contents.split(',')

    decoded = base64.b64decode(content_string)
    des = io.StringIO(decoded.decode('utf-8'))
    data = np.loadtxt(des)
    x = data[::2]
    y = data[1::2]
    data = np.transpose([x,y])
    #for l in des.split("\r"):
        #print(l)
    #    pass
    output = make_graph(data)
    return output



@app.callback(Output('output-data-upload', 'children'),
              [Input('upload-data', 'contents')],
              [State('upload-data', 'filename'),
               State('upload-data', 'last_modified')])

def update_output(list_of_contents, list_of_names, list_of_dates):
    if list_of_contents is not None:
        children = [
            parse_contents(c, n, d) for c, n, d in
            zip(list_of_contents, list_of_names, list_of_dates)]
        return children


def make_graph(data):
    output = dcc.Graph(
        id='example-graph',
        figure={
            'data': [
                {'x': data[:,0], 'y': data[:,1], 'mode': 'lines'},
            ],
            'layout': {
                'title': 'Mass Spectrum'
            }
        }
    )
    return output

if __name__ == '__main__':
    app.run_server(debug=True)