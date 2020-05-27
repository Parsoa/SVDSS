import io
import os
import re
import pwd
import sys
import copy
import json
import math
import time
import random
import argparse

import plotly.offline as plotly
import plotly.graph_objs as graph_objs
import plotly.figure_factory as ff

# ============================================================================================================================ #
# Plotly helpers
# ============================================================================================================================ #

def histogram(x, name, path, x_label, y_label, step = 1):
    data = [graph_objs.Histogram(x = _x, xbins = dict(start = min(_x), size = step, end = max(_x) + 1)) for _x in x]
    layout = graph_objs.Layout(title = name, xaxis = dict(title = x_label), yaxis = dict(title = y_label))
    figure = graph_objs.Figure(data = data, layout = layout)
    plotly.plot(figure, filename = os.path.join(path, 'histogram_' + name + '.html'), auto_open = False)
    pass

def scatter(x, y, name, path, x_label, y_label):
    data = [graph_objs.Scatter(x = x, y = y, mode = 'lines')]
    layout = graph_objs.Layout(title = name, xaxis = dict(title = x_label), yaxis = dict(title = y_label))
    figure = graph_objs.Figure(data = data, layout = layout)
    plotly.plot(figure, filename = os.path.join(path, 'scatter_' + name + '.html'), auto_open = False)
    pass

def bar(x, ys, name, path, x_label, y_label):
    #data = []
    #for y in ys:
    #    data.append(graph_objs.Bar(x = x, y = y))
    #layout = graph_objs.Layout(title = name, xaxis = dict(title = x_label), yaxis = dict(title = y_label))
    #figure = graph_objs.Figure(data = data, layout = layout)
    #plotly.plot(figure, filename = os.path.join(path, 'bar_' + name + '.html'), auto_open = False)
    pass

def violin(x, y, name, path, x_label, y_label):
    data = []
    s = sorted(set(x))
    for i in s:
        trace = {
                "x": [x[j] for j, _ in enumerate(y) if x[j] == i],
                "y": [y[j] for j, _ in enumerate(y) if x[j] == i],
                "name": i,
                "type": 'violin',
                "box": {
                    "visible": True
                },
                "meanline": {
                    "visible": True
                }
            }
        data.append(trace)
    fig = {
        "data": data,
        "layout" : {
            "title": "",
            "yaxis": {
                "zeroline": False,
            }
        }
    }
    #plotly.plot(fig, filename = os.path.join(path, 'violin_' + name + '.html'), auto_open = False)
    pass

