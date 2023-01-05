# import mpld3
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors
import plotly.tools as tls
# from xml.etree import ElementTree as ET
import numpy as np
import lxml.etree as ET
from io import StringIO
import unidec.tools as ud

luminance_cutoff = 135


def write_to_html(html_str, outfile):
    print(outfile)
    Html_file = open(outfile, "a")
    Html_file.write(html_str)
    Html_file.close()


'''
def fig_to_html(fig, outfile):
    html_str = mpld3.fig_to_html(fig, no_extras=True)
    write_to_html(html_str, outfile)'''


def fig_to_html_plotly(fig, outfile=None):
    xlabel = fig.gca().get_xlabel()
    ylabel = fig.gca().get_ylabel()
    plotly_fig = tls.mpl_to_plotly(fig)
    # choose the figure font
    font_dict = dict(family='Arial',
                     size=26,
                     color='black'
                     )
    # general figure formatting
    plotly_fig.update_layout(font=font_dict,  # font formatting
                             plot_bgcolor='white',  # background color
                             width=600,  # figure width
                             height=400,  # figure height
                             margin=dict(r=20, t=20, b=10)  # remove white space
                             )
    # x and y-axis formatting
    plotly_fig.update_yaxes(title_text=ylabel,  # axis label
                            title_font=font_dict,
                            showline=True,  # add line at x=0
                            linecolor='black',  # line color
                            linewidth=2.4,  # line size
                            ticks='outside',  # ticks outside axis
                            tickfont=font_dict,  # tick label font
                            # mirror='allticks',  # add ticks to top/right axes
                            tickwidth=2.4,  # tick width
                            tickcolor='black',  # tick color
                            )
    plotly_fig.update_xaxes(title_text=xlabel,
                            title_font=font_dict,
                            showline=True,
                            showticklabels=True,
                            linecolor='black',
                            linewidth=2.4,
                            ticks='outside',
                            tickfont=font_dict,
                            # mirror='allticks',
                            tickwidth=2.4,
                            tickcolor='black',
                            )

    html_str = plotly_fig.to_html(full_html=False)
    if outfile is not None:
        write_to_html(html_str, outfile)
    return html_str

'''
def wrap_to_grid(inputlist, outfile=None):
    # Wrap a list of strings in a grid
    grid = ET.Element("div")
    grid.set("class", "row")
    for i, item in enumerate(inputlist):
        div = ET.Element("div")
        div.set("class", "column")
        parser = ET.HTMLParser()
        xmlitem = ET.parse(StringIO(item), parser)
        xmlelement = xmlitem.getroot().find("body").find("div")  # .getchildren()[0].getchildren()[0]
        if xmlelement is None:
            xmlelement = xmlitem.getroot()
        div.append(xmlelement)
        grid.append(div)
    grid_str = ET.tostring(grid, encoding='unicode')
    print(grid_str)
    if outfile is not None:
        write_to_html(grid_str, outfile)
    return grid_str


# write a function to strip out the <div> tags from the html string
def strip_divs(html_str):
    # strip out the div tags
    html_str = html_str.replace('<div>', '')
    html_str = html_str.replace('</div>', '')
    return html_str
'''

def array_to_html(array, outfile=None, cols=None, rows=None, colors=None):
    df = pd.DataFrame(array, columns=cols, index=rows)
    return df_to_html(df, outfile, colors=colors)


def df_to_html(df, outfile=None, colors=None):
    html_str = df.to_html()
    if colors is not None:
        for i, color in enumerate(colors):
            hexcolor = matplotlib.colors.to_hex(color)
            try:
                luminance = ud.get_luminance(color*255, type=2)
            except:
                luminance = 255

            if luminance < luminance_cutoff:
                textcolor = 'white'
            else:
                textcolor = 'black'
            print("Colors:", color, luminance, textcolor)
            html_str = html_str.replace('<tr>', '<tr style="background-color: %s; color: %s">'
                                        % (hexcolor, textcolor), 1)
    if outfile is not None:
        write_to_html(html_str, outfile)
    return html_str


def html_title(outtitle, outfile=None):
    # CSS styling
    style = ET.Element('style')
    style.text = "header {background-color: #0C234B;}\n"
    style.text += "h1 {color: #e8a219; text-align:left; margin:0; padding:10}\n"
    style.text += "h2 {color: #AB0520; text-align:left; margin:0; padding:10}\n"
    style.text += "body {margin:0; padding:0; font-family: \"Helvetica Neue\", Helvetica, Arial, sans-serif;}"
    style.text += "table {border-collapse: collapse; margin:25px; padding:0}\n"
    style.text += "th {text-align:left; background-color:#ADD8E6;; color:black;}\n"
    style.text += "tr:nth-child(even) {background-color: #f2f2f2;}\n"
    style.text += ".grid-container {display:grid;} \n"
    style.text += ".row {display:flex;} \n"
    style_str = ET.tostring(style, encoding='unicode')
    html_str = style_str

    # Header
    head = ET.Element("head")
    title = ET.Element("title")
    title.text = str(outtitle)
    head.append(title)
    headerstring = ET.tostring(head, encoding='unicode')
    html_str += headerstring

    # Body
    body = ET.Element("body")
    header = ET.Element("header")
    h1 = ET.Element("h1")
    h1.text = "UniDec Report"
    header.append(h1)
    h2 = ET.Element("h2")
    h2.text = "File Name: " + str(outtitle)
    header.append(h2)
    body.append(header)
    bodystring = ET.tostring(body, encoding='unicode')
    html_str += bodystring

    if outfile is not None:
        write_to_html(str(html_str), outfile)


def html_open(outfile):
    Html_file = open(outfile, "w")
    Html_file.write("<html>\n")
    Html_file.close()


def html_close(outfile):
    html_str = "</body>\n</html>\n"
    write_to_html(html_str, outfile)


if __name__ == "__main__":
    path = "C:\\Python\\UniDec3\\unidec\\bin\\Example Data\\BSA_unidecfiles"
    import os

    os.chdir(path)

    fig = plt.figure()
    ax = plt.plot([1, 2, 3, 4, 5], [2, 5, 6, 3, 7])
    # plt.show()
    outfile = "test.html"
    html_open(outfile)
    html_title("Test File", outfile)

    colors = ["red", "green", "blue", "yellow", "orange"]

    s2 = array_to_html(np.random.random((5, 5)), outfile, cols=["a", "b", "c", "d", "e"], colors=colors)
    s1 = fig_to_html_plotly(fig, outfile)
    # s1 = fig_to_html_plotly(fig, outfile)

    # wrap_to_grid([s2, s1], outfile)

    html_close(outfile)

    opencommand = "start \"\" "
    os.system(opencommand + "\"" + outfile + "\"")
