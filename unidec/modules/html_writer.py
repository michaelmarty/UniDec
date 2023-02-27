import mpld3
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors
# import plotly.tools as tls
# from xml.etree import ElementTree as ET
import numpy as np
import lxml.etree as ET
from lxml.html.clean import Cleaner
from io import StringIO, BytesIO
import io
import unidec.tools as ud
import base64
import webbrowser
import os

luminance_cutoff = 135


def write_to_html(html_str, outfile, mode="a"):
    # print(outfile)
    html_file = io.open(outfile, mode, encoding='utf-8')
    html_file.write(html_str)
    html_file.close()


def fig_to_html_mpld3(fig, outfile):
    html_str = mpld3.fig_to_html(fig, no_extras=True)
    write_to_html(html_str, outfile)

'''
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
    return html_str'''


def wrap_to_grid(inputlist, outfile=None):
    # Wrap a list of strings in a grid
    grid = ET.Element("div")
    grid.set("class", "grid-container")
    for i, row in enumerate(inputlist):
        rowdiv = ET.SubElement(grid, "div")
        rowdiv.set("class", "row")
        for j, item in enumerate(row):
            div = ET.Element("div")
            div.set("class", "column")
            # div.text = item

            parser = ET.HTMLParser()
            try:
                item = item.encode()
            except Exception:
                pass
            xmlitem = ET.parse(BytesIO(item), parser)
            xmlelement = xmlitem.getroot().find("body")
            if xmlelement is None:
                xmlelement = xmlitem.getroot()
            rowdiv.append(xmlelement)
        grid.append(rowdiv)
    grid_str = ET.tostring(grid, encoding='unicode')

    if outfile is not None:
        write_to_html(grid_str, outfile)
    return grid_str


def array_to_html(array, outfile=None, cols=None, rows=None, colors=None, index=True):
    df = pd.DataFrame(array, columns=cols, index=rows)
    return df_to_html(df, outfile, colors=colors, index=index)


def df_to_html(df, outfile=None, colors=None, index=True):
    html_str = df.to_html(index=index)
    if colors is not None:
        for i, color in enumerate(colors):
            hexcolor = matplotlib.colors.to_hex(color)
            rgbcolor = matplotlib.colors.to_rgb(hexcolor)
            # print(rgbcolor)
            try:
                luminance = ud.get_luminance(np.array(rgbcolor) * 255, type=2)
            except Exception:
                luminance = 255

            if luminance < luminance_cutoff:
                textcolor = 'white'
            else:
                textcolor = 'black'
            # print("Colors:", color, luminance, textcolor)
            html_str = html_str.replace('<tr>', '<tr style="background-color: %s; color: %s">'
                                        % (hexcolor, textcolor), 1)
    if outfile is not None:
        write_to_html(html_str, outfile)
    return html_str


def html_title(outtitle, outfile=None):
    # CSS styling
    style = ET.Element('style')
    style.text = "header {background-color: #0C234B;}\n"
    style.text += "h1 {color: #e8a219; text-align:left; margin:0; padding:10px}\n"
    style.text += "h2 {color: #AB0520; text-align:left; margin:0; padding:10px}\n"
    style.text += "h3 {color: #AB0520; text-align:left; margin:0; padding:10px}\n"
    style.text += "body {margin:0; padding:0; font-family: \"Helvetica Neue\", Helvetica, Arial, sans-serif;}"
    style.text += "table {border-collapse: collapse; margin:25px; padding:0}\n"
    style.text += "th {text-align:left; background-color:#ADD8E6;; color:black;}\n"
    style.text += "tr:nth-child(even) {background-color: #f2f2f2;}\n"
    style.text += ".grid-container {display:grid; margin:25px;} \n"
    style.text += ".row {display:flex;} \n"
    style_str = ET.tostring(style, encoding='unicode')
    html_str = style_str

    try:
        # Split file name title into file and directory
        dirname = os.path.split(outtitle)[0]
        fname = os.path.split(outtitle)[1]
    except Exception as e:
        print("Error splitting file name:", e)
        dirname = ""
        fname = outtitle

    # Header
    head = ET.Element("head")
    title = ET.Element("title")
    title.text = "UniDec Report: " + str(fname)
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
    h2.text = "File Name: " + str(fname)
    h3 = ET.Element("h3")
    h3.text = "Directory: " + str(dirname)
    header.append(h2)
    header.append(h3)
    body.append(header)
    bodystring = ET.tostring(body, encoding='unicode')
    html_str += bodystring

    if outfile is not None:
        write_to_html(html_str, outfile)
    return html_str


# Function to create HTML collapsible from text
def to_html_collapsible(text, title="Click to expand", canopen=True, outfile=None, htmltext=False):
    """
    Create HTML collapsible from text.
    :param text: Text to be displayed in collapsible
    :param title: Title of collapsible
    :param canopen: Whether collapsible can be opened or not
    :param outfile: File output
    :param htmltext:
    :return:
    """
    # CSS styling
    style = ET.Element('style')
    style.text = ".collapsible {background-color: #0C234B; color: #e8a219; cursor: pointer; " \
                 "padding: 18px; width: 100%; border: none; text-align: left; outline: none; font-size: 15px;}\n"
    style.text += ".active, .collapsible:hover {background-color: #AB0520;}\n"
    style.text += ".collapsible:after {content: '+'; color: white; " \
                  "font-weight: bold; float: right; margin-left: 5px;}\n"
    style.text += ".active:after {content: '-';}\n"
    style.text += ".content {padding: 0 18px; display: none; overflow: hidden; background-color: #f1f1f1;}\n"
    style_str = ET.tostring(style, encoding='unicode')
    html_str = style_str

    # Body
    body = ET.Element("body")
    button = ET.Element("button")
    button.set("class", "collapsible")
    button.text = title
    body.append(button)
    content = ET.Element("div")
    content.set("class", "content")
    if not htmltext:
        p = ET.Element("p")
        p.text = text
        content.append(p)
    else:
        # Parse incoming HTML table from text
        parser = ET.HTMLParser()
        try:
            text = text.encode()
        except Exception:
            pass
        xmlitem = ET.parse(BytesIO(text), parser)
        content.append(xmlitem.getroot())
    body.append(content)
    bodystring = ET.tostring(body, encoding='unicode')
    html_str += bodystring

    if canopen:
        html_str += "<script>var coll = document.getElementsByClassName(\"collapsible\");\n"
        html_str += "var i;\n"
        html_str += "for (i = 0; i < coll.length; i++) {\n"
        html_str += "coll[i].addEventListener(\"click\", function() {\n"
        html_str += "this.classList.toggle(\"active\");\n"
        html_str += "var content = this.nextElementSibling;\n"
        html_str += "if (content.style.display === \"block\") {\n"
        html_str += "content.style.display = \"none\";\n"
        html_str += "} else {\n"
        html_str += "content.style.display = \"block\";\n"
        html_str += "}\n"
        html_str += "});\n"
        html_str += "}\n"
        html_str += "</script>"

    if outfile:
        write_to_html(html_str, outfile)
    return html_str


# Function to create an html table from a python dictionary
def dict_to_html(indict, outfile=None):
    html_str = "<table>\n"
    for key, value in indict.items():
        html_str += "<tr><td>" + str(key) + "</td><td>" + str(value) + "</td></tr>\n"
    html_str += "</table>\n"

    if outfile is not None:
        write_to_html(html_str, outfile)
    return html_str


# Function to embed png image from bite string in html
def png_to_html(png_str, outfile=None):
    png_str = base64.b64encode(png_str)
    png_str = png_str.decode("utf-8")
    html_str = "<img src=\"data:image/png;base64," + png_str + "\" \"/>\n"

    if outfile is not None:
        write_to_html(html_str, outfile)
    return html_str


def html_open(outfile):
    html_file = open(outfile, "w")
    html_file.write("<html lang=\"en\">\n")
    html_file.close()


def html_close(outfile):
    html_str = "</html>\n"
    write_to_html(html_str, outfile)
    html_cleaner(outfile)


def html_cleaner(outfile):
    cleaner = Cleaner(scripts=False, javascript=False, style=False, inline_style=False, page_structure=False,
                      meta=False, embedded=False, links=False, comments=False, frames=False, forms=False,
                      safe_attrs_only=False, processing_instructions=False, annoying_tags=False,
                      remove_unknown_tags=False, )

    html_file = io.open(outfile, "r", encoding='utf-8')
    html_str = html_file.read()
    html_file.close()

    clean_html = cleaner.clean_html(html_str)

    write_to_html(clean_html, outfile, mode="w")


if __name__ == "__main__":
    path = "C:\\Python\\UniDec3\\unidec\\bin\\Example Data\\BSA_unidecfiles"
    import os

    os.chdir(path)

    figure = plt.figure()
    ax = plt.plot([1, 2, 3, 4, 5], [2, 5, 6, 3, 7])
    # plt.show()
    outfile_html = "test.html"
    html_open(outfile_html)
    html_title("Test File", outfile_html)

    colorslist = ["red", "green", "blue", "yellow", "orange"]

    s2 = array_to_html(np.random.random((5, 5)), outfile_html, cols=["a", "b", "c", "d", "e"], colors=colorslist)
    # s1 = fig_to_html_plotly(fig, outfile)
    # s1 = fig_to_html_plotly(fig, outfile)
    # wrap_to_grid([s2, s1], outfile)

    png = io.BytesIO()
    figure.savefig(png, format="png")
    png_str2 = png.getvalue()
    png_to_html(png_str2, outfile_html)

    dict_string = dict_to_html({"a": 1, "b": 2, "c": 3, "d": 4, "e": 5})
    to_html_collapsible(dict_string, title="UniDec Parameters", outfile=outfile_html, htmltext=True)

    html_close(outfile_html)

    html_cleaner(outfile_html)

    #opencommand = "start \"\" "
    #os.system(opencommand + "\"" + outfile_html + "\"")
    webbrowser.open(outfile_html)
