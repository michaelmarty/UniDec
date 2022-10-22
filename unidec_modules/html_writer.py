import mpld3
import pandas as pd
import matplotlib.pyplot as plt
import plotly.tools as tls


def write_to_html(html_str, outfile):
    print(outfile)
    Html_file = open(outfile, "a")
    Html_file.write(html_str)
    Html_file.close()


def fig_to_html(fig, outfile):
    html_str = mpld3.fig_to_html(fig, no_extras=True)
    write_to_html(html_str, outfile)


def fig_to_html_plotly(fig, outfile):
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
                     #mirror='allticks',  # add ticks to top/right axes
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
                     #mirror='allticks',
                     tickwidth=2.4,
                     tickcolor='black',
                     )

    html_str = plotly_fig.to_html(full_html=False)
    write_to_html(html_str, outfile)


def array_to_html(array, outfile, cols=None, rows=None):
    df = pd.DataFrame(array, columns=cols, index=rows)
    html_str = df.to_html()
    write_to_html(html_str, outfile)


def html_title(title, outfile):
    html_str = '<head>\n<title>\n' + str(title) + '\n</title>\n</head>\n'
    html_str += '<body>\n <h1>\n UniDec Report \n</h1>\n'
    html_str += '<h2>\nFile Name: ' + str(title) + '\n</h2>\n'

    write_to_html(html_str, outfile)


def html_open(outfile):
    Html_file = open(outfile, "w")
    Html_file.write("<html>\n")
    Html_file.close()


def html_close(outfile):
    html_str = "</body>\n</html>\n"
    write_to_html(html_str, outfile)


if __name__ == "__main__":
    path = "C:\\Python\\UniDec3\\unidec_bin\\Example Data\\BSA_unidecfiles"
    import os

    os.chdir(path)

    fig = plt.figure()
    ax = plt.plot([1, 2, 3, 4, 5], [2, 5, 6, 3, 7])
    # plt.show()

    fig_to_html_plotly(fig, "test.html")
