from unidec.modules import unidecstructure, peakstructure, plot1d, plot2d
from unidec import tools as ud
import numpy as np
import os
import time
from unidec.modules.html_writer import *

version = "6.0.0.b1"


def copy_config(config):
    # return deepcopy(config)
    return type("UniDecConfig", (object,), dict(config.__dict__))


class UniDecEngine:
    def __init__(self):
        """
        UniDec Engine Base

        Consists of three main subclasses: Config, DataContiner, Peaks
        
        Establishes some core shared functions

        :return: None
        """

        self.version = version

        print("\nUniDec Engine v." + self.version)
        self.config = None
        self.config_history = []
        self.config_count = 0
        self.data = None
        self.pks = None
        self.olg = None
        self.matchcounts = None
        self.altmasses = None
        self.altindexes = None
        self.matchindexes = None
        self.matchlist = None
        self.initialize()

    def initialize(self):
        """
        Initialize Config, DataContainer, and Peaks
        :return: None
        """
        self.config = unidecstructure.UniDecConfig()
        self.clear_history()
        self.config.initialize_system_paths()
        self.reset_config()
        self.data = unidecstructure.DataContainer()
        self.pks = peakstructure.Peaks()
        self.olg = unidecstructure.OligomerContainer()

    def reset_config(self):
        """
        Resets UniDec config to default. Should not reset paths or filenames.
        :return: None
        """
        self.config.initialize()
        self.update_history()

    def load_config(self, f_name):
        """
        Import UniDec Configuration File
        :param f_name: File name
        :return: None
        """
        if f_name is not None:
            self.config.config_import(f_name)
            self.update_history()
        else:
            print("Load Config Error: No file provided.")

    def export_config(self, f_name=None):
        """
        Export UniDec Configuration File
        :param f_name: File name, Default of None will using config.confname
        :return: None
        """
        self.update_history()
        if f_name is not None:
            self.config.config_export(f_name)
        else:
            self.config.config_export(self.config.confname)

    def save_default(self):
        """
        Saves existing config in default location set at self.config.defaultconfig
        :return: None
        """
        self.config.config_export(self.config.defaultconfig)

    def load_default(self):
        """
        Loads config from default location set at self.config.defaultconfig
        :return: None
        """
        self.config.config_import(self.config.defaultconfig)

    def write_hdf5(self):
        self.update_history()
        self.config.write_hdf5(self.config.hdf_file)
        # self.data.write_hdf5(self.config.hdf_file)

    def read_hdf5(self):
        self.config.read_hdf5(self.config.hdf_file)
        self.update_history()
        # self.data.read_hdf5(self.config.hdf_file)

    def update_history(self):
        # print "Update"
        try:
            if self.config_count > 0 and self.config.check_new(self.config_history[len(self.config_history) - 1]):
                self.config_history.append(copy_config(self.config))
                self.config_count = len(self.config_history)
                # print "Updated History", self.config_count
            elif self.config_count == 0:
                self.clear_history()
                # else:
                # print "No changes"
        except Exception as e:
            self.clear_history()
        # print self.config_count
        pass

    def clear_history(self):
        self.config_history = [copy_config(self.config)]
        self.config_count = 1
        pass

    def undo(self):
        if self.config_count > 1:
            self.config_count -= 1
            new = self.config_history[self.config_count - 1]
            old = self.config
            for item in new.__dict__:
                try:
                    old.__dict__[item] = new.__dict__[item]
                except KeyError as e:
                    print(e)
        pass

    def redo(self):
        # print(self.config_count, len(self.config_history))
        if self.config_count < len(self.config_history):
            self.config_count += 1
            new = self.config_history[self.config_count - 1]
            old = self.config
            for item in new.__dict__:
                try:
                    old.__dict__[item] = new.__dict__[item]
                except KeyError as e:
                    print(e)
        pass

    def get_auto_peak_width(self, set_it=True):
        try:
            fwhm, psfun, mid = ud.auto_peak_width(self.data.data2)
            self.config.automzsig = fwhm
            self.config.autopsfun = psfun
            if set_it:
                self.config.psfun = psfun
                self.config.mzsig = fwhm
            print("Automatic Peak Width:", fwhm)
        except Exception as e:
            print("Failed Automatic Peak Width:", e)
            print(self.data.data2)

    def check_badness(self):
        """
        Check for problematic variables, such as upper bounds less than lower bounds and raise warning if found.
        :return:
        """
        badness, warning = self.config.check_badness()
        if warning != "":
            print(warning)
        return badness

    def auto_polarity(self, path=None, importer=None):
        if path is None:
            path = self.config.filename
        self.config.polarity = ud.get_polarity(path, importer=importer)
        if self.config.polarity == "Positive":
            self.config.adductmass = np.abs(self.config.adductmass)
            print("Adduct Mass:", self.config.adductmass)
        elif self.config.polarity == "Negative":
            self.config.adductmass = -1 * np.abs(self.config.adductmass)
            print("Adduct Mass:", self.config.adductmass)

    def linear_regression_peaks(self):
        print("Starting Linear Regression using a repeating mass of:", self.config.molig)
        fit = [0, 0]
        rsquared = 0
        if self.config.molig != 0:
            x = []
            y = []
            z = []
            for p in self.pks.peaks:
                if p.ignore == 0:
                    y.append(p.mass)
                    z.append(p.height)
                    mnum = np.floor(p.mass / self.config.molig)
                    x.append(mnum)

            x = np.array(x)
            y = np.array(y)
            z = np.array(z)

            fit = np.polyfit(x, y, 1, w=z ** 2)
            slope = fit[0]
            intercept = fit[1]

            fitdat = x * slope + intercept

            sse = np.sum((fitdat - y) ** 2)
            denom = np.sum((y - np.mean(y)) ** 2)
            rsquared = 1 - ud.safedivide1(sse, denom)

            print("Slope:", fit[0], "Intercept:", fit[1], "R-Squared:", rsquared)

            residuals = np.abs(fitdat - y)
            cutoff = self.config.molig * 0.1
            boo1 = residuals > cutoff
            if np.any(boo1):
                print("Removing outliers with residuals greater than:", cutoff)
                print(residuals)
                boo2 = residuals < cutoff
                fit = np.polyfit(x[boo2], y[boo2], 1, w=np.array(z[boo2]) ** 2)
                slope = fit[0]
                intercept = fit[1]
                fitdat = x[boo2] * slope + intercept
                sse = np.sum((fitdat - y[boo2]) ** 2)
                denom = np.sum((y[boo2] - np.mean([boo2])) ** 2)
                rsquared = 1 - ud.safedivide1(sse, denom)
                print("New Slope:", fit[0], "New Intercept:", fit[1], "R-Squared:", rsquared)

        else:
            print("Need to set the mass difference/mass of oligomer")
        return fit, rsquared

    def polydispersity_index(self, e=None):
        pdi = ud.polydispersity_index(self.data.massdat)
        print(pdi)

    def oxidation_analysis(self, e=None):
        data = self.data.massdat

        pdata = []
        for p in self.pks.peaks:
            if p.ignore == 0:
                pdata.append([p.mass, p.height])
        pdata = np.array(pdata)
        maxindex = np.argmax(pdata[:, 1])
        maxmass = pdata[maxindex][0]
        nox = np.arange(0, 4)
        oxmasses = maxmass + nox * 16
        areas = []
        for i, m in enumerate(oxmasses):
            low = m + self.config.integratelb
            high = m + self.config.integrateub
            area, intdat = ud.integrate(data, low, high)  # Peak Area
            # area = ud.data_extract(data, (high+low)/2., window=(high-low)/2., extract_method=1) # Peak Height
            areas.append(area)
            print("Peak:", m, "Number of Oxidations:", nox[i], "Range:", low, "to", high, "Local Max:", area)
        areas = np.array(areas)
        try:
            areas /= np.amax(areas)
        except Exception as e:
            pass
        print("Relative Heights:", areas)

        totox = np.sum(nox * areas) / np.sum(areas)
        print("Total Oxidations:", totox)

        return np.append(areas, totox)

    def load_ofile(self, file=None):
        if file is None:
            file = self.config.ofile

        if os.path.isfile(file):
            self.config.oligomerlist = unidecstructure.ofile_reader(file)
            if self.config.oligomerlist.shape == (4,) or self.config.oligomerlist.shape == (5,):
                self.config.oligomerlist = np.array([self.config.oligomerlist])
        else:
            print("Ofile Not Found: ", file)

    def make_oligomers(self, isolated=False, oligomerlist=None, minsites=None, maxsites=None):
        if oligomerlist is None:
            oligomerlist = self.config.oligomerlist
        self.olg.make_oligomers(isolated=isolated, oligomerlist=oligomerlist, minsites=minsites, maxsites=maxsites)

    def match(self, tolerance=100, isolated=False, glyco=False, minsites=None, maxsites=None):
        if ud.isempty(self.olg.oligomasslist):
            self.make_oligomers(minsites=minsites, maxsites=maxsites, isolated=isolated)
        if glyco:
            self.olg.pair_glyco()
        print(self.olg.oligomerlist)
        self.matchlist, self.matchindexes = ud.match(self.pks, self.olg.oligomasslist, self.olg.oligonames,
                                                     self.config.oligomerlist, tolerance=tolerance, return_numbers=True)

    def get_alts(self, tolerance=100):
        self.altmasses, self.altindexes, self.matchcounts = self.olg.get_alts(self.pks, tolerance)

    def get_summed_match_intensities(self, index, alts=True, normmode=0, probarray=None, get_2d=False):
        # To translate the index into oligomer number
        basenumber = int(self.config.oligomerlist[index, 2])

        # If alts are considered, divide up the intensity proporionately into different potential matches
        # Either uses a probability array to tell which is more probable or assumes all are equally likely
        # Note: it would be possible to assume a peak width and set likelihood based on error relative to PW
        if alts:
            # Merge the alternate indexes together into a single giant list
            altindexes = np.concatenate(self.altindexes)
            # Pull out only the ones with the right index
            snumbers = altindexes[:, index]
            sunique = np.unique(snumbers)

            # Allows input of probability array to filter the potential matches to get the most likely
            if probarray is not None:
                # Translate the probability data into indexes
                nvals = probarray[:, 0]
                probs = np.zeros(np.amax(sunique) + 1)
                for i, s in enumerate(sunique):
                    if s + basenumber in nvals:
                        j = np.argwhere(nvals == s + basenumber)[0][0]
                        probs[s] = probarray[j, 1]
                # Create a large array with probabilities matching the indexes
                probs2 = probs[snumbers]

                # Set the intensities based on the peak height and the normalized probability of each peak
                intensities = []
                index = 0
                for i, p in enumerate(self.pks.peaks):
                    n = self.matchcounts[i]
                    subprobs2 = probs2[index:index + n]
                    if np.sum(subprobs2) == 0:
                        print("No matches compatible with probability for peak:", p.mass)
                        intensities.append(np.ones(n) * p.height * subprobs2)
                    else:
                        intensities.append(np.ones(n) * p.height * subprobs2 / np.sum(subprobs2))
                    index += n
                intensities = np.concatenate(intensities)

            else:
                intensities = []
                for i, p in enumerate(self.pks.peaks):
                    intensities.append(np.ones(self.matchcounts[i]) * p.height / self.matchcounts[i])
                intensities = np.concatenate(intensities)

        else:
            # If alts are not considered, just use the matches and assume they have all the peak height
            intensities = np.array([p.height for p in self.pks.peaks])
            snumbers = self.matchindexes[:, index]
            sunique = np.unique(snumbers)

        sint = []
        psint = []
        for i, s in enumerate(sunique):
            b1 = snumbers == s
            subset = intensities[b1]
            sumval = np.sum(subset)
            sint.append(sumval)
            if get_2d:
                index = 0
                peaksums = []
                for j, p in enumerate(self.pks.peaks):
                    mask = np.zeros_like(intensities)
                    if alts:
                        n = self.matchcounts[j]
                        mask[index:index + n] = 1
                        index += n
                    else:
                        mask[j] = 1
                    peaksum = np.sum(intensities[b1 * mask.astype(bool)])
                    peaksums.append(peaksum)

                peaksums = np.array(peaksums)
                # if np.amax(peaksums) != 0:
                #    peaksums = peaksums / np.sum(peaksums)
                psint.append(peaksums)

        sint = np.array(sint)
        if normmode == 1:
            sint /= np.sum(sint)
        if normmode == 2:
            sint /= np.max(sint)
        if get_2d:
            psint = np.array(psint)
            if np.amax(psint) != 0:
                psint /= np.amax(psint)
            return np.transpose([sunique + basenumber, sint]), psint
        else:
            return np.transpose([sunique + basenumber, sint])

    def load_peaks(self, pfile):
        pdata = np.loadtxt(pfile, delimiter=",", usecols=(0, 1))
        self.pks = peakstructure.Peaks()
        self.pks.add_peaks(pdata)

    def makeplot6(self, plot=None, pks=None, show="height", config=None):
        """
        Plots bar chart of peak heights or areas in self.view.plot6.
        :param plot: plot object to use. If none, will create one
        :param pks: peak structure to use. If none, will use self.pks
        :param show: What parameter to plot
        "height" will plot p.height for p in self.eng.pks.peaks
        "integral" will plot p.integral
        :param config: config object to use. If none, will use self.config
        :return: plot object
        """
        if config is None:
            config = self.config
        if pks is None:
            pks = self.pks
        if plot is None:
            plot = plot1d.Plot1dBase()

        if config.batchflag == 0:
            if pks.plen > 0:
                num = 0
                ints = []
                cols = []
                labs = []
                marks = []
                for i, p in enumerate(pks.peaks):
                    if p.ignore == 0:
                        num += 1
                        if show == "height":
                            ints.append(p.height)
                        elif show == "integral":
                            ints.append(p.integral)
                        else:
                            ints.append(0)
                        cols.append(p.color)
                        labs.append(p.label)
                        marks.append(p.marker)
                indexes = list(range(0, num))
                plot.barplottop(indexes, ints, labs, cols, "Species", "Intensity",
                                "Peak Intensities", repaint=False)
                for i in indexes:
                    plot.plotadddot(i, ints[i], cols[i], marks[i])
            plot.repaint()
        return plot

    def makeplot3(self, plot=None, data=None, config=None):
        """
        Plot m/z vs charge grid.
        :param plot: Plot object to use
        :param data: Data to plot
        :param config: Config object to use
        :return: plot object
        """
        if plot is None:
            plot = plot2d.Plot2dBase()
        if config is None:
            config = self.config
        if data is None:
            data = self.data.mzgrid
        if config.batchflag == 0:
            tstart = time.perf_counter()
            plot.contourplot(data, config)
            print("Plot 3: %.2gs" % (time.perf_counter() - tstart))
        return plot

    def makeplot5(self, plot=None, xdata=None, ydata=None, zdata=None, config=None):
        """
        Plot mass vs charge grid
        :param plot: Plot object to use
        :param xdata: x data (2D array) default is eng.data.massdat
        :param ydata: y data (1D array) default is eng.data.ztab
        :param zdata: z data (1D array) default is eng.data.massgrid
        :param config: Config object to use
        :return: plot object
        """
        if plot is None:
            plot = plot2d.Plot2dBase()
        if xdata is None:
            xdata = self.data.massdat
        if ydata is None:
            ydata = self.data.ztab
        if zdata is None:
            zdata = self.data.massgrid
        if config is None:
            config = self.config

        if self.config.batchflag == 0:
            tstart = time.perf_counter()
            plot.contourplot(xvals=xdata[:, 0], yvals=ydata, zgrid=zdata, config=config, title="Mass vs. Charge",
                             test_kda=True)
            print("Plot 5: %.2gs" % (time.perf_counter() - tstart))
        return plot

    def makeplot2(self, plot=None, data=None, pks=None, config=None):
        """
        Plot mass data and peaks if possible in self.view.plot2
        :param plot: Plot object. Default is None, which will set plot to creating a new plot
        :param data: 2D data with mass in the first column and intensity in the second. Default is self.data.massdat
        :param pks: Peaks object. Default is self.pks
        :param config: Config objet. Default is self.config
        :return plot: The Plot object
        """
        if plot is None:
            plot = plot1d.Plot1dBase()
        if data is None:
            data = self.data.massdat
        if pks is None:
            pks = self.pks
        if config is None:
            config = self.config
        if config.batchflag == 0 and data.shape[1] == 2 and len(data) >= 2:
            tstart = time.perf_counter()
            plot.plotrefreshtop(data[:, 0], data[:, 1],
                                "Zero-charge Mass Spectrum", "Mass (Da)",
                                "Intensity", "Mass Distribution", config, test_kda=True,
                                nopaint=True)
            if pks.plen > 0:
                for p in pks.peaks:
                    if p.ignore == 0:
                        plot.plotadddot(p.mass, p.height, p.color, p.marker)
            plot.repaint()
            tend = time.perf_counter()
            print("Plot 2: %.2gs" % (tend - tstart))
        if data.shape[1] != 2 or len(data) < 2:
            print("Data Too Small. Adjust parameters.", data)
        return plot

    def makeplot4(self, plot=None, data=None, pks=None, config=None):
        """
        Plots isolated peaks against the data in self.view.plot4.
        Will plot dots at peak positions.
        If possible, will plot full isolated spectra.
        :param plot: Plot object. Default is None, which will set plot to creating a new plot
        :param data: 2D data with mass in the first column and intensity in the second. Default is self.data.massdat
        :param pks: Peaks object. Default is self.pks
        :param config: Config objet. Default is self.config
        :return plot: The Plot object
        """
        if plot is None:
            plot = plot1d.Plot1dBase()
        if data is None:
            data = self.data.data2
        if pks is None:
            pks = self.pks
        if config is None:
            config = self.config
        if config.batchflag == 0:
            tstart = time.perf_counter()
            # This plots the normal 1D mass spectrum
            plot.plotrefreshtop(data[:, 0], data[:, 1],
                                "Data with Offset Isolated Species", "m/z (Th)",
                                "Normalized and Offset Intensity", "Data", config, nopaint=True)
            num = 0
            # Corrections for if Isotope mode is on
            if config.isotopemode == 1:
                try:
                    stickmax = np.amax(np.array([p.stickdat for p in pks.peaks]))
                except (AttributeError, ValueError):
                    stickmax = 1.0
            else:
                stickmax = 1.0
            # Loop through each peak
            for i, p in enumerate(pks.peaks):
                # Check if the peak is ignored
                if p.ignore == 0:
                    # Check if the mztabs are empty
                    if (not ud.isempty(p.mztab)) and (not ud.isempty(p.mztab2)):
                        mztab = np.array(p.mztab)
                        mztab2 = np.array(p.mztab2)
                        maxval = np.amax(mztab[:, 1])
                        # Filter all peaks where the deconvolved intensity is above the relative threshold
                        b1 = mztab[:, 1] > config.peakplotthresh * maxval
                        # Plot the filtered peaks as dots on the spectrum
                        plot.plotadddot(mztab2[b1, 0], mztab2[b1, 1], p.color, p.marker)
                    # Check if convolved data is present
                    if not ud.isempty(p.stickdat):
                        # Plot the offset reconvolved data from the isolated species
                        plot.plotadd(data[:, 0], np.array(p.stickdat) / stickmax - (
                                num + 1) * config.separation, p.color, "useless label")
                    num += 1
            plot.repaint()
            tend = time.perf_counter()
            print("Plot 4: %.2gs" % (tend - tstart))
            return plot

    def gen_html_report(self, event=None, plots=None, interactive=False):
        """
        Generate an HTML report of the current UniDec run.
        :param event: Unused Event
        :param plots: List of plots to include in the report. Must be 2D with row and column format.
        :param interactive: If True, will include interactive plots. Default is False.
        :return: None
        """
        outfile = self.config.outfname + "_report.html"
        html_open(outfile)
        html_title(self.config.filename, outfile)

        peaks_df = self.pks.to_df()
        colors = [p.color for p in self.pks.peaks]

        if len(peaks_df) > 0:
            df_to_html(peaks_df, outfile, colors=colors)

        # array_to_html(np.transpose(self.matchlist), outfile,
        #              cols=["Measured Mass", "Theoretical Mass", "Error", "Match Name"])

        if plots is None:
            plot = self.makeplot2()
            plot2 = self.makeplot4()
            # plot5 = self.makeplot5()
            # plot3 = self.makeplot3()
            # plot6 = self.makeplot6()
            plots = [[plot, plot2]]  # , [plot5, plot3], [plot6, None]]
        if len(np.shape(plots)) != 2 and len(plots) > 0:
            try:
                # Reshape 1D array to 2D with 2 columns
                plots = np.reshape(plots, (int(len(plots) / 2), 2))
            except Exception:
                plots = np.reshape(plots[:-1], (int(len(plots[:-1]) / 2), 2))
                lastrow = np.array([plots[-1], None])
                plots = np.vstack((plots, lastrow))

        svg_grid = []
        figure_list = []
        for row in plots:
            svg_row = []
            goodrow = False
            for c in row:
                if c is not None and c.flag:
                    goodrow = True
                    if c.is2d:
                        png_str = c.get_png()
                        png_html = png_to_html(png_str)
                        svg_row.append(png_html)
                    else:
                        svg_row.append(c.get_svg())
                    figure_list.append(c.figure)
                else:
                    svg_row.append("<p></p>")
            if goodrow:
                svg_grid.append(svg_row)

        svg_grid_string = wrap_to_grid(svg_grid, outfile)

        if interactive:
            for f in figure_list:
                try:
                    fig_to_html_plotly(f, outfile)
                except Exception:
                    pass

        try:
            spectra_df = self.data.attrs_to_df()
            if len(spectra_df) > 0:
                spectra_df.drop(["beta", "error", "iterations", "psig", "rsquared", "zsig", "mzsig", "length_mass",
                                 "length_mz", "time"], axis=1, inplace=True)
                colors2 = self.data.get_colors()
                df_to_html(spectra_df, outfile, colors=colors2)
        except Exception:
            pass

        config_dict = self.config.get_config_dict()
        config_htmlstring = dict_to_html(config_dict)
        to_html_collapsible(config_htmlstring, title="UniDec Parameters", outfile=outfile, htmltext=True)

        html_close(outfile)

        os.system(self.config.opencommand + "\"" + outfile + "\"")
