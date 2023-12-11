import numpy as np
import matplotlib.pyplot as plt
import os
import scipy.optimize as opt
# import matplotlib.cm as cm
import matplotlib as mpl
from scipy.stats import f
import networkx as nx
import re
import warnings
from copy import deepcopy
import sys
import time
import scipy.special as special
import multiprocessing

__author__ = 'Michael.Marty'

warnings.filterwarnings('ignore', 'The iteration is not making good progress')


# TODO: Clean up this code...
# TODO: Allow simple monomer/dimer with no ligands


def safedivide(a, b):
    c = deepcopy(a)
    if b != 0:
        c = a / b
    # c[b!=0]=a[b!=0]/b[b!=0]
    return c


def fixlists(list1, list2):
    unique = np.sort(np.unique(np.array(list1)))
    index = np.arange(0, len(unique))
    d = {}
    for i in index:
        d.update({unique[i]: i})
    for i in range(0, len(list1)):
        list1[i] = d[list1[i]]
    for i in range(0, len(list2)):
        for j in range(0, len(list2[i])):
            list2[i][j] = d[list2[i][j]]
    return list1, list2


def sortarray(a):
    shape = a.shape
    a = [list(a[i]) for i in range(0, shape[0])]
    return np.array(sorted(a, key=lambda row: row[0]))


def CompareModels(e1, e2, p1, p2, n, thresh):
    plvl = 0
    if e1 > e2 and p1 >= p2:
        plvl = 1
        model = 1
    if e1 < e2 and p1 <= p2:
        plvl = 1
        model = 0
    if e1 > e2 and p1 < p2:
        v1 = p2 - p1
        v2 = n - p2
        fval = (e1 - e2) / (v1 * (e2 / v2))
        plvl = f.cdf(fval, v1, v2)
        if plvl >= thresh:
            model = 1
        else:
            model = 0
    if e2 > e1 and p2 < p1:
        v1 = p1 - p2
        v2 = n - p1
        fval = (e2 - e1) / (v1 * (e1 / v2))
        plvl = f.cdf(fval, v1, v2)
        if plvl >= thresh:
            model = 0
        else:
            model = 1
    if e1 == e2 and p1 > p2:
        plvl = 1
        model = 1
    if e1 == e2 and p1 < p2:
        plvl = 1
        model = 0
    if e1 == e2 and p1 == p2:
        plvl = 1
        model = 0
        print("Models identical")
    if p1 > n or p2 > n:
        plvl = 1
        print("Model Too Big")
    return plvl, model


def Ftest(e1, e2, p1, p2, n):
    dof1 = n - p1
    dof2 = n - p2
    plvl = 0
    if e1 > e2 and p1 == p2:
        v1 = e2 / float(dof2)
        v2 = e1 / float(dof1)
        fval = (v2 - v1) / v1
        plvl = f.cdf(fval, p2 - 1, dof2)
    if e2 > e1 and p2 == p1:
        v1 = e1 / float(dof1)
        v2 = e2 / float(dof2)
        fval = (v2 - v1) / v1
        plvl = f.cdf(fval, p1 - 1, dof1)
    if e1 == e2 and p1 == p2:
        plvl = 0
    if p1 != p2:
        print("Error p1!=p2")
        exit()
    return plvl


def findend(reactions, end):
    for i in range(0, len(reactions)):
        if np.all(reactions[i, 2] == end):
            return i
    return -1


def findpaths(reactions):
    paths = []
    for i in range(0, len(reactions)):
        path = []
        flag = 0
        end = reactions[i, 2]
        while flag == 0:
            react = findend(reactions, end)
            if react != -1:
                path.append(react)
                end = reactions[react, 0]
            else:
                flag = 1
        paths.append(path[::-1])
    return paths


def GetError(data, pfree, lfree, ureact, prottab, ligtab, paths, kds, weights, nodelist, nfactors):
    intgrid = MakeGrid(pfree, lfree, ureact, prottab, ligtab, paths, kds, nfactors)[0]
    extract = np.array([intgrid[row[0], row[1]] for row in nodelist])
    if np.any(extract < 0.) or np.any(kds < 0.):
        # extract=extract*0
        return sys.maxsize * weights
    summed = np.sum(extract)
    if summed > 0:
        return (weights * (extract / summed - data)) ** 2
    else:
        return (weights * (extract - data)) ** 2


def MinFree(ptot, ltot, ureact, prottab, ligtab, paths, kds, pguess, lguess, nfactors):
    pfree = pguess
    lfree = lguess
    p0 = np.array([pfree, lfree])
    # fit=opt.leastsq(GetFree,p0,args=(ptot,ltot,ureact,prottab,ligtab,paths,kds))
    fit = opt.fsolve(GetFree, p0, args=(ptot, ltot, ureact, prottab, ligtab, paths, kds, nfactors))
    pfree = np.clip(fit[0], 0, ptot)
    lfree = np.clip(fit[1], 0, ltot)
    return pfree, lfree


def GetFree(pO, ptot, ltot, ureact, prottab, ligtab, paths, kds, nfactors):
    pfree = pO[0]
    lfree = pO[1]
    intgrid, sumprot, sumlig = MakeGrid(pfree, lfree, ureact, prottab, ligtab, paths, kds, nfactors)
    return [(ptot - sumprot) ** 1, (ltot - sumlig) ** 1]


def MakeGrid(pfree, lfree, ureact, prottab, ligtab, paths, kds, nfactors):
    intgrid = np.zeros_like(prottab)
    intgrid[1, 0] = pfree
    try:
        intgrid[0, 1] = lfree
    except:
        pass

    h = 1
    for i in range(0, len(ureact)):
        if len(ureact[i, 2]) > 2:
            nump = ureact[i, 2, 0]
            numl = ureact[i, 2, 1]
            numa = ureact[i, 2, 2]
            if numa == 2 and nump == 0 and numl == 0:
                # print(kds[i])
                h = kds[i]
    for i in range(0, len(ureact)):
        nump = ureact[i, 2, 0]
        numl = ureact[i, 2, 1]
        denom = np.prod([kds[j] for j in paths[i]])
        if denom != 0:
            intgrid[nump, numl] += (pfree ** nump) * (lfree ** (numl * h)) / denom
    if nfactors is not None:
        intgrid = nfactors * intgrid
    sumprot = np.sum(prottab * intgrid)
    sumlig = np.sum(ligtab * intgrid)
    return intgrid, sumprot, sumlig


def draw_graph_structure(graph1, graph2, kdmap, ax=None):
    # extract nodes from graph
    # print(graph1)
    nodes = set([n1 for n1, n2 in graph1] + [n2 for n1, n2 in graph1])
    # print(nodes)
    # create networkx graph
    G = nx.Graph()
    G2 = nx.Graph()
    p = {}
    r = re.compile("([a-zA-Z]+)([0-9]+)")
    # add nodes
    for node in nodes:
        if "*" in node:
            split1 = node.split("*")
            split2 = split1[0].split("L")
            split3 = split2[0].split("P")
            z = int(split3[1])  # P
            x = int(split2[1])  # L
            y = int(split1[1])  # *
            if z < 1:
                if x > 1:
                    coord = [-0.5, x]
                else:
                    coord = [-0.5, y]
            else:
                if x < 1:
                    coord = [0, 1]
                else:
                    coord = [x, y]
        else:
            split1 = node.split("L")
            split2 = split1[0].split("P")
            coord = [int(split1[1]), int(split2[1])]

        # print(node)

        G.add_node(node)
        G2.add_node(node)
        p.update({node: coord})

    # add edges
    e = {}
    for edge in graph1:
        G.add_edge(edge[0], edge[1])
    for edge in graph2:
        G2.add_edge(edge[0], edge[1])
    # edge labels
    for i in range(0, len(graph1)):
        edge = graph1[i]
        e.update({(edge[0], edge[1]): "KD%d" % kdmap[i]})

    bbox = dict(boxstyle='round', ec=(1.0, 1.0, 1.0), fc=(0.0, 0.0, 0.0), )
    # draw graph
    if ax is None:
        plt.figure(figsize=(10, 10))
    pos = nx.shell_layout(G)
    nx.draw_networkx_nodes(G, p, node_size=1500, ax=ax)
    nx.draw_networkx_edges(G, p, width=15, ax=ax)
    nx.draw_networkx_edges(G2, p, width=10, edge_color="b", ax=ax)
    nx.draw_networkx_labels(G, p, ax=ax)
    nx.draw_networkx_edge_labels(G, p, edge_labels=e, rotate=False, bbox=bbox, font_color="w", ax=ax)
    # limits=plt.axis('off')
    # show graph
    if ax is None:
        plt.show()
    else:
        pass
        # ax.show()


def draw_graph(graph1, graph2, kds, errors, kdmap, header, ax=None, ulabel=""):
    # extract nodes from graph
    nodes = set([n1 for n1, n2 in graph1] + [n2 for n1, n2 in graph1])

    # create networkx graph
    G = nx.Graph()
    G2 = nx.Graph()
    p = {}
    r = re.compile("([a-zA-Z]+)([0-9]+)")
    # add nodes
    for node in nodes:
        if "*" in node:
            split1 = node.split("*")
            split2 = split1[0].split("L")
            split3 = split2[0].split("P")
            z = int(split3[1])  # P
            x = int(split2[1])  # L
            y = int(split1[1])  # *
            if z < 1:
                if x > 1:
                    coord = [-0.5, x]
                else:
                    coord = [-0.5, y]
            else:
                if x < 1:
                    coord = [0, 1]
                else:
                    coord = [x, y]
        else:
            split1 = node.split("L")
            split2 = split1[0].split("P")
            coord = [int(split1[1]), int(split2[1])]

        G.add_node(node)
        G2.add_node(node)
        p.update({node: coord})

    # add edges
    e = {}
    for edge in graph1:
        G.add_edge(edge[0], edge[1])

    for i in range(0, len(graph1)):
        edge = graph1[i]
        string1 = "KD%d=%.3g" % (kdmap[i], kds[kdmap[i]])
        if errors is not None:
            string2 = "\n\u00B1 %.3g" % errors[kdmap[i]]
            stringfin = string1 + string2
        else:
            stringfin = string1

        if stringfin != "":
            stringfin += " " + ulabel
        e.update({(edge[0], edge[1]): stringfin})

    for edge in graph2:
        G2.add_edge(edge[0], edge[1])
    bbox = dict(boxstyle='round', ec=(1.0, 1.0, 1.0), fc=(0.0, 0.0, 0.0), )
    # draw graph
    pos = nx.shell_layout(G)
    nx.draw_networkx_nodes(G, p, node_size=1500, ax=ax)
    nx.draw_networkx_edges(G, p, width=15, ax=ax)
    nx.draw_networkx_edges(G2, p, width=10, edge_color="b", ax=ax)
    nx.draw_networkx_labels(G, p, ax=ax)
    nx.draw_networkx_edge_labels(G, p, edge_labels=e, rotate=False, bbox=bbox, font_color="w", ax=ax)
    # limits=plt.axis('off')
    if ax is None:
        plt.savefig(header + "_graph.png")
        # show graph
        # plt.show()


def MinFreeError(kds, data, weights, kdargs):
    out = [MinFree(kdargs.pconc[i], kdargs.lconc[i], kdargs.ureact, kdargs.nprottab, kdargs.nligtab, kdargs.paths, kds,
                   kdargs.pfrees[i], kdargs.lfrees[i], kdargs.nfactors) for i in range(0, len(kdargs.lconc))]
    out = np.array(out)
    kdargs.pfrees = out[:, 0]
    kdargs.lfrees = out[:, 1]
    if np.any(kdargs.pfrees == 0) or np.any(kdargs.lfrees[1:] == 0):
        out = [
            MinFree(kdargs.pconc[i], kdargs.lconc[i], kdargs.ureact, kdargs.nprottab, kdargs.nligtab, kdargs.paths, kds,
                    kdargs.pfrees[i], kdargs.lfrees[i], kdargs.nfactors) for i in range(0, len(kdargs.lconc))]
        out = np.array(out)
        kdargs.pfrees = out[:, 0]
        kdargs.lfrees = out[:, 1]
    errors = np.ravel([GetError(data[:, i], kdargs.pfrees[i], kdargs.lfrees[i], kdargs.ureact, kdargs.nprottab,
                                kdargs.nligtab, kdargs.paths, kds, weights[:, i], kdargs.nodelist, kdargs.nfactors) for
                       i in range(0, len(kdargs.lconc))])
    return errors


def GetDegenKD(kds, kdargs):
    if len(kdargs.degen > 0):
        dktot = []
        for i in range(0, len(kdargs.lconc)):
            out = [
                MinFree(kdargs.pconc[i], kdargs.lconc[i], kdargs.ureact, kdargs.nprottab, kdargs.nligtab, kdargs.paths,
                        kds, kdargs.pfrees[i], kdargs.lfrees[i], kdargs.nfactors) for i in range(0, len(kdargs.lconc))]
            out = np.array(out)
            pfrees = out[:, 0]
            lfrees = out[:, 1]
            degenkd = []
            grid = MakeGrid(pfrees[i], lfrees[i], kdargs.ureact, kdargs.nprottab, kdargs.nligtab, kdargs.paths, kds,
                            kdargs.nfactors)[0]
            for react in kdargs.degen[:, 1]:
                if grid[react[2, 0], react[2, 1]] != 0:
                    val = grid[react[0, 0], react[0, 1]] * grid[react[1, 0], react[1, 1]] / grid[
                        react[2, 0], react[2, 1]]
                else:
                    val = 0
                degenkd.append(val)
            dktot.append(degenkd)
        dktot = np.array(dktot)
        degenkd = np.mean(dktot, axis=0)
    else:
        degenkd = []
    return degenkd


def Minimize(data, kdargs, **kwargs):
    if "weights" in kwargs:
        weights = kwargs["weights"]
    else:
        weights = kdargs.weights
    fit = opt.leastsq(MinFreeError, kdargs.kds, args=(data, weights, kdargs))
    degenfit = GetDegenKD(fit[0], kdargs)
    return np.array([fit[0], degenfit])


def BootMin(data, kdargs, weights):
    return np.concatenate(Minimize(data, kdargs, weights=weights))


def BootMinWorker(queue, results_queue):
    while True:
        try:
            out = queue.get(False)
        except Exception as e:
            break
        data = out[0]
        kdargs = out[1]
        weights = out[2]
        results = np.concatenate(Minimize(data, kdargs, weights=weights))
        results_queue.put(results)
    return


def make_graph_agg(reactions):
    graph = [("P" + str(int(reactions[i, 0, 0])) + "L" + str(int(reactions[i, 0, 1])) + "*" + str(
        int(reactions[i, 0, 2])),
              "P" + str(int(reactions[i, 2, 0])) + "L" + str(int(reactions[i, 2, 1])) + "*" + str(
                  int(reactions[i, 2, 2]))) for i in
             range(0, len(reactions))]
    return graph


class kdstruct:
    def __init__(self):
        self.pconc = []
        self.lconc = []
        self.ureact = []
        self.nprottab = []
        self.nligtab = []
        self.paths = []
        self.pfrees = []
        self.lfrees = []
        self.weights = []
        self.nodelist = []
        self.kds = []
        self.degen = []
        self.maxsites = None
        self.nfactors = []


class KDmodel:
    def __init__(self, data, pconc, lconc, nodelist=None, header=None, numtotprot=0, numtotlig=0, removeoutliers=False,
                 plot1=None, plot2=None, plot3=None, bootnum=1, maxsites=0, maxligagg=1, hill=False, label="",
                 cmap='rainbow', **kwargs):
        self.outlierflag = removeoutliers
        self.cmap = cmap
        self.label = label
        self.plot1 = plot1
        self.plot2 = plot2
        self.plot3 = plot3
        self.protflag = 0
        self.ligflag = 0
        self.mode = 0
        self.header = header
        if self.header is None:
            self.header = os.getcwd()
        self.maxsites = maxsites
        self.kdargs = kdstruct()
        self.randfit = None
        self.bootnum = bootnum
        self.ligaggmode = False
        self.hill = hill
        if self.hill:
            self.ligaggmode = True

        # Setting up the model
        if numtotprot == 0:
            numtotprot = 1
            numtotlig = len(data) - 1
        self.numtotprot = numtotprot
        self.numtotlig = numtotlig
        self.nprot = np.arange(0, self.numtotprot + 1, dtype=float)
        self.nlig = np.arange(0, self.numtotlig + 1, dtype=float)
        self.kdargs.nprottab, self.kdargs.nligtab = np.meshgrid(self.nprot, self.nlig, indexing="ij")

        if nodelist is None:
            nprot = np.ones(numtotlig + 1)
            nodelist = np.transpose([nprot, self.nlig]).astype(int)
        self.kdargs.nodelist = np.array(nodelist)

        # Getting in the experimenatl data
        self.data = np.array(data)
        self.kdargs.pconc = np.array(pconc)
        self.kdargs.lconc = np.array(lconc)
        self.kdargs.weights = np.ones_like(self.data)
        # Normalize each line to a sum of 1
        for i in range(0, len(self.data[0])):
            self.data[:, i] = self.data[:, i] / np.sum(self.data[:, i])

        # initial guess for pfree and lfree
        self.kdargs.pfrees = pconc
        self.kdargs.lfrees = lconc

        self.fixedligmodel = None
        self.fixedprotmodel = None
        self.fixedaggmodel = None
        self.findmodelflag = False
        if "prot" in kwargs:
            if "agg" in kwargs:
                print("Note: The agg or prot model keywords could be undefined depending on the mode")
            key = kwargs["prot"]
            # Key Should be free, parallel, one
            if key == "test":
                self.fixedprotmodel = None
                self.findmodelflag = True
            else:
                self.fixedprotmodel = key
        if "lig" in kwargs:
            key = kwargs["lig"]
            # Key Should be free, parallel, one
            if key == "test":
                self.fixedligmodel = None
                self.findmodelflag = True
            else:
                self.fixedligmodel = key
        if "agg" in kwargs:
            key = kwargs["agg"]
            # Key Should be free, parallel, one
            if key == "test":
                self.fixedaggmodel = None
                self.findmodelflag = True
            else:
                self.fixedaggmodel = key

        self.maxligagg = maxligagg
        if self.maxligagg == -1:
            self.maxligagg = self.numtotlig
        if self.maxligagg > 1:
            self.ligaggmode = True
        if self.ligaggmode and self.numtotprot > 1:
            print("ERROR: Simultaneous free protein and free ligand aggregation are not supported")
            print("\tSet either maxligagg or numtotprot to 1.")
            raise Exception

        if self.ligaggmode:
            print("Ligand Aggregation Mode")
            if "preview" in kwargs:
                preview = kwargs["preview"]
                if preview:
                    self.SetupModelAgg(**kwargs)
                    sys.exit()
            self.plotflag = 1
            self.RunKDFit(**kwargs)

        else:
            if "preview" in kwargs:
                preview = kwargs["preview"]
                if preview:
                    self.SetupModel(**kwargs)
                    sys.exit()

            if self.findmodelflag:
                self.plotflag = 0
                self.protmodel, self.ligmodel = self.FindBestModel(self.fixedprotmodel, self.fixedligmodel)
                self.plotflag = 1
                self.RunKDFit(prot=self.protmodel, lig=self.ligmodel)
            else:
                self.plotflag = 1
                self.RunKDFit(**kwargs)

    def FindBestModel(self, fixedprotmodel, fixedligmodel):
        print("\n\nAutomatic Determination of Best Model\n\n")
        tests = ["free", "parallel", "series", "one"]
        if self.numtotprot == 1:
            ptest = ["free"]
            ltest = ["free", "one"]
        else:
            ptest = tests
            ltest = tests
        if self.numtotlig < 1:
            ltest = ["free"]
            ptest = ["free", "one"]

        if fixedprotmodel is not None:
            ptest = [fixedprotmodel]
        if fixedligmodel is not None:
            ltest = [fixedligmodel]

        tests = [ptest, ltest]
        testsflat = []
        [[testsflat.append([i, j]) for i in tests[0]] for j in tests[1]]
        testsflat = np.array(testsflat)

        fits = np.array([self.RunKDFit(prot=test[0], lig=test[1]) for test in testsflat])
        print("SSE and #Parameters")
        print(fits)

        sortindex = np.argsort(fits, axis=0)[::-1, 0]
        fits = np.concatenate([fits.transpose(), np.transpose(testsflat)]).transpose()
        fitssorted = np.empty_like(fits)
        for i in range(0, len(fits)):
            fitssorted[i] = fits[sortindex[i]]
        compares = np.array([CompareModels(float(fitssorted[i][0]), float(fitssorted[i + 1][0]),
                                           float(fitssorted[i][1]), float(fitssorted[i + 1][1]),
                                           len(np.ravel(self.data)), 0.99) for i in range(0, len(fits) - 1)])
        test = 0
        print("P Values: ")
        print(compares)
        while test < len(compares) and compares[test][1] == 1:
            test += 1
        bestmodel = fitssorted[test][2:]
        print("\nBest Model is:\n\tProtein Model: ", bestmodel[0], "\n\tLigand Model: ", bestmodel[1])
        return bestmodel

    def RunKDFit(self, **kwargs):

        if "plot" in kwargs:
            self.plotflag = kwargs["plot"]

        # Setup Model
        if self.ligaggmode:
            self.SetupModelAgg(**kwargs)
            # Initial guess for KD
            self.kdargs.kds = np.array([np.mean(self.kdargs.lconc) for i in range(0, self.numkd)])
            if self.hill:
                self.kdargs.kds[0] = 1
        else:
            self.SetupModel(**kwargs)
            # Initial guess for KD
            self.kdargs.kds = np.array([np.mean(self.kdargs.lconc) for i in range(0, self.numkd)])

        if np.mean(self.kdargs.lconc) == 0:
            self.kdargs.kds = np.array([np.mean(self.kdargs.pconc)/2. for i in range(0, self.numkd)])
        #print(self.kdargs.kds)

        if np.any(self.kdargs.kds == 0):
            print("Error in guessing KDs", self.kdargs.kds)
            self.kdargs.kds = np.array([np.mean(self.kdargs.pconc) for i in range(0, self.numkd)])
            if np.any(self.kdargs.kds == 0):
                print("Error in guessing KDs again", self.kdargs.kds)
                self.kdargs.kds = np.array([1 for i in range(0, self.numkd)])

        if "model" in kwargs:
            if kwargs["model"] == 1:
                self.stddevs = None  # np.zeros_like(np.concatenate(self.fit))
                self.fit = np.array([self.kdargs.kds, GetDegenKD(self.kdargs.kds, self.kdargs)])
                self.MakeFitGrid()
                self.GraphPlot(ax=self.plot2)

                plt.figure()
                ax = plt.subplot(121)
                self.PlotTrace(topax=ax)

                self.fit = np.array([self.kdargs.kds, GetDegenKD(self.kdargs.kds, self.kdargs)])
                self.MakeFitGrid()
                ax = plt.subplot(122)
                self.PlotTrace(topax=ax)
                mng = plt.get_current_fig_manager()
                mng.window.state('zoomed')
                plt.show()
                exit()

        # Minimization
        self.fit = Minimize(self.data, self.kdargs)
        self.kdargs.kds = self.fit[0]
        self.degenkds = self.fit[1]

        # Results
        self.MakeFitGrid()

        self.residuals = np.ravel((self.fitgrid - self.data) ** 1)
        self.stdres = np.std(self.residuals)
        # print "Std Dev of Residuals: ",self.stdres
        self.error = np.sum((self.fitgrid - self.data) ** 2)
        print("Sum of Squared Error: ", self.error)
        # print "Number of Parameters: ",self.numkd
        # print "Number of Data Points: ",self.data.size
        print("Degrees of Freedom: ", self.data.size - self.numkd)

        if self.plotflag == 1:
            self.RunBootstrap()
            self.GraphPlot(ax=self.plot2)
            self.PlotTrace(topax=self.plot1)
            self.PlotHist(topax=self.plot3)
        return self.Return()

    def Return(self):
        return [self.error, self.numkd]

    def SetupModel(self, **kwargs):
        # Define path favoring ligand addition->0 or protein addition->1
        if "prot" in kwargs:
            key = kwargs["prot"]
            # Key Should be free, parallel, one
            if key == "free":
                self.mode = 0
                self.protflag = 0
            elif key == "parallel":
                self.mode = 1
                self.protflag = 3
            elif key == "one":
                self.mode = 1
                self.protflag = 2
            elif key == "series":
                self.protflag = 1
                self.mode = 1

        if "lig" in kwargs:
            key = kwargs["lig"]
            # Key Should be free, parallel, one
            if key == "free":
                self.ligflag = 0
            elif key == "parallel":
                self.ligflag = 3
            elif key == "one":
                self.mode = 0
                self.ligflag = 2
            elif key == "series":
                self.ligflag = 1

        # Defining the potential nodes and reactions
        # Note: the nodes are species that are experimentally measured.
        self.reactions = []
        # self.nodes = []
        for i in range(1, self.numtotprot + 1):
            for j in self.nlig:
                if j + 1 <= self.numtotlig and [i, j] != [0, 0] and (
                        (j + 1 <= i * self.maxsites and j <= i * self.maxsites) or self.maxsites == 0):
                    # self.nodes.append(tuple([int(i), int(j)]))
                    # self.nodes.append(tuple([int(i), int(j + 1)]))
                    self.reactions.append([[int(i), int(j)], [0, 1], [int(i), int(j + 1)]])
                if i + 1 <= self.numtotprot and [i, j] != [0, 0] and (j <= i * self.maxsites or self.maxsites == 0):
                    # self.nodes.append(tuple([int(i), int(j)]))
                    # self.nodes.append(tuple([int(i + 1), int(j)]))
                    self.reactions.append([[int(i), int(j)], [1, 0], [int(i + 1), int(j)]])
        # self.nodes = np.vstack([tuple(row) for row in self.nodes])
        # self.nodes = sortarray(self.nodes)
        self.reactions = np.array(self.reactions)

        # Fishing out degenerate KD values
        self.kdargs.degen = []
        for i in self.reactions:
            for j in self.reactions:
                if np.all(i[[2]] == j[[2]]):
                    if self.mode == 0:
                        if np.all(i[[1]] == [0, 1]) and np.all(j[[1]] == [1, 0]):
                            self.kdargs.degen.append([i, j])
                    else:
                        if np.all(i[[1]] == [1, 0]) and np.all(j[[1]] == [0, 1]):
                            self.kdargs.degen.append([i, j])
        self.kdargs.degen = np.array(self.kdargs.degen)

        # Deleting degenate KD values to get number of unique KDs
        self.kdargs.ureact = []
        for i in self.reactions:
            if len(self.kdargs.degen) != 0:
                if not np.any([np.all(i == j) for j in self.kdargs.degen[:, 1]]):
                    self.kdargs.ureact.append(i)
            else:
                self.kdargs.ureact.append(i)
        self.kdargs.ureact = np.array(self.kdargs.ureact)
        self.kdmap = list(range(0, len(self.kdargs.ureact)))

        # Finding the path to each KD. Will be used for calculating the grid.
        self.kdargs.paths = findpaths(self.kdargs.ureact)
        self.ModPaths()
        self.kdmap, self.kdargs.paths = fixlists(self.kdmap, self.kdargs.paths)
        self.numkd = len(np.unique(np.concatenate(self.kdargs.paths)))
        if len(self.kdargs.degen > 0):
            self.reactions = np.concatenate([self.kdargs.ureact, self.kdargs.degen[:, 1]])
            self.kdmap = np.concatenate(
                [np.array(self.kdmap), np.arange(self.numkd, self.numkd + len(self.kdargs.degen[:, 1]))])

        # Set up graph data for plot
        self.graph = [("P" + str(int(self.reactions[i, 0, 0])) + "L" + str(int(self.reactions[i, 0, 1])),
                       "P" + str(int(self.reactions[i, 2, 0])) + "L" + str(int(self.reactions[i, 2, 1]))) for i in
                      range(0, len(self.reactions))]
        self.ugraph = [("P" + str(int(self.kdargs.ureact[i, 0, 0])) + "L" + str(int(self.kdargs.ureact[i, 0, 1])),
                        "P" + str(int(self.kdargs.ureact[i, 2, 0])) + "L" + str(int(self.kdargs.ureact[i, 2, 1]))) for i
                       in range(0, len(self.kdargs.ureact))]
        self.nodenames = ["P" + str(int(self.kdargs.nodelist[i, 0])) + "L" + str(int(self.kdargs.nodelist[i, 1])) for i
                          in
                          range(0, len(self.kdargs.nodelist))]

        # Calculate Statistical Factors
        if self.maxsites != 0 and self.maxsites is not None:
            self.kdargs.maxsites = self.maxsites
            narr = np.array(self.kdargs.nprottab * self.maxsites)
            iarr = self.kdargs.nligtab
            niarr = np.clip(narr - iarr, 0, sys.maxsize)
            self.kdargs.nfactors = special.factorial(narr) / (
                    special.factorial(niarr) * special.factorial(iarr))
        else:
            self.kdargs.nfactors = None
        # Quick Plot to test structure. Note: exits the program without fitting.
        if "preview" in kwargs:
            preview = kwargs["preview"]
            if preview:
                draw_graph_structure(self.graph, self.ugraph, self.kdmap, ax=self.plot2)

    def ModPaths(self):
        # All ligand binding reactions have the same KD
        if self.ligflag == 2:
            for i in range(0, len(self.kdargs.paths)):
                if (len(self.kdargs.paths[i]) == 1 and np.all(
                        self.kdargs.ureact[self.kdargs.paths[i][0], 1] == [0, 1])):
                    onepath = self.kdargs.paths[i][0]
            for i in range(0, len(self.kdargs.paths)):
                for j in range(0, len(self.kdargs.paths[i])):
                    if np.all(self.kdargs.ureact[self.kdargs.paths[i][j], 1] == [0, 1]):
                        self.kdmap[self.kdargs.paths[i][j]] = onepath
                        self.kdargs.paths[i][j] = onepath
        # All protein binding reactions have the same KD
        if self.protflag == 2:
            for i in range(0, len(self.kdargs.paths)):
                if (len(self.kdargs.paths[i]) == 1 and np.all(
                        self.kdargs.ureact[self.kdargs.paths[i][0], 1] == [1, 0])):
                    onepath = self.kdargs.paths[i][0]
            for i in range(0, len(self.kdargs.paths)):
                for j in range(0, len(self.kdargs.paths[i])):
                    if np.all(self.kdargs.ureact[self.kdargs.paths[i][j], 1] == [1, 0]):
                        self.kdmap[self.kdargs.paths[i][j]] = onepath
                        self.kdargs.paths[i][j] = onepath
        # For a given ligand number, all protein binding reactions have the same KD
        if self.protflag == 1:
            ulig = []
            ufirst = []
            for i in range(0, len(self.kdargs.paths)):
                for j in range(0, len(self.kdargs.paths[i])):
                    if np.all(self.kdargs.ureact[self.kdargs.paths[i][j], 1] == [1, 0]):
                        testlig = self.kdargs.ureact[self.kdargs.paths[i][j], 0][1]
                        if testlig not in ulig:
                            ulig.append(testlig)
                            ufirst.append(self.kdargs.paths[i][j])
            for i in range(0, len(self.kdargs.paths)):
                for j in range(0, len(self.kdargs.paths[i])):
                    if np.all(self.kdargs.ureact[self.kdargs.paths[i][j], 1] == [1, 0]):
                        testlig = self.kdargs.ureact[self.kdargs.paths[i][j], 0][1]
                        onepath = ufirst[ulig.index(testlig)]
                        self.kdmap[self.kdargs.paths[i][j]] = onepath
                        self.kdargs.paths[i][j] = onepath
        # For a given protein association reaction, the KD is the same across all ligand bound states
        if self.protflag == 3:
            ulig = []
            ufirst = []
            for i in range(0, len(self.kdargs.paths)):
                for j in range(0, len(self.kdargs.paths[i])):
                    if np.all(self.kdargs.ureact[self.kdargs.paths[i][j], 1] == [1, 0]):
                        testlig = self.kdargs.ureact[self.kdargs.paths[i][j], 0][0]
                        if testlig not in ulig:
                            ulig.append(testlig)
                            ufirst.append(self.kdargs.paths[i][j])
            for i in range(0, len(self.kdargs.paths)):
                for j in range(0, len(self.kdargs.paths[i])):
                    if np.all(self.kdargs.ureact[self.kdargs.paths[i][j], 1] == [1, 0]):
                        testlig = self.kdargs.ureact[self.kdargs.paths[i][j], 0][0]
                        onepath = ufirst[ulig.index(testlig)]
                        self.kdmap[self.kdargs.paths[i][j]] = onepath
                        self.kdargs.paths[i][j] = onepath
        # For a given protein number, all ligand binding reactions have the same KD
        if self.ligflag == 1:
            ulig = []
            ufirst = []
            for i in range(0, len(self.kdargs.paths)):
                for j in range(0, len(self.kdargs.paths[i])):
                    if np.all(self.kdargs.ureact[self.kdargs.paths[i][j], 1] == [0, 1]):
                        testlig = self.kdargs.ureact[self.kdargs.paths[i][j], 0][0]
                        if testlig not in ulig:
                            ulig.append(testlig)
                            ufirst.append(self.kdargs.paths[i][j])
            for i in range(0, len(self.kdargs.paths)):
                for j in range(0, len(self.kdargs.paths[i])):
                    if np.all(self.kdargs.ureact[self.kdargs.paths[i][j], 1] == [0, 1]):
                        testlig = self.kdargs.ureact[self.kdargs.paths[i][j], 0][0]
                        onepath = ufirst[ulig.index(testlig)]
                        self.kdmap[self.kdargs.paths[i][j]] = onepath
                        self.kdargs.paths[i][j] = onepath
        # For a given ligand binding reaction, the KD is the same across all protein states
        if self.ligflag == 3:
            ulig = []
            ufirst = []
            for i in range(0, len(self.kdargs.paths)):
                for j in range(0, len(self.kdargs.paths[i])):
                    if np.all(self.kdargs.ureact[self.kdargs.paths[i][j], 1] == [0, 1]):
                        testlig = self.kdargs.ureact[self.kdargs.paths[i][j], 0][1]
                        if testlig not in ulig:
                            ulig.append(testlig)
                            ufirst.append(self.kdargs.paths[i][j])
            for i in range(0, len(self.kdargs.paths)):
                for j in range(0, len(self.kdargs.paths[i])):
                    if np.all(self.kdargs.ureact[self.kdargs.paths[i][j], 1] == [0, 1]):
                        testlig = self.kdargs.ureact[self.kdargs.paths[i][j], 0][1]
                        onepath = ufirst[ulig.index(testlig)]
                        self.kdmap[self.kdargs.paths[i][j]] = onepath
                        self.kdargs.paths[i][j] = onepath

    def MakeFitGrid(self):
        kdargs = self.kdargs
        # kdargs.kds=[1E-8,1.0E-7,1.0E-5]
        out = [MinFree(kdargs.pconc[i], kdargs.lconc[i], kdargs.ureact, kdargs.nprottab, kdargs.nligtab, kdargs.paths,
                       kdargs.kds, kdargs.pfrees[i], kdargs.lfrees[i], kdargs.nfactors) for i in
               range(0, len(kdargs.lconc))]
        out = np.array(out)
        kdargs.pfrees = out[:, 0]
        kdargs.lfrees = out[:, 1]
        if np.any(kdargs.pfrees == 0) or np.any(kdargs.lfrees[1:] == 0):
            out = [
                MinFree(kdargs.pconc[i], kdargs.lconc[i], kdargs.ureact, kdargs.nprottab, kdargs.nligtab, kdargs.paths,
                        kdargs.kds, kdargs.pfrees[i], kdargs.lfrees[i], kdargs.nfactors) for i in
                range(0, len(kdargs.lconc))]
            out = np.array(out)
            kdargs.pfrees = out[:, 0]
            kdargs.lfrees = out[:, 1]

        out2 = [[MakeGrid(kdargs.pfrees[i], kdargs.lfrees[i], kdargs.ureact, kdargs.nprottab, kdargs.nligtab,
                          kdargs.paths, kdargs.kds, kdargs.nfactors)[0][row[0], row[1]] for row in kdargs.nodelist] for
                i in range(0, len(kdargs.lconc))]
        self.fitgrid = np.transpose(
            np.array([safedivide(np.array(out2[i]), np.sum(out2[i])) for i in range(0, len(kdargs.lconc))]))

    def OutlierTest(self):
        kdarr = self.randfit
        logarr = np.log10(kdarr)
        logmean = np.mean(logarr, axis=0)
        logstd = np.std(logarr, axis=0)
        logdiff = np.abs(logarr - logmean)
        logbool = logdiff < 2.5 * logstd
        normbool = np.all(np.transpose(logbool), axis=0)
        removed = self.randfit[np.logical_not(normbool)]
        if len(removed) > 0:
            print("Removing Outliers: ")
            print(removed)

            # print "Remaining Bootstrap Fits: "
            # print self.randfit[normbool]
        out = self.randfit[normbool]
        self.randfit = out

    def Bootstrap(self, std, numpts):
        # Makes an array of random mask arrays that weight the error calculation. Sampling by replacement
        hists = [np.reshape(np.histogram(
            np.array(np.random.random_integers(0, len(np.ravel(self.data)) - 1, size=len(np.ravel(self.data)))),
            bins=list(range(0, len(np.ravel(self.data)) + 1)))[0], np.shape(self.data)) for i in range(0, numpts)]
        print("Bootstrapping number of tests: ", numpts)
        startt = time.perf_counter()
        # TODO: Fix parallel processing in built version

        try:
            self.randfit = []
            queue = multiprocessing.Queue()
            results_queue = multiprocessing.Queue()
            [queue.put((self.data, self.kdargs, hists[i])) for i in range(0, numpts)]
            workers = [multiprocessing.Process(target=BootMinWorker, args=(queue, results_queue)) for i in
                       range(0, multiprocessing.cpu_count())]
            for p in workers:
                p.start()
            count = 0
            while True:
                if count == numpts:
                    break
                try:
                    out = results_queue.get(False)
                    self.randfit.append(out)
                    count += 1
                except Exception as e:
                    pass
            for p in workers:
                p.join()
            self.randfit = np.array(self.randfit)

        except Exception as e:
            print("Parallel Failed, using sequential...", e)

            self.randfit = np.array([BootMin(self.data, self.kdargs, hists[i]) for i in range(0, numpts)])

        print("Evaluation Time: ", time.perf_counter() - startt)
        if self.outlierflag:
            self.OutlierTest()
        means = np.mean(self.randfit, axis=0)
        stddevs = np.std(self.randfit, axis=0)
        return means, stddevs

    def RunBootstrap(self):
        # Find Confidence Interval
        if self.bootnum > 0:
            self.meanfit, self.stddevs = self.Bootstrap(self.stdres, self.bootnum)
            np.savetxt(self.header + "_kds.txt", np.array([np.concatenate(self.fit), self.stddevs, self.meanfit]))
        else:
            self.meanfit = None
            self.stddevs = None
            np.savetxt(self.header + "_kds.txt", np.array([np.concatenate(self.fit)]))

    def GraphPlot(self, ax=None):
        draw_graph(self.graph, self.ugraph, np.concatenate(self.fit), self.stddevs, self.kdmap, self.header, ax=ax,
                   ulabel=self.label)

    def PlotTrace(self, topax=None):
        dims = self.data.shape
        # colormap = cm.get_cmap(self.cmap, dims[0])
        colormap = mpl.colormaps[self.cmap].resampled(dims[0])
        colors = colormap(np.arange(dims[0]))
        if topax is None:
            fig = plt.figure()
            ax = fig.add_subplot(111)
        else:
            ax = topax
        for i in range(0, dims[0]):
            # print data[i]
            # print fitgrid[:,i]
            xvals = self.kdargs.lconc
            if np.all(xvals == 0):
                xvals = self.kdargs.pconc
            try:
                ax.plot(xvals, self.data[i], 'o', color=colors[i], linestyle="None", label=self.nodenames[i])
                ax.plot(xvals, self.fitgrid[i], '--', color=colors[i])
            except Exception as e:
                ax.plot(xvals, self.data[i], 'o', color=colors[i], linestyle="None", label="Unknown")
                ax.plot(xvals, self.fitgrid[i], '--', color=colors[i])

        ax.set_xlabel(self.label)
        if topax is None:
            plt.legend(bbox_to_anchor=(1.0, 0.5, 0.25, 0.5), loc=1, numpoints=1)
            plt.tight_layout(rect=(0.05, 0.05, 1.05, 1))
            plt.xlabel("Ligand Concentration")
            plt.ylabel("Relative Intensity")
            # plt.savefig(header+".png")

            plt.show()

    def PlotHist(self, topax=None):
        if topax is not None and self.randfit is not None:
            plot = topax
            if self.maxsites == 0:
                xlab = "Log Macro KD"
            else:
                xlab = "Log Micro KD"
            plot.histogram(np.log10(self.randfit.transpose()), labels=self.kdmap, xlab=xlab, ylab="histogram",
                           title="KD Distribution")

    def SetupModelAgg(self, **kwargs):
        # Define path favoring ligand addition->0 or protein addition->1
        if "agg" in kwargs:
            key = kwargs["agg"]
            # Key Should be free, parallel, one
            if key == "free":
                self.mode = 0
                self.protflag = 0
            elif key == "parallel":
                self.mode = 1
                self.protflag = 3
            elif key == "one":
                self.mode = 1
                self.protflag = 2
            elif key == "series":
                self.protflag = 1
                self.mode = 1

        if "lig" in kwargs:
            key = kwargs["lig"]
            # Key Should be free, parallel, one
            if key == "free":
                self.ligflag = 0
            elif key == "parallel":
                self.ligflag = 3
            elif key == "one":
                self.mode = 0
                self.ligflag = 2
            elif key == "series":
                self.ligflag = 1

        # Defining the potential nodes and reactions
        # Note: the nodes are species that are experimentally measured.
        self.reactions = []

        if self.hill:
            self.reactions.append([[0, 0, 1], [0, 0, 1], [0, 0, 2]])
        for i in range(1, self.maxligagg + 1):
            i = int(i)
            # if i + 1 <= self.maxligagg:
            #    self.reactions.append([[0, i, i], [0, 1, 1], [0, i + 1, i + 1]])
            for j in self.nlig:
                j = int(j)
                if i + 1 <= self.maxligagg and j >= i + 1:
                    self.reactions.append([[1, j, i], [0, 0, 1], [1, j, i + 1]])
                if j + 1 <= self.numtotlig and (j + 1 <= self.maxsites or self.maxsites == 0) and (
                        j >= i or i == 1):
                    self.reactions.append([[1, j, i], [0, 1, 0], [1, j + 1, i]])
        self.reactions = np.array(self.reactions)
        # print(self.reactions)

        # Fishing out degenerate KD values
        self.kdargs.degen = []
        for i in self.reactions:
            for j in self.reactions:
                if np.all(i[[2]] == j[[2]]):
                    if self.mode == 0:
                        if np.all(i[[1]] == [0, 1, 0]) and np.all(j[[1]] == [0, 0, 1]):
                            self.kdargs.degen.append([i, j])
                    else:
                        if np.all(i[[1]] == [0, 0, 1]) and np.all(j[[1]] == [0, 1, 0]):
                            self.kdargs.degen.append([i, j])
        self.kdargs.degen = np.array(self.kdargs.degen)
        # print("Degen:", self.kdargs.degen)

        # Deleting degenate KD values to get number of unique KDs
        self.kdargs.ureact = []
        for i in self.reactions:
            if len(self.kdargs.degen) != 0:
                if not np.any([np.all(i == j) for j in self.kdargs.degen[:, 1]]):
                    self.kdargs.ureact.append(i)
            else:
                self.kdargs.ureact.append(i)
        self.kdargs.ureact = np.array(self.kdargs.ureact)
        self.kdmap = list(range(0, len(self.kdargs.ureact)))

        # Finding the path to each KD. Will be used for calculating the grid.
        self.kdargs.paths = findpaths(self.kdargs.ureact)
        # print(self.kdargs.paths)
        self.ModPathsAgg()
        self.kdmap, self.kdargs.paths = fixlists(self.kdmap, self.kdargs.paths)
        self.numkd = len(np.unique(np.concatenate(self.kdargs.paths)))
        if len(self.kdargs.degen > 0):
            self.reactions = np.concatenate([self.kdargs.ureact, self.kdargs.degen[:, 1]])
            self.kdmap = np.concatenate(
                [np.array(self.kdmap), np.arange(self.numkd, self.numkd + len(self.kdargs.degen[:, 1]))])
        # Set up graph data for plot
        self.graph = make_graph_agg(self.reactions)
        self.ugraph = make_graph_agg(self.kdargs.ureact)
        self.nodenames = ["P" + str(int(self.kdargs.nodelist[i, 0])) + "L" + str(int(self.kdargs.nodelist[i, 1])) for i
                          in
                          range(0, len(self.kdargs.nodelist))]

        # Calculate Statistical Factors
        if self.maxsites != 0 and self.maxsites is not None:
            self.kdargs.maxsites = self.maxsites
            narr = np.array(self.kdargs.nprottab * self.maxsites)
            iarr = self.kdargs.nligtab
            niarr = np.clip(narr - iarr, 0, sys.maxsize)
            self.kdargs.nfactors = special.factorial(narr) / (
                    special.factorial(niarr) * special.factorial(iarr))
        else:
            self.kdargs.nfactors = None

        # Quick Plot to test structure. Note: exits the program without fitting.
        if "preview" in kwargs:
            preview = kwargs["preview"]
            if preview:
                draw_graph_structure(self.graph, self.ugraph, self.kdmap, ax=self.plot2)
                sys.exit()

    def ModPathsAgg(self):
        # All ligand binding reactions have the same KD
        if self.ligflag == 2:
            for i in range(0, len(self.kdargs.paths)):
                if (len(self.kdargs.paths[i]) == 1 and np.all(
                        self.kdargs.ureact[self.kdargs.paths[i][0], 1] == [0, 1, 0])):
                    onepath = self.kdargs.paths[i][0]
            for i in range(0, len(self.kdargs.paths)):
                for j in range(0, len(self.kdargs.paths[i])):
                    if np.all(self.kdargs.ureact[self.kdargs.paths[i][j], 1] == [0, 1, 0]):
                        self.kdmap[self.kdargs.paths[i][j]] = onepath
                        self.kdargs.paths[i][j] = onepath
        # All protein binding reactions have the same KD
        if self.protflag == 2:
            for i in range(0, len(self.kdargs.paths)):
                if (len(self.kdargs.paths[i]) == 1 and np.all(
                        self.kdargs.ureact[self.kdargs.paths[i][0], 1] == [0, 0, 1])):
                    onepath = self.kdargs.paths[i][0]
            for i in range(0, len(self.kdargs.paths)):
                for j in range(0, len(self.kdargs.paths[i])):
                    if np.all(self.kdargs.ureact[self.kdargs.paths[i][j], 1] == [0, 0, 1]):
                        self.kdmap[self.kdargs.paths[i][j]] = onepath
                        self.kdargs.paths[i][j] = onepath
        # For a given ligand number, all protein binding reactions have the same KD
        if self.protflag == 1:
            ulig = []
            ufirst = []
            for i in range(0, len(self.kdargs.paths)):
                for j in range(0, len(self.kdargs.paths[i])):
                    if np.all(self.kdargs.ureact[self.kdargs.paths[i][j], 1] == [0, 0, 1]):
                        testlig = self.kdargs.ureact[self.kdargs.paths[i][j], 1][1]
                        if testlig not in ulig:
                            ulig.append(testlig)
                            ufirst.append(self.kdargs.paths[i][j])
            for i in range(0, len(self.kdargs.paths)):
                for j in range(0, len(self.kdargs.paths[i])):
                    if np.all(self.kdargs.ureact[self.kdargs.paths[i][j], 1] == [0, 0, 1]):
                        testlig = self.kdargs.ureact[self.kdargs.paths[i][j], 1][1]
                        onepath = ufirst[ulig.index(testlig)]
                        self.kdmap[self.kdargs.paths[i][j]] = onepath
                        self.kdargs.paths[i][j] = onepath
        # For a given protein association reaction, the KD is the same across all ligand bound states
        if self.protflag == 3:
            ulig = []
            ufirst = []
            for i in range(0, len(self.kdargs.paths)):
                for j in range(0, len(self.kdargs.paths[i])):
                    if np.all(self.kdargs.ureact[self.kdargs.paths[i][j], 1] == [0, 0, 1]):
                        testlig = self.kdargs.ureact[self.kdargs.paths[i][j], 2][0]
                        if testlig not in ulig:
                            ulig.append(testlig)
                            ufirst.append(self.kdargs.paths[i][j])
            for i in range(0, len(self.kdargs.paths)):
                for j in range(0, len(self.kdargs.paths[i])):
                    if np.all(self.kdargs.ureact[self.kdargs.paths[i][j], 1] == [0, 0, 1]):
                        testlig = self.kdargs.ureact[self.kdargs.paths[i][j], 2][0]
                        onepath = ufirst[ulig.index(testlig)]
                        self.kdmap[self.kdargs.paths[i][j]] = onepath
                        self.kdargs.paths[i][j] = onepath
        # For a given protein number, all ligand binding reactions have the same KD
        if self.ligflag == 1:
            ulig = []
            ufirst = []
            for i in range(0, len(self.kdargs.paths)):
                for j in range(0, len(self.kdargs.paths[i])):
                    if np.all(self.kdargs.ureact[self.kdargs.paths[i][j], 1] == [0, 1, 0]):
                        testlig = self.kdargs.ureact[self.kdargs.paths[i][j], 0][0]
                        if testlig not in ulig:
                            ulig.append(testlig)
                            ufirst.append(self.kdargs.paths[i][j])
            for i in range(0, len(self.kdargs.paths)):
                for j in range(0, len(self.kdargs.paths[i])):
                    if np.all(self.kdargs.ureact[self.kdargs.paths[i][j], 1] == [0, 1, 0]):
                        testlig = self.kdargs.ureact[self.kdargs.paths[i][j], 0][0]
                        onepath = ufirst[ulig.index(testlig)]
                        self.kdmap[self.kdargs.paths[i][j]] = onepath
                        self.kdargs.paths[i][j] = onepath
        # For a given ligand binding reaction, the KD is the same across all protein states
        if self.ligflag == 3:
            ulig = []
            ufirst = []
            for i in range(0, len(self.kdargs.paths)):
                for j in range(0, len(self.kdargs.paths[i])):
                    if np.all(self.kdargs.ureact[self.kdargs.paths[i][j], 1] == [0, 1, 0]):
                        testlig = self.kdargs.ureact[self.kdargs.paths[i][j], 2][1]
                        if testlig not in ulig:
                            ulig.append(testlig)
                            ufirst.append(self.kdargs.paths[i][j])
            for i in range(0, len(self.kdargs.paths)):
                for j in range(0, len(self.kdargs.paths[i])):
                    if np.all(self.kdargs.ureact[self.kdargs.paths[i][j], 1] == [0, 1, 0]):
                        testlig = self.kdargs.ureact[self.kdargs.paths[i][j], 2][1]
                        onepath = ufirst[ulig.index(testlig)]
                        self.kdmap[self.kdargs.paths[i][j]] = onepath
                        self.kdargs.paths[i][j] = onepath


if __name__ == "__main__":
    multiprocessing.freeze_support()
    switch = 3
    if switch == 0:
        path = "C:\\cprog\\Tim\\raw"
        file_name = "man7.inp"
        os.chdir(path)
        data = np.loadtxt(file_name)
        data = np.transpose(data)
        pconc = data[0]
        lconc = data[1]
        data = data[2:]
        header = file_name.split(".")[0]
        nprot = 1
        nlig = 4
        maxsites = 4
        nodelist = [[1, j] for j in range(0, nlig + 1)]
    elif switch == 1:
        path = "C:\\cprog\\UniFit"
        file_name = "trimersim.txt"
        os.chdir(path)
        data = np.loadtxt(file_name)
        data = np.transpose(data)
        pconc = data[0]
        lconc = data[1]
        data = data[2:6]
        header = file_name.split(".")[0]
        nprot = 1
        nlig = 3
        maxsites = 3
        nodelist = [[1, j] for j in range(0, nlig + 1)]
    elif switch == 2:
        path = "C:\\cprog\\UniFit"
        file_name = "dimersim.txt"
        os.chdir(path)
        data = np.loadtxt(file_name)
        data = np.transpose(data)
        pconc = data[0]
        lconc = data[1]
        data = data[2:5]
        header = file_name.split(".")[0]
        nprot = 1
        nlig = 2
        maxsites = 2
        nodelist = [[1, j] for j in range(0, nlig + 1)]
    else:
        path = "C:\\cprog\\PaperData\\Ben_1D_2D"
        path = "Z:\\mtmarty\\Archive\\PostDoc\\LargeBackup\\Archive-Done\\PaperData\\Ben_1D_2D"
        file_name = "PeptideBinding1D.txt"
        os.chdir(path)
        data = np.loadtxt(file_name)
        dims = data.shape
        pconc = [5. for i in range(0, dims[1])]
        lconc = [4., 8., 16., 32., 64., 128.]
        header = file_name.split(".")[0]
        nprot = 2
        nlig = 2
        maxsites = 1
        nodelist = [[1, 0], [1, 1], [2, 0], [2, 1], [2, 2]]

    # KDmodel(None, None, None, None, None, 3, 3, preview=True, prot="free", lig="free")

    fit = KDmodel(data, pconc, lconc, nodelist, header, nprot, nlig, bootnum=0, maxsites=maxsites, prot="free",
                  lig="parallel")
    # fit.RunBootstrap()
