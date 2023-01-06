import numpy as np
from matplotlib.transforms import Bbox, TransformedBbox


def GetMaxes(axes, xmin=None, xmax=None):
    yvals = []
    xvals = []
    for line in axes.lines:
        ydat = np.array(line.get_ydata())
        xdat = line.get_xdata()

        if xmin is not None and xmax is not None:
            bool1 = np.all(np.array([xdat > xmin, xdat < xmax]), axis=0)
            ydat = ydat[bool1]
            line.set_clip_on(True)
        else:
            pass
            # line.set_clip_on(True)
        try:
            yvals.append([np.amin(ydat), np.amax(ydat)])
            xvals.append([np.amin(xdat), np.amax(xdat)])
        except Exception as e:
            pass

    for p in axes.collections:
        try:
            xys = np.array(p.get_paths()[0].vertices)
        except IndexError:
            break
        offsets = p.get_offsets()
        if len(offsets) > 1:
            xys = np.array(offsets)
        ydat = xys[:, 1]
        xdat = xys[:, 0]
        if xmin is not None and xmax is not None:
            bool1 = np.all(np.array([xdat > xmin, xdat < xmax]), axis=0)
            ydat = ydat[bool1]
        try:
            yvals.append([np.amin(ydat), np.amax(ydat)])
            xvals.append([np.amin(xdat), np.amax(xdat)])
        except Exception as e:
            pass

    for patch in axes.patches:
        try:
            if (patch.get_width()) and (patch.get_height()):
                vertices = patch.get_path().vertices
                if vertices.size > 0:
                    xys = np.array(patch.get_patch_transform().transform(vertices))
                    ydat = xys[:, 1]
                    xdat = xys[:, 0]
                    if xmin is not None and xmax is not None:
                        bool1 = np.all(np.array([xdat > xmin, xdat < xmax]), axis=0)
                        ydat = ydat[bool1]
                    try:
                        yvals.append([np.amin(ydat), np.amax(ydat)])
                        xvals.append([np.amin(xdat), np.amax(xdat)])
                    except Exception as e:
                        pass
        except Exception as e:
            try:
                xys = patch.xy
                ydat = xys[:, 1]
                xdat = xys[:, 0]
                if xmin is not None and xmax is not None:
                    bool1 = np.all(np.array([xdat > xmin, xdat < xmax]), axis=0)
                    ydat = ydat[bool1]

                yvals.append([np.amin(ydat), np.amax(ydat)])
                xvals.append([np.amin(xdat), np.amax(xdat)])
            except Exception as e:
                pass

    for t in axes.texts:
        x, y = t.get_position()
        y = y * 1.01
        if xmin is not None and xmax is not None:
            if xmax > x > xmin:
                t.set_visible(True)
                yvals.append([y, y])
                xvals.append([x, x])
            else:
                t.set_visible(False)
        else:
            yvals.append([y, y])
            xvals.append([x, x])

    if len(yvals) != 0 and len(xvals) != 0:
        ymin = np.amin(np.ravel(yvals))
        ymax = np.amax(np.ravel(yvals))
        if xmin is None or xmax is None:
            xmin = np.amin(np.ravel(xvals))
            xmax = np.amax(np.ravel(xvals))

        if xmin > xmax: xmin, xmax = xmax, xmin
        if ymin > ymax: ymin, ymax = ymax, ymin
        if xmin == xmax: xmax = xmin * 1.0001
        if ymin == ymax: ymax = ymin * 1.0001

        out = [xmin, ymin, xmax, ymax]
    else:
        out = axes.dataLim.bounds

    # set_clipping(axes)

    return out


def set_clipping(axes):
    # print('Set Clipping')
    # clip_all(axes)
    # clip_none(axes)
    special_clip(axes)


def clip_all(axes):
    oldpts = axes.bbox.get_points()
    newbbox = Bbox(oldpts)
    axes.set_clip_box(newbbox)
    axes.set_clip_on(True)

    for line in axes.lines:
        line.set_clip_box(newbbox)
        line.set_clip_on(True)
        pass

    for o in axes.findobj():
        o.set_clip_box(newbbox)
        o.set_clip_on(True)


def special_clip(axes):
    if axes.get_clip_on():
        axes.set_clip_on(True)

        # oldpts = axes.bbox.get_points()
        # oldpts = axes.bbox.get_points()
        # oldpts[1, 1] *= 1.05
        # oldpts[0, 1] = 0.95
        # newbbox = Bbox(oldpts)
        newbbox = TransformedBbox(Bbox([[0.00, -0.05], [1.00, 1.05]]),
                                  axes.transAxes)
        # print(axes.bbox.get_points(), newbbox)
        for o in axes.lines:
            # label = o.properties()["label"]
            marker = o.properties()["marker"]
            if marker != "None":  # "child" in label:
                o.set_clip_box(newbbox)
                pass

    else:
        clip_none(axes)
        '''
        oldpts = axes.bbox.get_points()
        oldpts[1, 1] *= 1.05
        newbbox = Bbox(oldpts)
        for o in axes.lines:
            o.set_clip_box(newbbox)
            o.set_clip_on(True)
            
            # label = o.properties()["label"]
            marker = o.properties()["marker"]
            if marker != "None":  # "child" in label:
                o.set_clip_box(newbbox)
                o.set_clip_on(True)
            else:
                o.set_clip_on(False)
            pass'''


def clip_none(axes):
    axes.set_clip_on(False)

    for line in axes.lines:
        line.set_clip_on(False)

    for o in axes.findobj():
        o.set_clip_on(False)


def ResetVisible(axes):
    set_clipping(axes)
    for t in axes.texts:
        t.set_visible(True)


def GetStart(axes):
    outputs = np.array([GetMaxes(axis) for axis in axes])
    # print "Outputs",outputs
    xmin = np.amin(outputs[:, 0])
    ymin = np.amin(outputs[:, 1])
    xmax = np.amax(outputs[:, 2])
    ymax = np.amax(outputs[:, 3])

    if xmin > xmax: xmin, xmax = xmax, xmin
    if ymin > ymax: ymin, ymax = ymax, ymin
    if xmin == xmax: xmax = xmin * 1.0001
    if ymin == ymax: ymax = ymin * 1.0001

    out = [xmin, ymin, xmax, ymax]
    return out


class ZoomCommon:
    def __init__(self):
        self.active = True
        self.lines = []
        self.texts = []
        self.lflag = 0
        self.axes = []
        self.canvas = None
        self.background = None
        self.useblit = True

    def switch_label(self):
        self.lflag = (self.lflag + 1) % 2
        print("Plot Labels Switched to:", self.lflag)

    def label(self):
        self.texts = []
        self.lines = []
        for line in self.axes[0].lines:
            ydat = np.array(line.get_ydata())
            xdat = line.get_xdata()

            x0, x1 = self.axes[0].get_xlim()

            bool1 = x0 < xdat
            bool2 = xdat < x1
            bool3 = np.logical_and(bool1, bool2)

            cutdat = np.transpose([xdat[bool3], ydat[bool3]])
            if len(cutdat) > 0:
                ymin = np.amin(cutdat[:, 1])
                ymax = np.amax(cutdat[:, 1])

                window = len(cutdat) / 25.
                thresh = 0.05 * (ymax - ymin) + ymin

                peaks = self.peakdetect(cutdat, window=window, threshold=thresh)

                for p in peaks:
                    val = round(p[0], 3)
                    text = self.axes[0].text(p[0], p[1] + 0.1 * ymax, str(val), horizontalalignment="center",
                                             verticalalignment="top")
                    line = self.axes[0].vlines(p[0], p[1], p[1] + 0.05 * ymax, linestyles="dashed")
                    self.texts.append(text)
                    self.lines.append(line)
            self.canvas.draw()
            return

    def kill_labels(self):
        try:
            for t in self.texts:
                t.remove()
            for l in self.lines:
                l.remove()
            self.texts = []
            self.lines = []
            self.canvas.draw()
        except:
            pass
        if self.lflag == 1:
            self.label()

    def peakdetect(self, data, window=10., threshold=0.):
        peaks = []
        length = len(data)
        # maxval = np.amax(data[:, 1])
        for i in range(1, length):
            if data[i, 1] > threshold:  # Note this is different
                start = i - window
                end = i + window
                if start < 0:
                    start = 0
                if end > length:
                    end = length
                testmax = np.amax(data[int(start):int(end) + 1, 1])
                if data[i, 1] == testmax and data[i, 1] != data[i - 1, 1]:
                    peaks.append([data[i, 0], data[i, 1]])
        return np.array(peaks)

    def set_active(self, active):
        """ Use this to activate / deactivate the RectangleSelector

            from your program with an boolean variable 'active'.
        """
        self.active = active

    def get_active(self):
        """ to get status of active mode (boolean variable)"""
        return self.active

    def set_manual(self, xmin, xmax):
        for axes in self.axes:
            axes.set_xlim((xmin, xmax))
            if True:
                xmin, ymin, xmax, ymax = GetMaxes(axes, xmin=xmin, xmax=xmax)
                axes.set_ylim((ymin, ymax))
        self.canvas.draw()

    def set_manual_y(self, ymin, ymax):
        for axes in self.axes:
            axes.set_ylim((ymin, ymax))
        self.canvas.draw()

    def update_background(self, event):
        """force an update of the background"""
        if self.useblit:
            self.background = self.canvas.copy_from_bbox(self.canvas.figure.bbox)
