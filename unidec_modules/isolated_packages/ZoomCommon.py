import numpy as np
from matplotlib.transforms import Bbox
import matplotlib


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
            #line.set_clip_on(True)
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
            if x < xmax and x > xmin:
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
        if xmin == None or xmax == None:
            xmin = np.amin(np.ravel(xvals))
            xmax = np.amax(np.ravel(xvals))

        if xmin > xmax: xmin, xmax = xmax, xmin
        if ymin > ymax: ymin, ymax = ymax, ymin
        if xmin == xmax: xmax = xmin * 1.0001
        if ymin == ymax: ymax = ymin * 1.0001

        out = [xmin, ymin, xmax, ymax]
    else:
        out = axes.dataLim.bounds

    set_clipping(axes)

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

        oldpts = axes.bbox.get_points()
        oldpts[1, 1] *= 1.05
        newbbox = Bbox(oldpts)
        for o in axes.lines:
            #label = o.properties()["label"]
            marker = o.properties()["marker"]
            if marker != "None": #"child" in label:
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
