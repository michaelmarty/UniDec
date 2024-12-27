from matplotlib.transforms import blended_transform_factory
from matplotlib.patches import Rectangle

from unidec.modules.plotting.ZoomCommon import *


class ZoomSpan(ZoomCommon):
    """
    Expansion of matplotlib embed in wx example by John Bender and Edward
    Abraham, see https://www.scipy.org/Matplotlib_figure_in_a_wx_panel

    This version allows the user to zoom in on the figure using either
    a span selector or a box selector. You can also set a persistent span
    selector that acts as cursor references on top of whatever is plotted

    ZoomSpan based on matplotlib.widgets.SpanSelector
    CursorSpan based on matplotlib.widgets.SpanSelector
    BoxZoom based on matplotlib.widgets.RectangleSelector

    Brian J. Soher, Duke University, 20 October, 2010


    Select a min/max range of the x or y axes for a matplotlib Axes

    Example usage:

      axes = subplot(111)
      axes.plot(x,y)

      def onselect(vmin, vmax):
          print vmin, vmax
      span = ZoomSpan(axes, onselect, 'horizontal')

      onmove_callback is an optional callback that will be called on mouse move
      with the span range

    """

    def __init__(self, axes, onselect,
                 parent=None,
                 minspan=None,
                 useblit=False,
                 rectprops=None,
                 onmove_callback=None):
        """
        Create a span selector in axes.  When a selection is made, clear
        the span and call onselect with

          onselect(vmin, vmax)

        If minspan is not None, ignore events smaller than minspan

        The span rect is drawn with rectprops; default
          rectprops = dict(facecolor='red', alpha=0.5)

        set the visible attribute to False if you want to turn off
        the functionality of the span selector


        """
        super(ZoomCommon, self).__init__()
        if rectprops is None:
            rectprops = dict(facecolor='yellow', alpha=0.2)

        self.parent = parent
        self.pad = 0.0001
        self.axes = None
        self.canvas = None
        self.visible = True
        self.cids = []

        self.rect = []
        self.background = None
        self.pressv = None

        self.rectprops = rectprops
        self.onselect = onselect
        self.onmove_callback = onmove_callback
        self.useblit = useblit
        self.minspan = minspan

        # Needed when dragging out of axes
        self.buttonDown = False
        self.prev = (0, 0)

        self.new_axes(axes)
        self.data_lims = GetStart(self.axes)

    def new_axes(self, axes):
        self.axes = axes
        if self.canvas is not axes[0].figure.canvas:
            for cid in self.cids:
                self.canvas.mpl_disconnect(cid)

            self.canvas = axes[0].figure.canvas

            self.cids.append(self.canvas.mpl_connect('motion_notify_event', self.onmove))
            self.cids.append(self.canvas.mpl_connect('button_press_event', self.press))
            self.cids.append(self.canvas.mpl_connect('button_release_event', self.release))
            self.cids.append(self.canvas.mpl_connect('draw_event', self.update_background))

        for axes in self.axes:
            trans = blended_transform_factory(axes.transData, axes.transAxes)
            self.rect.append(Rectangle((0, 0), 0, 1,
                                       transform=trans,
                                       visible=False,
                                       **self.rectprops))

        if not self.useblit:
            for axes, rect in zip(self.axes, self.rect):
                axes.add_patch(rect)

    def ignore(self, event):
        """return True if event should be ignored"""
        return event.inaxes not in self.axes or not self.visible or event.button != 1

    def press(self, event):
        """on button press event"""
        if self.ignore(event): return
        self.buttonDown = True

        for axes, rect in zip(self.axes, self.rect):
            if rect in axes.patches:
                axes.patches.remove(rect)
                self.canvas.draw()

        for rect in self.rect:
            rect.set_visible(self.visible)
        self.pressv = event.xdata
        return False

    def release(self, event):
        """on button release event"""
        if self.pressv is None or (self.ignore(event) and not self.buttonDown): return
        self.buttonDown = False

        for rect in self.rect:
            rect.set_visible(False)

        # left-click in place resets the x-axis
        if event.xdata == self.pressv:
            # x0,y0,x1,y1=GetMaxes(event.inaxes)
            x0, y0, x1, y1 = self.data_lims
            for axes in self.axes:
                axes.set_xlim(x0, x1)
                axes.set_ylim(y0, y1 + y1 * 0.03)
                ResetVisible(axes)
            self.canvas.draw()
            return

        vmin = self.pressv
        vmax = event.xdata or self.prev[0]

        if vmin > vmax: vmin, vmax = vmax, vmin
        span = vmax - vmin
        if self.minspan is not None and span < self.minspan: return

        for axes in self.axes:
            # axes.set_xlim((self.pressv, event.xdata))
            axes.set_xlim((vmin, vmax))
            # Autoscale Y
            xmin, ymin, xmax, ymax = GetMaxes(axes, xmin=vmin, xmax=vmax)
            axes.set_ylim((ymin, ymax))
        self.canvas.draw()

        value = 0.0
        if self.onselect is not None and event.inaxes.lines != []:
            # gather the values to report in a selection event
            value = []
            x0, y0, x1, y1 = event.inaxes.dataLim.bounds
            dat = event.inaxes.lines[0].get_ydata()
            npts = len(dat)
            indx = int(round((npts - 1) * (event.xdata - x0) / (x1 - x0)))
            if indx > (npts - 1): indx = npts - 1
            if indx < 0: indx = 0
            for line in event.inaxes.lines:
                dat = line.get_ydata()
                if indx < len(dat):
                    value.append(dat[indx])
            if value == []: value = 0.0

            self.onselect(vmin, vmax, value, 0, 0)  # zeros are for consistency with box zoom

        self.pressv = None
        return False

    def update(self):
        """draw using newfangled blit or oldfangled draw depending on useblit"""
        if self.useblit:
            if self.background is not None:
                self.canvas.restore_region(self.background)
            for axes, rect in zip(self.axes, self.rect):
                axes.draw_artist(rect)
            self.canvas.blit(self.canvas.figure.bbox)
        else:
            self.canvas.draw_idle()

        return False

    def onmove(self, event):
        self.on_newxy(event.xdata, event.ydata)

        #'on motion notify event'
        if self.pressv is None or self.ignore(event):
            return
        x, y = event.xdata, event.ydata
        self.prev = x, y

        minv, maxv = x, self.pressv
        if minv > maxv: minv, maxv = maxv, minv
        for rect in self.rect:
            rect.set_x(minv)
            rect.set_width(maxv - minv)

        value = 0.0
        if self.onmove_callback is not None and event.inaxes.lines != []:
            vmin = self.pressv
            value = []
            vmax = event.xdata or self.prev[0]
            x0, y0, x1, y1 = event.inaxes.dataLim.bounds
            dat = event.inaxes.lines[0].get_ydata()
            npts = len(dat)
            indx = int(round((npts - 1) * (event.xdata - x0) / (x1 - x0)))
            if indx > (npts - 1): indx = npts - 1
            if indx < 0: indx = 0
            for line in event.inaxes.lines:
                dat = line.get_ydata()
                if indx < len(dat):
                    value.append(dat[indx])
            if value == []: value = 0.0

            if vmin > vmax: vmin, vmax = vmax, vmin
            self.onmove_callback(vmin, vmax, value, 0, 0)  # zeros are for consistency with box zoom

        self.update()
        return False

    def zoomout(self, event=None):
        # print("Zoom Out")
        # x0,y0,x1,y1=GetMaxes(event.inaxes)
        # print GetMaxes(event.inaxes)

        xmin, ymin, xmax, ymax = self.data_lims
        if xmin > xmax: xmin, xmax = xmax, xmin
        if ymin > ymax: ymin, ymax = ymax, ymin
        # assure that x and y values are not equal
        if xmin == xmax: xmax = xmin * 1.01
        if ymin == ymax: ymax = ymin * 1.01

        # Check if a zoom out is necessary
        zoomout = False
        for axes in self.axes:
            if axes.get_xlim() != (xmin, xmax) and axes.get_ylim() != (ymin, ymax):
                zoomout = True

                xspan = xmax - xmin
                yspan = ymax - ymin
                if xmin - xspan * self.pad != xmax + xspan * self.pad:
                    axes.set_xlim(xmin - xspan * self.pad, xmax + xspan * self.pad)
                axes.set_ylim(ymin - yspan * self.pad, ymax + yspan * self.pad)
                ResetVisible(axes)
        self.canvas.draw()
