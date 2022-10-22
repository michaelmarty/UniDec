from matplotlib.transforms import blended_transform_factory
from matplotlib.patches import Rectangle
# from pubsub import setupkwargs
from pubsub import pub
from unidec_modules.isolated_packages.ZoomCommon import *

class NoZoomSpan:
    """
    Expansion of matplotlib embed in wx example by John Bender and Edward
    Abraham, see http://www.scipy.org/Matplotlib_figure_in_a_wx_panel

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
        if rectprops is None:
            rectprops = dict(facecolor='yellow', alpha=0.2)

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
                                       edgecolor="blue",
                                       **self.rectprops))

        if not self.useblit:
            for axes, rect in zip(self.axes, self.rect):
                axes.add_patch(rect)

    def update_background(self, event):
        'force an update of the background'
        if self.useblit:
            self.background = self.canvas.copy_from_bbox(self.canvas.figure.bbox)

    def ignore(self, event):
        'return True if event should be ignored'
        return event.inaxes not in self.axes or not self.visible or event.button != 1

    def press(self, event):
        'on button press event'
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

        'on button release event'
        if self.pressv is None or (self.ignore(event) and not self.buttonDown): return
        self.buttonDown = False

        vmin = self.pressv
        vmax = event.xdata or self.prev[0]

        if vmin > vmax: vmin, vmax = vmax, vmin
        span = vmax - vmin

        print(vmin, vmax, span)
        pub.sendMessage('scans_selected', min=vmin, max=vmax)
        if self.minspan is not None and span <= self.minspan:
            self.canvas.draw()
            for rect in self.rect:
                rect.set_x(vmin)
                rect.set_width(0.01)
            self.update()
            return

        self.pressv = None

        return False

    def update(self):
        'draw using newfangled blit or oldfangled draw depending on useblit'
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
        pub.sendMessage('newxy', xpos=event.xdata, ypos=event.ydata)
        'on motion notify event'
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
