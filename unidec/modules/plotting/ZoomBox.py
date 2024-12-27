from matplotlib.patches import Rectangle
import wx
from unidec.modules.plotting.ZoomCommon import *


class ZoomBox(ZoomCommon):
    """
    Select a min/max range of the x axes for a matplotlib Axes

    Example usage::

        from matplotlib.widgets import  RectangleSelector
        from pylab import *

        def onselect(xmin, xmax, value, ymin, ymax):
          'eclick and erelease are matplotlib events at press and release'
          print ' x,y min position : (%f, %f)' % (xmin, ymin)
          print ' x,y max position   : (%f, %f)' % (xmax, ymax)
          print ' used button   : ', eclick.button

        def toggle_selector(event):
            print ' Key pressed.'
            if event.key in ['Q', 'q'] and toggle_selector.RS.active:
                print ' RectangleSelector deactivated.'
                toggle_selector.RS.set_active(False)
            if event.key in ['A', 'a'] and not toggle_selector.RS.active:
                print ' RectangleSelector activated.'
                toggle_selector.RS.set_active(True)

        x = arange(100)/(99.0)
        y = sin(x)
        fig = figure
        axes = subplot(111)
        axes.plot(x,y)

        toggle_selector.RS = ZoomBox(axes, onselect, drawtype='line')
        connect('key_press_event', toggle_selector)
        show()
    """

    def __init__(self, axes, onselect, parent=None, groups=None, drawtype='box',
                 minspanx=None,
                 minspany=None,
                 useblit=False,
                 lineprops=None,
                 rectprops=None,
                 onmove_callback=None,
                 spancoords='data',
                 button=None,
                 data_lims=None, swoop=False,
                 integrate=0, smash=0, pad=0.0001):

        """
        Create a selector in axes.  When a selection is made, clear
        the span and call onselect with

          onselect(pos_1, pos_2)

        and clear the drawn box/line. There pos_i are arrays of length 2
        containing the x- and y-coordinate.

        If minspanx is not None then events smaller than minspanx
        in x direction are ignored(it's the same for y).

        The rect is drawn with rectprops; default
          rectprops = dict(facecolor='red', edgecolor = 'black',
                           alpha=0.5, fill=False)

        The line is drawn with lineprops; default
          lineprops = dict(color='black', linestyle='-',
                           linewidth = 2, alpha=0.5)

        Use type if you want the mouse to draw a line, a box or nothing
        between click and actual position ny setting

        drawtype = 'line', drawtype='box' or drawtype = 'none'.

        spancoords is one of 'data' or 'pixels'.  If 'data', minspanx
        and minspanx will be interpreted in the same coordinates as
        the x and y axis, if 'pixels', they are in pixels

        button is a list of integers indicating which mouse buttons should
        be used for rectangle selection.  You can also specify a single
        integer if only a single button is desired.  Default is None, which
        does not limit which button can be used.
        Note, typically:
         1 = left mouse button
         2 = center mouse button (scroll wheel)
         3 = right mouse button
        """
        super(ZoomCommon, self).__init__()
        self.crossoverpercent = 0.06
        self.pad = pad
        self.parent = parent

        self.axes = None
        self.canvas = None
        self.visible = True
        self.cids = []
        self.mzz = None
        self.mzz2 = None
        self.mzcurve = None
        self.sarray = None
        self.swoop = swoop

        self.lflag = 0

        self.active = True  # for activation / deactivation
        self.to_draw = []
        self.background = None

        self.onselect = onselect
        self.onmove_callback = onmove_callback

        self.useblit = useblit
        self.minspanx = minspanx
        self.minspany = minspany

        self.integrate = integrate
        self.smash = smash
        self.comparemode = False
        self.comparexvals = []
        self.compareyvals = []

        if button is None or isinstance(button, list):
            self.validButtons = button
        elif isinstance(button, int):
            self.validButtons = [button]

        assert (spancoords in ('data', 'pixels'))

        self.spancoords = spancoords
        self.eventpress = None  # will save the data (position at mouseclick)
        self.eventrelease = None  # will save the data (pos. at mouserelease)

        self.new_axes(axes, rectprops)
        if data_lims is None:
            self.data_lims = GetStart(self.axes)
        else:
            self.data_lims = data_lims

        self.initialize()
        # self.zoomout()

    def initialize(self):
        # print("Data Lims:", self.data_lims)
        xmin, ymin, xmax, ymax = self.data_lims
        if xmin > xmax: xmin, xmax = xmax, xmin
        if ymin > ymax: ymin, ymax = ymax, ymin
        # assure that x and y values are not equal
        if xmin == xmax: xmax = xmin * 1.0001
        if ymin == ymax: ymax = ymin * 1.0001
        xspan = xmax - xmin
        yspan = ymax - ymin
        if xspan != 0:
            self.set_xlim(xmin - xspan * self.pad, xmax + xspan * self.pad, draw=False)
            # print("1")
        if yspan != 0:
            self.set_ylim(ymin - yspan * self.pad, ymax + yspan * self.pad, draw=False)

    def set_clipping(self):
        for axes in self.axes:
            ResetVisible(axes)
        # self.canvas.draw()

    def new_axes(self, axes, rectprops=None):
        self.axes = axes
        if self.canvas is not axes[0].figure.canvas:
            for cid in self.cids:
                self.canvas.mpl_disconnect(cid)

            self.canvas = axes[0].figure.canvas
            self.cids.append(self.canvas.mpl_connect('motion_notify_event', self.onmove))
            self.cids.append(self.canvas.mpl_connect('button_press_event', self.press))
            self.cids.append(self.canvas.mpl_connect('button_release_event', self.release))
            self.cids.append(self.canvas.mpl_connect('draw_event', self.update_background))
            # self.cids.append(self.canvas.mpl_connect('motion_notify_event', self.onmove))

        if rectprops is None:
            rectprops = dict(facecolor='white',
                             edgecolor='black',
                             alpha=0.5,
                             fill=False)
        self.rectprops = rectprops

        for axes in self.axes:
            self.to_draw.append(Rectangle((0, 0), 0, 1, visible=False, **self.rectprops))

        for axes, to_draw in zip(self.axes, self.to_draw):
            axes.add_patch(to_draw)

    def ignore(self, event):
        """return True if event should be ignored"""
        # If ZoomBox is not active :
        if not self.active:
            return True

        # If canvas was locked
        if not self.canvas.widgetlock.available(self):
            return True

        # Only do selection if event was triggered with a desired button
        if self.validButtons is not None:
            if event.button not in self.validButtons:
                if event.button == 3:
                    self.right_click(event)
                return True

        # If no button pressed yet or if it was out of the axes, ignore
        if self.eventpress is None:
            return event.inaxes not in self.axes

        # If a button pressed, check if the release-button is the same
        return (event.inaxes not in self.axes or
                event.button != self.eventpress.button)

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
        # Register a click if zoomout was not necessary
        if not zoomout and event is not None:
            if event.button == 1 and self.smash == 1:
                self.left_click(event.xdata, event.ydata)

        xspan = xmax - xmin
        yspan = ymax - ymin
        if xmin - xspan * self.pad != xmax + xspan * self.pad:
            self.set_xlim(xmin - xspan * self.pad, xmax + xspan * self.pad, draw=False)
            # print("2")
        self.set_ylim(ymin - yspan * self.pad, ymax + yspan * self.pad, draw=False)

        for axes in self.axes:
            ResetVisible(axes)
        self.kill_labels()
        self.canvas.draw()

    def press(self, event):
        """on button press event"""
        # Is the correct button pressed within the correct axes?
        if self.ignore(event): return

        self.buttonDown = True

        # make the drawed box/line visible get the click-coordinates,
        # button, ...
        for to_draw in self.to_draw:
            to_draw.set_visible(self.visible)
        self.eventpress = event
        if self.comparemode is True:
            self.comparexvals.append(event.xdata)
            self.compareyvals.append(event.ydata)

        if wx.GetKeyState(wx.WXK_SHIFT):
            # Get the current event x and y
            x, y = event.xdata, event.ydata
            self.mzz = [x, y]
            self.mzz2 = None
        elif wx.GetKeyState(wx.WXK_ESCAPE):
            # Get the current event x and y
            x, y = event.xdata, event.ydata
            self.mzz2 = [x, y]
            self.mzz = None
        else:
            self.mzz = None
            self.mzz2 = None
        return False

    def release(self, event):
        """on button release event"""
        if self.eventpress is None or (self.ignore(event) and not self.buttonDown):
            return
        # Do compare mode stuff
        self.buttonDown = False

        if self.sarray is not None and self.swoop:
            self.parent.on_swoop_drag(self.sarray)
            return

        if self.comparemode is True:
            if event.xdata is None:
                try:
                    self.comparexvals.append(self.prev[0])
                    self.compareyvals.append(self.prev[1])
                except Exception:
                    self.comparexvals.append(event.xdata)
                    self.compareyvals.append(event.ydata)
            else:
                self.comparexvals.append(event.xdata)
                self.compareyvals.append(event.ydata)
            self.eventpress = None  # reset the variables to their
            self.eventrelease = None  # inital values
            for to_draw in self.to_draw:
                to_draw.set_visible(False)
            return False
        else:

            # make the box/line invisible again
            for to_draw in self.to_draw:
                to_draw.set_visible(False)

            # left-click in place resets the x-axis or y-axis
            if self.eventpress.xdata == event.xdata and self.eventpress.ydata == event.ydata:

                if wx.GetKeyState(wx.WXK_CONTROL):
                    # Ignore the resize if the control key is down
                    if event.button == 1 and self.smash == 1:
                        self.left_click(event.xdata, event.ydata)
                    return
                self.zoomout(event)
                return

            self.canvas.draw()
            # release coordinates, button, ...
            self.eventrelease = event

            if self.spancoords == 'data':
                xmin, ymin = self.eventpress.xdata, self.eventpress.ydata
                # xmax, ymax = self.eventrelease.xdata, self.eventrelease.ydata

                # fix for if drag outside axes boundaries

                try:
                    offx = self.prev[0]
                    offy = self.prev[1]
                    xlims = self.axes[0].get_xlim()
                    ylims = self.axes[0].get_ylim()
                    if offx < xlims[0]:
                        offx = xlims[0]
                    elif offx > xlims[1]:
                        offx = xlims[1]

                    if offy < ylims[0]:
                        offy = ylims[0]
                    elif offy > ylims[1]:
                        offy = ylims[1]
                    '''
                    if np.abs(self.prev[0] - xlims[0]) < np.abs(self.prev[0] - xlims[1]):
                        # offx = xlims[0]
                        pass
                    else:
                        # offx = xlims[1]
                        pass

                    if np.abs(self.prev[1] - ylims[0]) < np.abs(self.prev[1] - ylims[1]):
                        # offy = ylims[0]
                        pass
                    else:
                        # offy = ylims[1]
                        pass'''
                except:
                    pass

                xmax, ymax = self.eventrelease.xdata or offx, self.eventrelease.ydata or offy
                # print(xmin, xmax)

            elif self.spancoords == 'pixels':
                xmin, ymin = self.eventpress.x, self.eventpress.y
                xmax, ymax = self.eventrelease.x, self.eventrelease.y
            else:
                raise ValueError('spancoords must be "data" or "pixels"')

            # assure that min<max values
            if xmin > xmax: xmin, xmax = xmax, xmin
            if ymin > ymax: ymin, ymax = ymax, ymin
            # assure that x and y values are not equal
            if xmin == xmax: xmax = xmin * 1.000001
            if ymin == ymax: ymax = ymin * 1.000001

            # Switch to span if a small delta y is used
            try:
                y0, y1 = event.inaxes.get_ylim()
            except Exception as e:
                y0, y1 = self.data_lims[1], self.data_lims[3]
                # print(e, xmin, xmax, ymin, ymax)
            if ymax - ymin < (y1 - y0) * self.crossoverpercent:
                # print(ymax,ymin,y0, y1, "Test")
                ymax = y1
                ymin = y0
                spanflag = 1
            else:
                spanflag = 0

            spanx = xmax - xmin
            spany = ymax - ymin
            xproblems = self.minspanx is not None and spanx < self.minspanx
            yproblems = self.minspany is not None and spany < self.minspany
            if xproblems or yproblems:
                """Box too small"""  # check if drawed distance (if it exists) is
                return  # not to small in neither x nor y-direction

            self.set_xlim(xmin, xmax, draw=False)
            # print("3")
            # print("xmin, xmax", xmin, xmax)
            # print("ymin, ymax", ymin, ymax)
            if spanflag:
                self.set_auto_ylim(xmin, xmax, draw=False)
            else:
                self.set_ylim(ymin, ymax, draw=False)

            self.set_clipping()
            self.kill_labels()
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

                self.onselect(xmin, xmax, value, ymin, ymax)  # zeros are for consistency with box zoom

            self.eventpress = None  # reset the variables to their
            self.eventrelease = None  # inital values
            return False

    def update(self):
        """draw using newfangled blit or oldfangled draw depending on useblit"""
        if self.useblit:
            if self.background is not None:
                self.canvas.restore_region(self.background)
            for axes, to_draw in zip(self.axes, self.to_draw):
                axes.draw_artist(to_draw)
            self.canvas.blit(self.canvas.figure.bbox)
        else:
            self.canvas.draw_idle()
        # self.set_clipping()
        return False

    def onmove(self, event):
        self.on_newxy(event.xdata, event.ydata)
        # newevent = NewXYEvent(xpos=event.xdata, ypos=event.ydata)
        # self.GetEventHandler().ProcessEvent(newevent)
        # 'on motion notify event if box/line is wanted'
        if self.eventpress is None or self.ignore(event):
            return
        x, y = event.xdata, event.ydata  # actual position (with
        #   (button still pressed)

        self.prev = x, y

        minx, maxx = self.eventpress.xdata, x  # click-x and actual mouse-x
        miny, maxy = self.eventpress.ydata, y  # click-y and actual mouse-y
        if minx > maxx: minx, maxx = maxx, minx  # get them in the right order
        if miny > maxy: miny, maxy = maxy, miny

        if self.mzz is not None and self.swoop:
            spany = maxy - miny
            minz = self.mzz[0] * self.mzz[1] / maxx
            maxz = self.mzz[0] * self.mzz[1] / minx
            zwidth = maxz - minz
            #print(self.mzz, spany, zwidth)
            if self.mzcurve is not None:
                self.mzcurve.remove()
            sarray = [self.mzz[0], self.mzz[1], spany, zwidth]
            self.mzcurve, self.sarray = self.parent.draw_mz_curve(sarray)
            self.canvas.draw()
            return
        elif self.mzz2 is not None and self.swoop:
            spany = maxy - miny
            roughmass = self.mzz2[0] * self.mzz2[1]
            minz = roughmass / maxx
            maxz = roughmass / minx
            minz = np.ceil(minz)
            maxz = np.floor(maxz)
            zwidth = maxz - minz
            #print(roughmass, spany, zwidth, minz, maxz)
            # print(self.mzz, spany, zwidth)
            if self.mzcurve is not None:
                self.mzcurve.remove()

            zmid = (minz + maxz) / 2
            mzmid = roughmass/zmid

            sarray = [mzmid, zmid, spany, zwidth]
            self.mzcurve, self.sarray = self.parent.draw_mz_curve(sarray)
            self.canvas.draw()
            return
        else:
            self.sarray = None

        # Changes from a yellow box to a colored line
        for axes in self.axes:
            y0, y1 = axes.get_ylim()
        if abs(maxy - miny) < abs(y1 - y0) * self.crossoverpercent:
            # print self.to_draw
            # print miny,maxy,y
            avg = (miny + maxy) / 2
            if y > miny:
                avg = miny
            else:
                avg = maxy
            miny = avg
            maxy = avg
            for to_draw in self.to_draw:
                to_draw.set_edgecolor('m')
                to_draw.set_alpha(0.9)
        else:
            for to_draw in self.to_draw:
                to_draw.set_edgecolor('k')
                to_draw.set_alpha(0.2)
        for to_draw in self.to_draw:
            to_draw.set_x(minx)  # set lower left of box
            to_draw.set_y(miny)
            to_draw.set_width(maxx - minx)  # set width and height of box
            to_draw.set_height(maxy - miny)

        value = 0.0
        if self.onmove_callback is not None and event.inaxes.lines != []:
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

            self.onmove_callback(minx, maxx, value, miny, maxy)  # zeros are for consistency with box zoom

        self.update()
        return False
