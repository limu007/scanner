from pyface.qt import QtGui, QtCore
import matplotlib

from traitsui.qt4.editor import Editor
from traitsui.qt4.basic_editor_factory import BasicEditorFactory
from traitsui.api  import Handler

from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt4 import NavigationToolbar2QT


class _MPLSimFigureEditor(Editor):
    scrollable  = True

    def init(self, parent):
        self.control = self._create_canvas(parent)
#		self.set_tooltip()

    def update_editor(self):
        pass

    def _create_canvas(self, parent):
        """ Create the MPL canvas. """
        # matplotlib commands to create a canvas
        frame = QtGui.QWidget()
        mpl_canvas = FigureCanvas(self.value)
        mpl_canvas.setParent(frame)
#		mpl_toolbar = NavigationToolbar2QT(mpl_canvas,frame)

        vbox = QtGui.QVBoxLayout()
        vbox.addWidget(mpl_canvas)
#		vbox.addWidget(mpl_toolbar)
        frame.setLayout(vbox)
        return frame

class MPLSimFigureEditor(BasicEditorFactory):
    klass = _MPLSimFigureEditor


class _MPLFigureEditor(Editor):
    scrollable  = True

    def init(self, parent):
        self.control = self._create_canvas(parent)
        self.set_tooltip()

    def update_editor(self):
        pass

    def _create_canvas(self, parent):
        """ Create the MPL canvas. """
        # matplotlib commands to create a canvas
        frame = QtGui.QWidget()
        mpl_canvas = FigureCanvas(self.value)
        mpl_canvas.setParent(frame)
        mpl_toolbar = NavigationToolbar2QT(mpl_canvas,frame)
        vbox = QtGui.QVBoxLayout()
        vbox.addWidget(mpl_canvas)
        vbox.addWidget(mpl_toolbar)
        frame.setLayout(vbox)
        return frame

from matplotlib.widgets import AxesWidget

class MPLFigureEditor(BasicEditorFactory):
    klass = _MPLFigureEditor

    def setup_mpl_events(self):
        self.image_axeswidget = AxesWidget(self.image_axes)
        self.image_axeswidget.connect_event('motion_notify_event', self.image_on_motion)
        #self.image_axeswidget.connect_event('figure_leave_event', self.on_cursor_leave)
        #self.image_axeswidget.connect_event('figure_enter_event', self.on_cursor_enter)
        wx.EVT_RIGHT_DOWN(self.image_figure.canvas, self.on_right_down)

    def on_right_down(self, event):
        if self.image_popup_menu is None:
            menu = wx.Menu()

    def image_on_motion(self, event):
        if event.xdata is None or event.ydata is None:
            return

#cid = fig.canvas.mpl_connect('button_press_event', onclick)

from traits.etsconfig.api import ETSConfig
ETSConfig.toolkit = 'qt4'

matplotlib.rcParams['backend.qt4']='PySide'
class MPLInitHandler(Handler):
    """Handler calls mpl_setup() to initialize mpl events"""

    def init(self, info):
        """This method gets called after the controls have all been
        created but before they are displayed.
        """
        info.object.mpl_setup()
        return True


from matplotlib.widgets import LassoSelector
from matplotlib.path import Path

import numpy as np 
class SelectFromCollection(object):
    """Select indices from a matplotlib collection using `LassoSelector`.

    Selected indices are saved in the `ind` attribute. This tool fades out the
    points that are not part of the selection (i.e., reduces their alpha
    values). If your collection has alpha < 1, this tool will permanently
    alter the alpha values.

    Note that this tool selects collection objects based on their *origins*
    (i.e., `offsets`).

    Parameters
    ----------
    ax : :class:`~matplotlib.axes.Axes`
        Axes to interact with.

    collection : :class:`matplotlib.collections.Collection` subclass
        Collection you want to select from.

    alpha_other : 0 <= float <= 1
        To highlight a selection, this tool sets all selected points to an
        alpha value of 1 and non-selec  ted points to `alpha_other`.
    """

    def __init__(self, ax, collection, alpha_other=0.3):
        

        self.canvas = ax.figure.canvas
        self.collection = collection
        self.alpha_other = alpha_other

        self.xys = collection.get_offsets() #works for scatter
        #self.xys = np.array(collection.get_data()).T # works for plot
        self.Npts = len(self.xys)

        # Ensure that we have separate colors for each object
        self.fc = collection.get_facecolors()
        if len(self.fc) == 0:
            raise ValueError('Collection must have a facecolor')
        elif len(self.fc) == 1:
            self.fc = np.tile(self.fc, (self.Npts, 1))

        self.lasso = LassoSelector(ax, onselect=self.onselect)
        self.ind = []

    def onselect(self, verts):
        path = Path(verts)
        self.ind = np.nonzero(path.contains_points(self.xys))[0]
        self.fc[:, -1] = self.alpha_other
        self.fc[self.ind, -1] = 1
        self.collection.set_facecolors(self.fc)
        self.canvas.draw_idle()

    def disconnect(self):
        self.lasso.disconnect_events()
        self.fc[:, -1] = 1
        self.collection.set_facecolors(self.fc)
        self.canvas.draw_idle()
