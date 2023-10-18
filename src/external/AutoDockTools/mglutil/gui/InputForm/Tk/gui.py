#############################################################################
#
# Author: Michel F. SANNER, Sophie I COON, Ruth HUEY
#
# Copyright: M. Sanner TSRI 2000
#
#############################################################################

# $Header: /mnt/raid/services/cvs/python/packages/share1.5/mglutil/gui/InputForm/Tk/gui.py,v 1.52.4.2 2016/05/12 23:06:49 mgltools Exp $
#
# $Id: gui.py,v 1.52.4.2 2016/05/12 23:06:49 mgltools Exp $
#      

#
"""
This module implements a bunch of objects useful for user interaction
"""
#from DejaVu.extendedSlider import ExtendedSlider
#from DejaVu.Slider import Slider
import tkinter
import tkinter.scrolledtext
import Pmw
import types
import warnings
from mglutil.util.packageFilePath import findFilePath
from mglutil.util.misc import ensureFontCase, isInstance

def autoPlaceWidget(master, relx=0.5, rely=0.3):
    """Auto center the window rooted at master"""
    
    if master.winfo_ismapped():
        m_width = master.winfo_width()
        m_height = master.winfo_height()
        m_x = master.winfo_rootx()
        m_y = master.winfo_rooty()
    else:
        m_width = master.winfo_screenwidth()
        m_height = master.winfo_screenheight()
        m_x = m_y = 0
    w_width = master.winfo_reqwidth()
    w_height = master.winfo_reqheight()
    x = m_x + (m_width - w_width) * relx
    y = m_y + (m_height - w_height) * rely
    if x+w_width > master.winfo_screenwidth():
        x = master.winfo_screenwidth() - w_width
    elif x < 0:
        x = 0
    if y+w_height > master.winfo_screenheight():
        y = master.winfo_screenheight() - w_height
    elif y < 0:
        y = 0
    master.geometry("+%d+%d" % (x, y))


class AutoPlaceWidget:
    """defines method _set_transient used to place a widget automatically
       to the center of the screen"""
    
    def _set_transient(self, master, relx=0.5, rely=0.3):
        widget = self.root
        widget.withdraw() # Remain invisible while we figure out the geometry
        widget.transient(master)
        widget.update_idletasks() # Actualize geometry information
        if master.winfo_ismapped():
            m_width = master.winfo_width()
            m_height = master.winfo_height()
            m_x = master.winfo_rootx()
            m_y = master.winfo_rooty()
        else:
            m_width = master.winfo_screenwidth()
            m_height = master.winfo_screenheight()
            m_x = m_y = 0
        w_width = widget.winfo_reqwidth()
        w_height = widget.winfo_reqheight()
        x = m_x + (m_width - w_width) * relx
        y = m_y + (m_height - w_height) * rely
        if x+w_width > master.winfo_screenwidth():
            x = master.winfo_screenwidth() - w_width
        elif x < 0:
            x = 0
        if y+w_height > master.winfo_screenheight():
            y = master.winfo_screenheight() - w_height
        elif y < 0:
            y = 0
        widget.geometry("+%d+%d" % (x, y))
        widget.deiconify() # Become visible at the desired location

##             self.info()


class ListChooser:
    """class to present a list of objects in a scrolled listbox.
    entries are given as a list of items associated with a description.
    double clicking on an entry display the associated comment"""
    
    def __init__(self, root, mode='single', title='Choose',
                 entries = None, width=20, defaultValue=None,
                 lbpackcfg={'fill':'x', 'expand':0}, lbwcfg={}):

        assert mode in ['single', 'browse', 'multiple', 'extended' ]
        # put everything inside a Frame
        self.top = tkinter.Frame(root, bd=4)


        # add a scrolled text widget to display comments
        self.comments = tkinter.scrolledtext.ScrolledText(self.top, width=width,
                                                  height=4)
        self.comments.pack(fill='both', expand=0)

        # add a title
        self.title = tkinter.Label(self.top,text=title)
        self.title.pack(side='top', anchor='n')

        # build the scroll list
        f = tkinter.Frame(self.top)
        self.mode = mode
        scrollbar = tkinter.Scrollbar(f, orient='vertical')
        lbwcfg['selectmode'] = mode
        lbwcfg['yscrollcommand'] = scrollbar.set
        self.lb = tkinter.Listbox(*(f, ), **lbwcfg)
        scrollbar.config(command=self.lb.yview)
        scrollbar.pack(side='right', fill='y')
        self.lb.pack(*(), **lbpackcfg)
        f.pack(fill='x', expand=0)
        self.hasInfo = 0

        self.entries = []
        # add entries
        if entries:
            for e in entries:
                if type(e)==list or type(e)==tuple:
                    self.add(e[0], e[1])
                else:
                    self.add(e)

        if self.hasInfo==0:
            self.comments.forget()

        if defaultValue:
            if self.hasInfo:
                assert defaultValue in [x[0] for x in entries]
                # get the entry of the listBox corresponding to the defaultValue.
                defaultEntry = list(filter(lambda x, df = defaultValue:
                                  x[0] == df, entries))[0]
            else:
                assert defaultValue in [x for x in entries]
                # get the entry of the listBox corresponding to the defaultValue.
                defaultEntry = list(filter(lambda x, df = defaultValue:
                                  x == df, entries))[0]
            # Get its index
            indexDefaultValue = entries.index(defaultEntry)
            # Set the the selection at this index using the select_set of the
            # listBox widget.
            self.lb.select_set(indexDefaultValue)

    def add(self, entry, comment=None):
        """add an entry to the list"""
        self.entries.append( (entry, comment) )
        if type(entry) is not bytes:
            entry = repr(entry)
        self.lb.insert('end', entry)
        if not comment is None:
            self.lb.bind("<ButtonRelease-1>", self.info, '+')
            self.hasInfo=1
            self.comments.pack()

            
    def insert(self, pos, entry, comment=None):
        """insert an entry to the list"""
        #assert pos < len(self.entries)
        if pos == 'end':
            self.entries.append( (entry, comment) )
        else:
            self.entries.insert(pos, (entry, comment))
        if type(entry) is not bytes:
            entry = repr(entry)
        self.lb.insert(pos, entry)


    def info(self, event=None):
        """dispaly information about the curent selection"""
        self.comments.delete(0.0, 'end')
        s=self.entries[list(map( int, self.lb.curselection()))[0]][1]
        self.comments.insert('end', s+'\n')


    def clear(self):
        """clear all entries"""
        #self.lb.delete(0.0, 'end')
        self.lb.delete(0, 'end')
        self.entries = []

        
    def clearComments(self):
        """clear all entries"""
        #self.comments.delete(0.0, 'end')
        self.comments.delete(0, 'end')

        
    def remove(self, entry):
        """remove an entry"""
        ind = [x[0] for x in self.entries].index(entry)
        self.entries.remove(self.entries[ind])
        self.lb.delete( ind )

    
    def get(self, event=None):
        """get the current selection"""
        res = []
        for ent in map( int, self.lb.curselection() ):
            res.append( self.entries[ ent ][0] )
        return res

    
    def select(self, first, last=None):
        """add an entry to the current selection"""
        if last and self.mode!='single':
            self.lb.select_set(first, last)
        else:
            self.lb.select_set(first)
            self.info()







from tkinter.simpledialog import SimpleDialog
import string
        
##  class InputFormItem:
##      """Class to describe things that can be put inside an inputForm"""

##      def __init__(self, tkwcfg={}, tkpackOptions={}, tkpackMode=None):

##          self.master = None
##          self.type = Tkinter.Entry
##          self.name = 'NoName'
##          self.tkwcfg = tkwcfg
##          self.tkpackOptions = tkpackOptions
##          if tkpackMode==None: self.tkpacker = Tkinter.grid
##          else: self.tkpackMode = tkpacker


##  class IFentry(InputFormItem):
##      pass
##  class IFradiobutton(InputFormItem):
##      pass
##  class IFcheckbutton(InputFormItem):
##      pass
##  class IFbutton(InputFormItem):
##      pass
##  class IFframe(InputFormItem):
##      pass
##  class IFscrolledText(InputFormItem):
##      pass

##  class IFlistChooser(InputFormItem):
##      def __init__(self, mode='single', title='choose one', entries=[],
##                   tkwcfg={}, tkpackOptions={}, tkpackMode=None):

##          InputFormItem.__init__(self, tkwcfg=tkwcfg,
##                                 tkpackOptions=tkpackOptions,
##                                 tkpackMode=tkpackMode)
##          self.mode = mode
##          self.title = title
##          self.entries = entries

##  class IFunctionButton(InputFormItem):
##      pass


from collections import UserList

class InputFormDescr(UserList):
    """
    InputFormDescr is a list of dictionary describing the widgets to be placed
    in an inputform.
    Such an object can be passed to an InputForm object to actually build an
    InputForm

    The dictionnary can contain the following key:value couples:
    * widgetType: This specifies the type of the widget. It is usually a class
        name that can be any of the following

        o any Tkinter widget (Tkinter.Button etc....)
        o any CustomizedWidgets widgets (ListChooser, Dial etc....)
        o any of the Pmw widgets (Pmw.EntryFields....)
    For examples go to any Pmv commands.
    The dictionary of a widget can hold the following (keys:value) couples:

    * name      : string to specify a name associated to the widget. This name
        is then used to create the value dictionary.

    * tooltip   : string describing the widget which will be displayed in
                  a balloon when the mouse is over it.
    * wcfg      : The widget configuration dictionary must contain all the
        couple keyword:value you want to use to describe your widget. 
        For a list of the option available to describe widgets:  
        o Tkinter widgets options
        http://www.pythonware.com/library/tkinter/introduction/index.htm
        o Pmw widgets options
        http://pmw.sourceforge.net/doc/refindex.html
        o Tix widgets options
        http://www.python.org/doc/current/lib/module-Tix.html
        o Customized widgets options:
        http://www.scripps.edu/~sanner/python/inputform/customizedWidgets.html

    * gridcfg    : The grid configuration dictionary specifies where in the
        form the current widget will be placed.
        It contains keyword:value where the keys are the option of the
        Tkinter grid manager. The value of the row and column can be negative
        if you want to place your widget several row or column before the last
        one. The position of the widget will be found by the function
        findGridPosition. If no row or column is specified it will place your
        widget after the last one according to the specified direction you
        want your inputform to be built.You can't leave an empty space between
        two widgets.

    * container  : dictionary specifying the current widget or part of the
        current widget as a container.
        o key :(string) name of the container. This name will then be used
        to specify their parent
        o value : instance of the container or a string that will be evaluated
        and the result of which will be the container.
        'container':{'Page1':'w.page('Page1')},
        example:
        Notebook widget with two pages Page1 and Page2
        'container':{'Page1':'w.page('Page1')',
                     'Page2':'w.page('Page2')'
        Frame
        'container':{'myFrame':'w'}

        (see colorCommands.py, APBSCommand.py, fileCommands for example)  
    * parent    : allows to specify in which container to pack the current
      widget.(see colorCommands.py, APBSCommand.py, fileCommands for example) 
        
    * componentcfg: dictionnary allowing the developer to describe the
        components of the current widget.
        o key   : is the name of the component
        o value : description of the widget
        (see APBSCommand for an example)

    * defaultValue: The given defaultValue is used to initialize the widget
        if it has a bound variable or if the widget has a set method.

    * type        :The keyword type is used to check the validity of the data
        the user gives (it is a string by default).

    * required    : (1 or 0) The keyword required is flag to specify
        whether or not this field is required.

    example:
    =========================================================================

        # Build a Tkinter.Button with the text 'hello', the text will be in
        # red, using the font times size 20 and bold, the background will
        # be blue when not activated and green when activated.
        # The foreground will be white when desactivated. The button's relief
        # will be sunken and its borderwidth will be 8 pixel.
        idf = InputFormDescr(title='Button Example')
        idf.append({'widgetType':Tkinter.Button,
                    'name':'helloBut',
                    'wcfg':{'text':'Hello','foreground':'White',
                            'background' : 'Blue',
                            'relief' : Tkinter.SUNKEN,
                            'font':(ensureFontCase('times'),20, 'bold'), 'borderwidth' : 8,
                            'activebackground' : 'Green'},
                     'gridcfg':{'sticky':'we'}})

    
    """
    def __init__(self, title='Form'):
        UserList.__init__(self)
        self.title = title
        self.entryByName = {} # reverse lookup to find an item using a name

class CallBackFunction:
    """
Class to allow to specify arguments to a callback function bound to a widget
"""

    def __init__(self, function, *args, **kw):
        self.function = function
        self.args = args
        self.kw = kw

    def __call__(self, *args, **kw):
        args = self.args + args
        kw.update(self.kw)
        self.function(*args, **kw)

def getkw(kw, name):
    """ function to get a dictionary entry and remove it from the dictionary"""

    v = kw.get(name)
    if name in list(kw.keys()): del kw[name]
    return v


class InputForm(AutoPlaceWidget):
    
    """
Class to create a inputForm to collect values from users and
return these values.
More information and examples can be found at
http://www.scripps.edu/~sanner/software/inputform/tableOfContent.html

required arguments:
   master  -- argument to specify the master of the widget to be created
              which means that the current form will be a 'slave' of the
              given master. MS april 2012: this seems to be the root window,
              i.e. top window. not really needed I think

   root    --  If None is specified then a Tkinter.TopLevel object is
               created. else it is the frame into which the form will be
               created and displayed

   descr   -- (instance of InputFormDescr) which is a list of dictionary
               describing the widgets to be added to the form.
               see the documentation string for InputFormDescr

optional arguments :

    modal  -- Flag specifying if the form is modal or not. When a form
      is modal it grabs the focus and only releases it when the
      form is dismissed. When a form is modal an OK and a CANCEL
      button will be automatically added to the form.
      (default = 1)

    blocking -- Flag specifying if the form is blocking or not. When set to
      1 the form is blocking and the calling code will be stopped until the
      form is dismissed. An OK and a CANCEL button will be automatically
      added to the form. (default = 0)

    defaultDirection -- ('row', 'col') specifies the direction in
      which widgets are gridded  into the form by default. (default='row')

    closeWithWindow -- Flag specifying whether or not the form should be
      minimized/maximized  when the master window is. (default=1)

    okCfg -- dictionnary specifying the configuration of the OK button.
             if a callback function is specified using the keyword
             command this callback will be added to the default callback
             Ok_cb

    cancelCfg -- dictionnary specifying the configuration of the CANCEL button
             if a callback function is specified using the keyword
             command this callback will be added to the default callback
             Ok_cb

    initFunc  -- specifies a function to initialize the form.

    onDestroy -- specifies a function to be called when using the close
       widget of a window.

    okcancel -- Boolean Flag to specify whether or not to create the OK and CANCEL
                button.
    scrolledFrame -- Flag when set to 1 the main frame is a scrollable frame
                     else it is static Frame (default 0)
    width -- specifies the width of the main frame (400)
    height -- specifies the height of the main frame. (200)
    help   -- specifies the web adress to a help page. If this is provided
              a Help (?) button will be created which will open a
              web browser to the given adress.

    If the list of argument of the constructor changes please update
    the showForm method of the ViewerFramework.VFCommand
    class to present a form for the user to fill out.
    After building the form the entryByName dictionary has been built
    and the InputFormDescr has been extended with a 'widget' entry..

    The form can be constructed to be modal, i.e. the form grabs the focus
    and requires the form to be dissmissed using the OK or CANCEL button
    before the user can proceed. (this is the default behavior)
    When the form is NOT modal, it can still block and require the user to
    dismiss the form. The difference with modal is that other forms remain
    active, allowing for instance to pick something in a nother window.
    The program however blocks until the form is dismissed.
    When both the modal and blocking are 0, the form is displayed and
    worksasynchronously. The OK and Cancel buttons are not displayed.
    The form descriptor should provide a button to dismiss the form. The
    call back function for that button should call the destroy method. In
    this mode the current value of widgets have to retrieved manually.
    This can be done either by accessiing the widget directly:
    ipf[i]['widget'].get() or ipf.entryByName['mywidget'][widget].get();
    or by calling the checkValues() method that returns a dictionary.
    You can also by giving to the argument defaultDirection the values
    'col' or 'row', choose the direction to build the whole inputForm.
    For example if row is choosen all the widget if not specifies
    otherwise will be added row after row..
    When the flag closeWidthWindow is set to 1 then the form created is
    closed simultaneously with the main window. It is set to 1 by default.
    """
    def __init__(self, master, root, descr, modal=1, blocking=0,
                 defaultDirection='row', closeWithWindow=1,
                 okCfg={'text':'OK'}, cancelCfg={'text':'Cancel'},
                 initFunc=None, onDestroy=None, okcancel=1,
                 scrolledFrame=0, width=None, height=None, okOnEnter=True,
                 help=None):
        
        assert isinstance(descr, InputFormDescr)

        # master is the container widget
        self.master = master
        self.root = root
        self.ownsRoot = False
        
        # Creates a Tkinter Toplevel widget if self.root is None.

        if root is None:
            self.root = tkinter.Toplevel(self.master)

        if isinstance(self.root, tkinter.Toplevel):
            self.ownsRoot = True # use to decide to iconify or forget
            
        # Setting the title of self.root with the one provided.
        if isinstance(self.root, tkinter.Tk) or \
           isinstance(self.root, tkinter.Toplevel):
            self.root.title(descr.title)


        self.modal = modal  # when true, form will grab focus AND block
        self.blocking = blocking  # when true and form is not modal,
                                  # make the form block
                                  
        self.okcancel = okcancel  # When true an OK/CANCEL button will be created.

        # Variable used to grid the given widgets
        self.lastCol = 0
        self.lastRow = 0
        self.lastUsedRow = 0
        self.lastUsedCol = 0
        self.defaultDirection = defaultDirection

        # dictionary widgetID_entries : key = id of the widget
        #                        value = entries corresponding
        self.widgetID_entries = {}
        self.descr = descr
        self.descr.form = self

        # width and height
        autoSize = False
        if width is None:
            width = 400
            autoSize = True
        if height is None:
            height = 200
            autoSize = True
        
        # MAIN FRAME
        self.mf = tkinter.Frame(self.root, bd=4)

        # form main frame
        if scrolledFrame:
            if autoSize:
                self.sf = Pmw.ScrolledFrame(self.mf)
            else:
                self.sf = Pmw.ScrolledFrame(self.mf, usehullsize = 1,
                                            hull_width = width,
                                            hull_height = height)

            self.sf.pack(padx=5, pady=3, fill='both', expand=0)
            self.f = self.sf.interior()
            
        else:
            self.f = tkinter.Frame(self.mf, bd=4, relief='groove',
                                   width=width, height=height)


        
        # Create the parentName dictionnary with the name of the parent
        # and the container object
        self.parentName = {'mainContainer':{'container':self.f, 'widgets':[]}}

        #if not hasattr(self, 'last'): self.last = {}

        # Keep track of all the notebooks widget contained in the form and
        # all the Pmw widgets.
        self.notebooks = []
        self.wwidgets=[]

        wnum = 0
        for e in descr:
            # build reverse lookup dictionary
            if 'name' in e:
                descr.entryByName[e['name']] = e
            else:
                descr.entryByName[wnum] = e
            wnum=wnum+1
            self.addEntry(e)

        # Call the initFunc method if not none.
        if not initFunc is None:
            initFunc(self)

        # Pack the frames.
        self.f.pack(fill='both', expand=0)
        self.mf.pack(fill='both', expand=0)

        # Build the ok and cancel buttons with the information
        # provided.
        if (self.blocking or self.modal) and self.okcancel:
            self.buildOkCancelButtons( okCfg, cancelCfg, help=help )
            if self.ownsRoot:
                self._set_transient(self.master)
                if okOnEnter:
                    self.root.bind("<Return>", self.OK_cb)
        else:
            if help:
                self.buildHelpButton(help)
        # FIXME: Need to find all the possibilities....
        # bind the [x] button of the inputform 
        if self.ownsRoot:
            if onDestroy is None:
                if (self.blocking or self.modal) and self.okcancel:
                    # self Cancel_cb is modal or blocking
                        self.root.protocol('WM_DELETE_WINDOW', self.Cancel_cb)
                else:
                    # to just withdraw
                    if isinstance(self.root, tkinter.Tk) or \
                        isinstance(self.root, tkinter.Toplevel):
                        self.root.protocol('WM_DELETE_WINDOW', self.withdraw)
            else:
                # or finally the specified function.
                self.root.protocol('WM_DELETE_WINDOW', onDestroy)
            if closeWithWindow :
                if hasattr(self.root, "transient"):
                    self.root.transient(self.master)
            #self._set_transient(self.master)

        ## Set the notebook widget to their natural size.
        for nb in self.notebooks:
            nb.setnaturalsize()
        ## Call the alignLabels method to align all the labels of the
        ## different Pmw widgets added to the form
        # apply(Pmw.alignlabels, (self.wwidgets,), {})

        if autoSize:
            self.autoSize()
        
    def autoSize(self):
        if not self.ownsRoot: return
        self.root.update() # force packing to happen to size are right
        # self.root can be Tkinter.Frame instance. Frames do not have geometry attribute.
        if hasattr(self.root, "geometry"):
            w = self.root.winfo_reqwidth()
            h = self.root.winfo_reqheight()
            self.root.geometry('%dx%d'%(w,h))


    ## ###########################################################
    ## # Methods handling the creation of each widgets
        
    def addEntry(self, entry = {}):
        """
        optional argument:
        entry -- ({}) Dictionary describing the widget to be added.
                 See the InputFormDescr documentation for the accepted
                 format.
                 
        This method creates each given widgets, find their position
        in the form, grid them, sets their default values etc...
        """
        w = None
        if 'gridcfg' not in entry: entry['gridcfg'] = {}

        # dealing with the fact that the given widget belongs to
        # another container than the default mainContainer.
        if 'parent' in entry and \
           entry['parent'] in self.parentName:
            parent = self.parentName[entry['parent']]['container']
            
        else:
            parent = self.f
            entry['parent'] = 'mainContainer'

        # widget type
        if 'widgetType' in entry: wtype = entry['widgetType']
        else: entry['widgetType'] = wtype = tkinter.Entry

        # InputForm
        if wtype == 'InputForm':
            w = self.addInputForm(parent, entry, wtype)
        # Container
        elif wtype == 'Container':
            w = self.addContainer(parent, entry, wtype)
            
        # FunctionCheckButton or FunctionRadiobutton:
        elif wtype in ['FunctionCheckbutton', 'FunctionRadiobutton']:
            w = self.addFunctionButton(parent, entry, wtype)

        # Save Button
        elif wtype == 'SaveButton':
            w = self.addSaveButton( parent, entry)

        # Open Button
        elif wtype == 'OpenButton':
            w = self.addOpenButton( self.f, entry)

        # ScrolledText
        elif wtype == 'ScrolledText':
            w = self.addScrolledText(parent, entry)

        # ListChooser:
        elif wtype ==  'ListChooser':
            w = self.addListChooser( parent, entry )

        # Entry
        elif issubclass( wtype, tkinter.Entry):
            w = self.addTkinterEntry( parent, entry )

        # other Tkinter Widget:
        elif issubclass( wtype, tkinter.Widget):
            if 'listtext' in entry:
                w = self.addGroupWidget(parent, entry)
            else:
                w = self.addTkinterWidget( parent, entry )

        # Python Mega Widgets  :
        elif issubclass(wtype, Pmw.MegaWidget):
            w = self.addPMegaWidget(parent, entry)
            self.wwidgets.append(w)

        # Pmw MegaArchetype
        elif issubclass(wtype, Pmw.MegaArchetype):
            w = self.addPMegaWidget(parent, entry)
            self.wwidgets.append(w)

        else:
            # You can add any customized widget if you pass the class in wtype
            if type(wtype) in [type, type]:
                if 'wcfg' not in entry:
                    entry['wcfg'] = {}
                #w = apply(wtype, (self.f,),entry['wcfg'])
                w = wtype(*(parent,), **entry['wcfg'])
                if hasattr(w, 'set') and 'defaultValue' in entry:
                    w.set(entry['defaultValue'])
                #w = wtype(self.f,entry['wcfg'])
                if 'gridcfg' not in entry:
                    entry['gridcfg']={}
            else:
                raise ValueError("%s unsupported widget type"%str(wtype))

        if w :
            # Once the widgets has been created
            # If a parent has been specified for this widgets
            if entry['parent'] in self.parentName and \
               'name' in entry:
                widgets = self.parentName[entry['parent']]['widgets']
                widgets.append(entry['name'])

            # The widget is a container widget
            if 'container' in entry:
                for key, cont in list(entry['container'].items()):
                    if type(cont) is bytes:
                        exec("self.parentName['%s'] = {'container':%s, 'widgets':[]}"%(key,cont))
                        # print self.parentName[key]
                    #elif type(cont) is types.InstanceType:
                    elif isInstance(cont) is True:
                        self.parentName[key]= {'container':cont, 'widgets':[]}

            # If tooltip string has been specified creating a Balloon
            # widget to display the tooltip
            if 'tooltip' in entry:
                balloon = Pmw.Balloon(parent, yoffset=30)
                try:
                    balloon.bind(w, entry['tooltip'])
                    entry['balloon'] = balloon
                except:
                    entry['balloon'] = None

            # Finding the position of the widget in the form
            gridcfg = self.findGridPosition(entry['gridcfg'])
            if 'sticky' not in gridcfg:
                gridcfg['sticky']='wens'
            if 'weight' in gridcfg:
                weight=gridcfg['weight']
                del gridcfg['weight']
            else:
                weight=1

            # Grid the widget.
            if hasattr(w, 'grid'):
                w.grid(*(), **gridcfg)

            # the 2 following lines mean that the frame self.f
            # will allow the widget w placed at gridcfg['row'],
            # gridcfg['column'] to resize with the corresponding
            # weight with itself only if the sticky option of the
            # grid method is specified.
            self.f.rowconfigure(gridcfg['row'], weight = weight)
            self.f.columnconfigure(gridcfg['column'], weight = weight)
            # Same for container widget
            if not parent == self.f:
                parent.rowconfigure(gridcfg['row'], weight = weight)
                parent.columnconfigure(gridcfg['column'], weight = weight)

        # Keep track of the different notebooks
        if wtype in [Pmw.NoteBook,]:
            self.notebooks.append(w)

        # Handle some widgets a particular way
        if wtype not in ['ListChooser','ScrolledText', tkinter.Entry]:
            entry['widget'] = w
            self.widgetID_entries[w] = entry
        else:
            self.widgetID_entries[entry['widget']] = entry

    def addContainer(self, f, entry, wtype):
        """
        Method handling container widgets
        """
        if 'wcfg' in entry: wcfg = entry['wcfg']
        if 'frameType' in entry and entry['frameType']=='scrolledFrame':
            parent = Pmw.ScrolledFrame(*(f,), **wcfg)
            self.parentName[entry['name']] = {'container':parent.interior(),
                                              'widgets':[]}
        else:
            parent = tkinter.Frame(*(f,), **wcfg)
            self.parentName[entry['name']] = {'container':parent,
                                              'widgets':[]}
        return parent

    def addInputForm(self, f, entry, wtype):
        """
        Method handling InputForm widgets
        """
        if 'wcfg' in entry: wcfg = entry['wcfg']
        else: wcfg = {}
        frame = tkinter.Frame(f)
        description = entry['description']
        wcfg['modal'] = 0
        wcfg['blocking'] = 0
        frame.form = InputForm(*( frame, description), **wcfg)
        return frame
    

    def addFunctionButton(self, f, e, wtype):
        """Add a FunctionCheckbutton or FunctionRadiobutton"""

        if 'wcfg' in e: wcfg = e['wcfg']
        else: wcfg = {}
        # if a variable is bound to the widget 
        if 'variable' not in e:
            if wtype=='FunctionCheckbutton': f.var=tkinter.IntVar()
            elif wtype=='FunctionRadiobutton': f.var=tkinter.StringVar()

            e['variable'] = f.var
            wcfg['variable'] = e['variable']
        else:
            wcfg['variable'] = e['variable']
            if 'defaultValue' in e:
                wcfg['variable'].set(e['defaultValue'])

        if 'text' in e: wcfg['text'] = e['text']
        if 'value' in e: wcfg['value'] = e['value']

        if wtype=='FunctionCheckbutton':
            w = tkinter.Checkbutton(*(f,), **wcfg)
        elif wtype=='FunctionRadiobutton':
            w = tkinter.Radiobutton(*(f,), **wcfg)
        w.bind("<Button-1>", self.FuncButton_cb, "+")

        return w

    def addScrolledText(self, f, e):
        """Method handling ScrolledText widgets"""
        f1 = tkinter.Frame(f)
        if 'label' in e:
            tkinter.Label(f1, text=e['label']).pack(anchor = 'w')
        w = tkinter.Text(f1)
        w.vbar = tkinter.Scrollbar(f1, name='vbar')
        w.vbar['width'] = 15
        w.vbar.pack(side = tkinter.LEFT, fill = 'y')
        if 'size' in e:
            w['width'] = e['size'][0]
            w['height'] = e['size'][1]
        else:
            w['width'] = 80
            w['height'] = 2
        if 'defaultValue' in e:
            w.insert(tkinter.END, e['defaultValue'])
        if 'writeButton' in e:
            wb = self.addWriteButton(f1)
            # Bind the descr to the Button write, will be used
            # by the write_cb function to write in the appropriate Text
            # zone.
            self.widgetID_entries[wb] = e
            wb.pack(side = 'right')

        if 'readButton' in e:
            rb = self.addReadButton(f1)
            # Bind the descr to the Button read, will be used
            # by the read_cb function to read in the appropriate Text
            # zone.
            self.widgetID_entries[rb] = e
            rb.pack(side = 'right')
        w.pack()
        w['yscrollcommand'] = w.vbar.set
        w.vbar['command'] = w.yview
        e['widget'] = w
        return f1

    def addReadButton(self,frame):
        """ Method to add a button to read text from file and add it in a
        scrolled Text or a text zone. """
        readBut = tkinter.Button(frame, text = 'Read')
        readBut.bind("<Button-1>", self.read_cb, "+")
        return readBut

    def addWriteButton(self,frame):
        """ Method to add a button to write the text contained in a scrolled
        text or a text zone in a file."""
        
        writeBut = tkinter.Button(frame, text = 'Write')
        writeBut.bind("<Button-1>", self.write_cb, "+")
        return writeBut

    def addSaveButton(self, frame, wdescr):
        """ Method to add a button or a checkbutton or a radiobutton  to open
        a Save
        fileBrowser.
        """
##          wdescr['widgetType'] = Tkinter.Button
        if 'typeofwidget' in wdescr:
            wdescr['widgetType'] = wdescr['typeofwidget']
        else:
            wdescr['widgetType'] = tkinter.Button
        saveBut = self.addTkinterWidget(frame, wdescr)
        wdescr['widgetType'] = 'SaveButton'
        saveBut.bind("<Button-1>", self.save_cb, "+")
        return saveBut

    def addOpenButton(self, frame, wdescr):
        """ Method to add a button , checkbutton or radiobutton to open file browser.
        """
        wdescr['widgetType'] = tkinter.Button
        openBut = self.addTkinterWidget(frame, wdescr)
        wdescr['widgetType'] = 'OpenButton'
        openBut.bind("<Button-1>", self.open_cb, "+")
        return openBut

    def addTkinterWidget(self, f, e):
        """
        Method handling any Tkinter widgets
        """
        if 'wcfg' in e: wcfg = e['wcfg']
        else: wcfg = {}
        if 'text' in e: wcfg['text'] = e['text']

        if 'variable' in e:
            wcfg['variable'] = e['variable']

        elif 'variable' not in wcfg and \
             (issubclass(e['widgetType'],tkinter.Radiobutton) or \
             issubclass(e['widgetType'],tkinter.Checkbutton)):
            e['variable'] = wcfg['variable'] = tkinter.StringVar()

        if 'textvariable' in e:
            wcfg['textvariable'] = e['textvariable']
            e['variable'] = e['textvariable']
                    
        if 'command' in e: wcfg['command'] = e['command']

        if 'value' in e: wcfg['value'] = e['value']

        w = e['widgetType'](*(f,), **wcfg)

        if 'defaultValue' in e:
            if  hasattr(e['widgetType'], 'set') :
                w.set(e['defaultValue'])
            elif hasattr(e['widgetType'], 'insert') \
                 and ('textvariable' not in e or \
                      'textvariable' not in wcfg):
                #Text case but not Entry
                w.insert(tkinter.END, e['defaultValue'])

            elif 'textvariable' in wcfg:
                #Entry widget case.
                wcfg['textvariable'].set(e['defaultValue'])

            elif 'variable' in wcfg:
                wcfg['variable'].set(e['defaultValue'])
        return w

    def addPMegaWidget(self, f, e):
        """
        Method handling any Python Mega Widgets (Pmw)
'-' character are replaced by '-' to avoid Pmw exceptions, the oppposite
translation is done in checkValues()
        """
        # - pass directly the python mega widget dictionary as an argument
        if 'wcfg' in e:
            wcfg = e['wcfg']
        else:
            wcfg = {}
        w = e['widgetType'](*(f,), **wcfg)
        #w.pack(fill = 'x', padx = 10, pady = 10)
        # This mega widget can have components to be added to it
        # This will replace the listtext and listoption.
        if 'componentcfg' in e:
            for compcfg in e['componentcfg']:
                w.add(*(compcfg['name'].replace('_','-'),), **compcfg['cfg'])

        elif 'listtext' in e:
            if 'listoption' in e and \
               type(e['listoption']) is dict:
                for name in e['listtext']:
                    if name in e['listoption']:
                        w.add(*(name.replace('_','-'),), **e['listoption'][name])
            else:
                for name in e['listtext']:
                    w.add(name.replace('_','-'))

##         if e.has_key('container'):
##             for key, cont in e['container'].items():
##                 if type(cont) is types.StringType:
##                     exec("self.parentName['%s'] = {'container':%s, 'widgets':[]}"%(key,cont))
##                     #print self.parentName[key]
##                 elif type(cont) is types.InstanceType:
##                     self.parentName[key]= {'container':cont, 'widgets':[]}

        if 'defaultValue' in e:
            if isinstance(w, Pmw.ScrolledText):
                w.settext(e['defaultValue'])
            elif isinstance(w, Pmw.ScrolledListBox) \
                     or isinstance(w, Pmw.ScrolledText):
                w.setvalue(e['defaultValue'])
            elif hasattr(w,'set'):
                w.set(e['defaultValue'])
            elif hasattr(w,'invoke'):
                try:
                    w.invoke(e['defaultValue'])
                except:
                    if hasattr(w,'selectitem'):
                        try:
                            w.selectitem(e['defaultValue'])
                        except:
                            print('Could not use %s as the defaultValue!'%e['defaultValue'])
                    else:
                        print('sorry cannot set to the defaultValue', e['defaultValue'])
            elif hasattr(w,'insert'):
                w.insert('end',e['defaultValue'])
        return w

    def addGroupWidget(self,f,e):
        """
        Method handling GroupWidgets
        """
        if 'wcfg' not in e:
            e['wcfg'] = {}

        if 'gridcfg' not in e:e['gridcfg'] = {}
        
        j = 0
        if 'direction' in e:
            direction = e['direction']
        else:
            direction = self.defaultDirection
        if 'groupedBy' in e:
            groupedBy=e['groupedBy']
        else:
            groupedBy = len(e['listtext'])
        number=len(e['listtext'])/groupedBy

        if len(e['listtext'])%groupedBy != 0:
            number=number+1

        if j==0:
            # For the first radiobutton, use findGridPosition to get
            # the widget's position.
            gridcfg = self.findGridPosition(e['gridcfg'])
            firstrow = e['gridcfg']['row']
            firstcolumn = e['gridcfg']['column']
        for n in range(number):
            if j!= 0:
                if direction == 'col':
                    self.lastUsedRow = firstrow + n
                    self.lastUsedCol = firstcolumn - 1
                else:
                    self.lastUsedRow = firstrow - 1
                    self.lastUsedCol = firstcolumn + n
                
            while j<(groupedBy*(n+1)) and j<len(e['listtext']):
                e['wcfg']['text'] = e['listtext'][j]
                e['wcfg']['value'] = e['listtext'][j]
                #w = apply(e['widgetType'], (f,), wcfg)
                if j!=0:
                    if direction == 'col':
                        e['gridcfg']['column'] = self.lastUsedCol+1
                        e['gridcfg']['row']= self.lastUsedRow
                        gridcfg = self.findGridPosition(e['gridcfg'])

                    else:
                        e['gridcfg']['row'] = self.lastUsedRow+1
                        e['gridcfg']['column'] = self.lastUsedCol
                        self.lastUsedRow = self.lastUsedRow + 1
                        gridcfg = self.findGridPosition(e['gridcfg'])
                
                w = self.addTkinterWidget(f,e)
                if w:
                    if e['parent'] in self.parentName and \
                       'name' in e:
                        widgets = self.parentName[e['parent']]['widgets']
                        widgets.append(e['name'])
                    
                w.grid(*(), **gridcfg)
                j=j+1

        if 'defaultValue' in e:
            e['wcfg']['variable'].set(e['defaultValue'])
        else:
            e['wcfg']['variable'].set(e['listtext'][0])
                
                
        return None


    def addTkinterEntry(self, f, e ):
        """
        add an Entry widget with a possible label
        """
        if 'textvariable' not in e['wcfg']:
            e['wcfg']['textvariable'] = tkinter.StringVar()

	command = getkw(e['wcfg'], 'command')
	eventType = getkw(e['wcfg'], 'eventType')
	commandArgs = getkw(e['wcfg'], 'commandArgs')
        label = getkw(e['wcfg'],'label')
        if label:
            f1 = tkinter.Frame(f)
            tkinter.Label(f1, text=label).grid(row=0,
                                                    column=0,sticky='e')
            w = self.addTkinterWidget(f1, e)
            w.grid(row=0, column=1, sticky='w')
            e['widget'] = w
            self.widgetID_entries[w] = e
            returnVal =  f1

        else:
            w = self.addTkinterWidget(f, e)
            e['widget'] = w
            self.widgetID_entries[w] = e
            returnVal = w
	
        if command:
            if eventType is None: eventType = '<Return>'
            if commandArgs:
                w.bind(eventType, CallBackFunction(command, commandArgs))
            else:
                w.bind(eventType, command)

	return returnVal
	    
    def addListChooser(self, f, entry={}):
        mode = 'single'
        if 'mode' in entry: mode = entry['mode']
        title = 'Title'
        if 'title' in entry: title = entry['title']
        width = 40
        if 'width' in entry: width = entry['width']
        if 'defaultValue' in entry:
            defaultValue = entry['defaultValue']
        else:
            defaultValue = None

        if 'lbpackcfg' in entry: lbpackcfg=entry['lbpackcfg']
        else: lbpackcfg= {'fill':'x', 'expand':0}

        if 'lbwcfg' in entry: lbwcfg=entry['lbwcfg']
        else: lbwcfg = {}
        
        lc = entry['widget'] = ListChooser(f, mode, title=title,
                                           entries = entry['entries'],
                                           width=width,
                                           defaultValue = defaultValue,
                                           lbpackcfg = lbpackcfg,
                                           lbwcfg = lbwcfg)
        if 'command' in entry:
            lc.lb.bind("<ButtonRelease-1>", entry['command'], '+')
        return lc.top

        

    ## #############################################################
    ## CALLBACK FUNCTIONS
    
    def open_cb(self, event = None):
        from ViewerFramework.VFGUI import fileOpenAsk
        self.releaseFocus()
        entry = self.widgetID_entries[event.widget]

        if 'idir' in entry: idir = entry['idir']
        else: idir = None

        if 'ifile' in entry: ifile = entry['ifile']
        else: ifile = None

        if 'title' in entry: title = entry['title']
        else: title = ''

        if 'types' in entry: types = entry['types']
        else: types = [('All types', '*.*')]
        
        newfile = fileOpenAsk(self.root, idir = idir, ifile = ifile,
                              types = types,
                              title = title)
        if 'callback' in entry and newfile:
            entry['callback'](newfile)
        entry['filename'] = newfile
        self.grabFocus()

        
    def save_cb(self, event = None):
        from ViewerFramework.VFGUI import fileSaveAsk
        self.releaseFocus()
        entry = self.widgetID_entries[event.widget]

        if 'idir' in entry: idir = entry['idir']
        else: idir = None

        if 'ifile' in entry: ifile = entry['ifile']
        else: ifile = None

        if 'title' in entry: title = entry['title']
        else: title = ''

        if 'types' in entry: types = entry['types']
        else: types = [('All types', '*.*')]
        
        newfile = fileSaveAsk(self.root,idir=idir,ifile=ifile,types= types,
                              title = title)
        if 'callback' in entry and newfile:
            entry['callback'](newfile)
        entry['filename'] = newfile
        self.grabFocus()


    def read_cb(self, event = None):
        from ViewerFramework.VFGUI import fileOpenAsk
        self.releaseFocus()
        entry = self.widgetID_entries[event.widget]
        if 'readFileType' in entry:
            fileType = entry['readFileType']
        else:
            fileType = [ ('Python Files', '*.py')]

        if 'idir' in entry: idir = entry['idir']
        else: idir = None

        if 'ifile' in entry: ifile = entry['ifile']
        else: ifile = None

        newfile = fileOpenAsk(self.root,idir=idir,ifile=ifile,
                              types = fileType,
                              title = 'Read function from file: ')
        if newfile != None:
            f = open(newfile, 'r')
            lines = f.read()
            entry['widget'].insert(tkinter.END,lines)
            f.close()
        self.grabFocus()


    def write_cb(self, event = None):
        from ViewerFramework.VFGUI import fileSaveAsk
        self.releaseFocus()
        entry = self.widgetID_entries[event.widget]
        if 'writeFileType' in entry:
            fileType = entry['writeFileType']
        else:
            fileType = [ ('Python Files', '*.py')]

        if 'idir' in entry: idir = entry['idir']
        else: idir = None

        if 'ifile' in entry: ifile = entry['ifile']
        else: ifile = None

        newfile = fileSaveAsk(self.root,idir=idir,ifile=ifile,types = fileType,
                              title = 'Save function in file: ')
        if newfile != None:
            f = open(newfile, 'w')
            f.write(entry['widget'].get(1.0, tkinter.END))
            f.close()
        self.grabFocus()


    def FuncButton_cb(self,event = None):
        """
        method used to display a window allowing the user to type in
        Python code
        """
        import tkinter.messagebox

        # get the widget's description
        assert 'variable' in self.widgetID_entries[event.widget]
        widget = self.widgetID_entries[event.widget]

        # code typein widget is modal
        self.releaseFocus()
        if 'defaultText' not in widget: widget['defaultText']=''
        if widget['widgetType']=='FunctionRadiobutton':
            widget['variable'].set(widget['text'])
        v = widget['variable'].get()
        
        # create a form description
        ifd = InputFormDescr(title = 'Define a function')
        if 'readFileType' not in widget:
            widget['readFileType']= [('Python Files','*.py')]
        if 'writeFileType' not in widget:
            widget['writeFileType']= [('Python Files','*.py')]
        
        if v or v == 'function' or v == 'mapfunction':
            if 'size' in widget:

                ifd.append({'name': 'UserFunction',
                            'label':widget['label'],
                            'widgetType':'ScrolledText',
                            'readFileType':widget['readFileType'],
                            'writeFileType': widget['writeFileType'],
                            'readButton':1,'writeButton':1,
                            'size':widget['size'],
                            'defaultValue':widget['defaultText']})
            else:
                ifd.append({'name': 'UserFunction',
                            'label':widget['label'],
                            'widgetType':'ScrolledText',
                            'readFileType':widget['readFileType'],
                            'writeFileType': widget['writeFileType'],
                            'readButton':1,'writeButton':1,
                            'defaultValue':widget['defaultText']})
        else:
           ifd.append({'name': 'UserFunction',
                       'label':'type in zone' ,
                       'widgetType':'ScrolledText',
                       'readFileType':widget['readFileType'],
                       'writeFileType': widget['writeFileType'],
                       'readButton':1,'writeButton':1,
                       'defaultValue':widget['defaultText']})

        window = InputForm(self.master, self.root, ifd)
        vals = window.go()
        self.grabFocus()
        if widget['widgetType']=='FunctionCheckbutton':
            widget['variable'].set('0')
                                   
        if len(vals)==0: return
        widget['userFunction'] = vals['UserFunction']
        function = evalString(vals['UserFunction'])
        
        if widget['userFunction'] is None:
            tkinter.messagebox.showwarning('Color by properties',
                                     'no function definied')
            return

        if 'command' in widget:
            values = widget['command'](function)
            widget['result'] = values
        else:
            widget['result'] = function

    def PMWnameToReal(self, e, name):
        # PMW does not accept names containing '_' for widgets
        # addPMegaWidget systematically replaced '_' in names by '-'
        # this functions is called by checkValues to fo the opposite
        # i.e. replace '-' by '_' if needed
        names = []
        if 'componentcfg' in e:
            names = list(map( lambda x, e=e: e['componentcfg']['name'], e))
        elif 'listtext' in e:
            if 'listoption' in e and \
               type(e['listoption']) is dict:
                for name in e['listtext']:
                    if name in e['listoption']:
                        names.append(name)
            else:
                names = e['listtext']

        # len of names can be 0, for instance the widget listing the MSMS
        # surfaces has an empty list for names
        if len(names)==0 or name in names:
            return name
        else:
            return name.replace('-','_')

    ## #####################################################################
    ## 

    def checkValues(self, container='all'):
        """
        optional argument:
        container -- ('all') specifies the name of a container 
                     
        Method which gets the values of each widgets in the given container
        It will return a dictionary where the key is the name of the widget
        and the value is the value of the widget.
        This method is called by the OK_cb when the OK/CANCEL buttons
        exist and has to called otherwise.
        If the widget doesn't implement a get method a special case has to be
        implemented in this method.
        """
        values = {}
        if container == 'all':
            contNames = list(self.parentName.keys())

        elif type(container) is bytes:
            contNames = [container,]
        else:
            contNames = container
        
        pN = self.parentName
        
        entryBN = self.descr.entryByName
        wNames = []
        for name in contNames:
            wNames = wNames + pN[name]['widgets']
        if len(wNames) == 0:
            print("WARNING: No widgets in these containers")
            
        for wName in wNames:
            wtype = entryBN[wName]['widgetType']
            e = entryBN[wName]
            if wtype == 'Container':
                continue

            elif wtype == 'InputForm':
                val = e['widget'].form.checkValues()
                
            elif wtype in['SaveButton','OpenButton']:
                if 'filename' in e:
                    val = e['filename']
                else:
                    val = None
            elif wtype in ['FunctionCheckbutton', 'FunctionRadiobutton']:
                if 'result' in e:
                    val = (e['text'], e['result'], e['userFunction'])
            # ScrolledText
            elif wtype in ['ScrolledText', tkinter.Text]:
                val = e['widget'].get(1.0, tkinter.END)

            # ListChooser:
            elif wtype ==  'ListChooser':
                val = e['widget'].get()

            # Tkinter Widget:
            elif issubclass( wtype, tkinter.Widget):
                if hasattr(e['widget'], 'get'):
                    val = e['widget'].get()
                elif 'wcfg' in e and 'variable' in e['wcfg']:
                    val = e['wcfg']['variable'].get()
                elif 'variable' in e:
                    # Here right now because all the descr don't have a wcfg.
                    val = e['variable'].get()
                else:
                    # Add this else for the widget with a
                    # name but no variable.
                    val = None

            elif issubclass(wtype, Pmw.ComboBox):
                getVal  = e['widget'].get()
                selVal = e['widget'].getcurselection()
                if not getVal in selVal:
                    val = (getVal,)
                else:
                    val = selVal
                        
            elif issubclass(wtype, Pmw.MegaWidget):
                if (hasattr(e['widget'],'get') and \
                   hasattr(e['widget'],'getcurselection')) or \
                   hasattr(e['widget'],'getcurselection'):
                    res=e['widget'].getcurselection()
                    val = []
                    if type(res) != bytes:
                        for name in res:
                            val.append( self.PMWnameToReal(e, name) )
                    else:
                        val = self.PMWnameToReal(e, res)
                       
                elif hasattr(e['widget'],'get'):
                    val = e['widget'].get()
                    val = self.PMWnameToReal(e, e['widget'].get())
                elif hasattr(e['widget'],'getint'):
                    val = e['widget'].getint()
                elif hasattr(e['widget'],'getstring'):
                    val = self.PMWnameToReal(e, e['widget'].getstring())
                else:
                    val=None
                
            elif hasattr(e['widget'],'get'):
                val = e['widget'].get()
            else:
                val = None
                    
            if val==None:
                continue

            if 'required' in e:
                if len(val) == 0:
                    self.releaseFocus()
                    t = 'missing value for field %s' % e['name']
                    
                    SimpleDialog(self.root, text=t,
                                 buttons=["Continue"],
                                 default=0,
                                 title="Missing required field").go()
                    self.grabFocus()
                    return 'Error'

            if 'type' in e:
                try:
                    val = e['type'](val)
                except:
                    self.releaseFocus()
                    t = 'invalid value for field %s' % e['name']
                    SimpleDialog(self.root, text=t,
                                 buttons=["Continue"],
                                 default=0,
                                 title="Invalid value").go()
                    self.grabFocus()
                    return 'Error'
            ok = 1
            if 'validateFunc' in e:
                if 'validateArgs' in e:
                    arguments = e['validateArgs']
                else:
                    arguments  = ()
                if 'err_msg' in e:
                    err_msg = e['err_msg']
                else:
                    err_msg  = 'invalid value for field %s' % e['name']
                    
                v =  e['validateFunc'](*(val,) + arguments)
                
                if not v: ok = 0

            elif 'validValues' in e:
                if not val in e['validValues']: ok = 0


            if ok:
                values[e['name']] = val
            else:
                self.releaseFocus()
                t=err_msg
                SimpleDialog(self.root, text=t, buttons=["Continue"],
                             default=0, title="Invalid value").go()
                self.grabFocus()
                return None
        return values


    def testForm(self, setWidget={}, container='all'):
        """This method can be called by tests instead of the go method, It will
        first set the widgets to the given values provided by the setWidgets dictionary
        then the checkValues method will be called and the val dictionary returned.
        """
        # Need to set the widgets to the given values if setWidget is
        # specified
        if setWidget:
            print("SETTING THE WIDGETS VALUE....")
            for name, val in list(setWidget.items()):
                if name in self.descr.entryByName:
                    widget = self.descr.entryByName[name]['widget']
                    if hasattr(widget, 'set'):
                        print('setting')
                        widget.set(val)
                    elif hasattr(widget, 'setvalue'):
                        print('setvalue')
                        widget.setvalue(val)
                    elif hasattr(widget, 'selectitem'):
                        print('selectitem')
                        widget.selectitem(val)
                    elif hasattr(widget, 'invoke'):
                        try:
                            print('invoking with val')
                            widget.invoke(val)
                        except:
                            print('invoking without')
                            widget.invoke()
                    else:
                        print("Could not set the value of ", name, val)
            print("=================================================")
        from time import sleep
        self.master.update()
        sleep(0.5)
        if hasattr(self, 'ok'):
            self.ok.invoke()
        else:
            self.values = self.checkValues(container)
        return self.values

##  Creates problems on MacOSX, th form gets lifted but not the combobox
##  pull down menus (see bindGeometry command)
##  Also is anoying with some other forms
##     def autoLift(self, event=None):
##         self.lift()
##         self.autoRaiseID = self.root.after(1000, self.autoLift)


    def go(self):
        """This method starts the inputform in modal mode"""
        # FIXME WHEN BLOCKING NOT NECESSARILY MODAL.
        if (self.modal or self.blocking) and self.okcancel == 1:
##            self.autoRaiseID = self.root.after(1000, self.autoLift)
            self.grabFocus()               # grabs the focus because modal
            self.root.mainloop()           # wait for the OK      
            self.releaseFocus()
##             if self.autoRaiseID is not None:
##                 self.root.after_cancel(self.autoRaiseID)
##                 self.autoRaiseID = None
                
            return self.values
        
    ## #########################################################3
    ## HELPER METHODS.
    def buildHelpButton(self, help, frame=None):
        """ Method which builds the HELP button.
        required argument:
        help == String providing the URL of the help. The callback function
                bound to the help button will open a webbrowser at this adress.
        optional argument
        frame -- (None) Specifies in which frame this button widget should be
                 added, If None provided then a new frame will be created.
        """
        if not isinstance(frame, tkinter.Frame):
            frame = tkinter.Frame(self.root, bd=4, relief='groove')
            frame.pack(side='bottom', fill='x', expand=0)
        ICONPATH = findFilePath('Icons', 'ViewerFramework')
        iconfile = os.path.join(ICONPATH,'32x32','info.gif')
        self.Icon = tkinter.PhotoImage(file=iconfile, master=self.master)
        self.help = tkinter.Button(*(frame,), **{'text':'?','image':self.Icon,'height':22,'width':24,
                           'command':CallBackFunction(self.showHelp_cb,
                                                      help=help)})
        self.help.grid(sticky='we', row=0, column=2)
        self.balloon = Pmw.Balloon(frame)
        self.balloon.bind(self.help, "Click to open help link:\n"+help)
        frame.rowconfigure(0, weight=1)
        frame.columnconfigure(0, weight=1)

    def showHelp_cb(self, help=None):
        """
        Callback function of the HELP button.
        It will open a webbrowser at the given help URL when specified.
        If it can't open the page then a warning message will be displayed.
        """
        import urllib.parse, webbrowser
        if help is None: return
        try:
            webbrowser.open(help)
        except:
            import warning
            warnings.warn("Error opening %s"%help)
        
        
    def buildOkCancelButtons(self, okCfg, cancelCfg, help=None):
        """
        This method creates an OK button and a CANCEL button.
        required arguments:
        okCfg == configuration dictionary of the OK button
        cancelCfg == configuration dictionay of the CANCEL button

        optional arguments:
        help == (None) if a help string is provided the buildHelpButton will be called
                with the f1 frame so that the 3 buttons are in the same frame.
        """
        from types import DictType
        # Frame for the buttons.
        #f1 = Tkinter.Frame(self.root, bd=4, relief='groove')
        f1 = tkinter.Frame(self.mf, bd=1, relief='groove')
        f1.pack(side='bottom', fill ='x', expand=0)
        self.okCancelButtonsMaster = f1
        # Ok button
        okopts = {'text':'OK', 'command':self.OK_cb}
        okFunc = None
        cancelFunc = None
        if type(okCfg) is DictType and not okCfg == {}:
            if 'command' in okCfg:
                okFunc = okCfg['command']
                del okCfg['command']
                okopts['command'] = CallBackFunction(self.OK_cb,
                                                     func=okFunc)
            okopts.update(okCfg)
        self.ok = tkinter.Button(*(f1,), **okopts)
        self.ok.grid(row = 0, sticky='we')

        # Cancel Button
        cancelopts = {'text':'Cancel', 'command':self.Cancel_cb}
        if type(cancelCfg) is DictType and not cancelCfg == {}:
            if 'command' in cancelCfg:
                cancelFunc = cancelCfg['command']
                del cancelCfg['command']
                cancelCfg['command'] = CallBackFunction(self.Cancel_cb,
                                                        func=cancelFunc)
            cancelopts.update(cancelCfg)
        self.cancel = tkinter.Button(*(f1,), **cancelopts)

        self.cancel.grid(row = 0, column = 1, sticky='we')
        if help:
            self.buildHelpButton(help, frame=f1)
        f1.rowconfigure(0, weight = 1)
        f1.columnconfigure(0, weight = 1)
        f1.columnconfigure(1, weight = 1)

    def OK_cb(self, event=None, func=None):
        """call back for OK button"""
        self.values = self.checkValues()
        if self.values=='Error':
            return
        self.root.quit()
        self.withdraw()
        # If a func is specified the func is then called
        if not func is None:
            func()

    def Cancel_cb(self, event = None, func=None):
        """call back for Cancel button"""
        #print 'ZZZZZZZ', event.widget, event.widget==self.root
        #print dir(event)
        #print event.type
        #print event.state, event.state=='??'
        #print event.width, event.height
        #print event.x, event.y
        #if event.state=='??': # focusOut from Thumbwheel not from frame
        #    return
        self.values = {}
        # If a func is specified the func is then called
        self.root.quit()
        self.withdraw()
        if not func is None:
            func()


    def findGridPosition(self, gridcfg):
        """find next available position in grid"""
        if 'col' in gridcfg and 'column' not in gridcfg:
            gridcfg['column'] = gridcfg['col']
            
        if 'column' in gridcfg and 'row' in gridcfg:
            requestedRow = gridcfg['row']
            if requestedRow < 0: requestedRow=self.lastRow+requestedRow
            if requestedRow < 0 or requestedRow >= self.lastRow:
                requestedRow = self.lastRow
                self.lastRow = self.lastRow+1              
            requestedCol = gridcfg['column']
            if requestedCol < 0: requestedCol=self.lastCol + requestedCol
            if requestedCol < 0 or requestedCol >= self.lastCol:
                requestedCol = self.lastCol
                self.lastCol = self.lastCol+1
 
        elif 'row' in gridcfg:
            requestedRow = gridcfg['row']
            self.lastUsedCol = self.lastUsedCol + 1
            requestedCol = self.lastUsedCol
            if requestedRow < 0: requestedRow=self.lastRow + requestedRow
            if requestedRow < 0 or requestedRow >= self.lastRow:
                requestedRow = self.lastRow
                self.lastRow = self.lastRow+1

        elif 'column' in gridcfg:
            requestedCol = gridcfg['column']
            self.lastUsedRow = self.lastUsedRow + 1
            requestedRow = self.lastUsedRow
            if requestedCol < 0: requestedCol = self.lastCol + requestedCol
            if requestedCol < 0 or requestedCol >= self.lastCol:
                requestedCol = self.lastCol
                self.lastCol = self.lastCol+1
                
        else:
            if self.defaultDirection == 'col':
                #requestedRow = self.lastUsedRow
                requestedRow = 0
                #self.lastRow = self.lastRow+1
		#NB: lastCol is always 1 more than biggest used
                requestedCol = self.lastCol
                self.lastCol = self.lastCol+1
            else:
                #requestedCol = self.lastUsedCol
                requestedCol =  0
                #self.lastCol = self.lastCol+1
		#NB: lastRow is always 1 more than biggest used
                requestedRow = self.lastRow
                self.lastRow = self.lastRow+1


        gridcfg['row'] = self.lastUsedRow = requestedRow
        gridcfg['column'] = self.lastUsedCol = requestedCol

        # make sure we point to first empty column and first empty line
        if self.lastRow <= requestedRow:
            self.lastRow = requestedRow + 1
        if self.lastCol <= requestedCol:
            self.lastCol = requestedCol + 1
            
	if 'columnspan' in gridcfg:
	    extraCol = gridcfg['columnspan']-1
	    self.lastCol = self.lastCol + extraCol
	    self.lastUsedCol = self.lastUsedCol + extraCol
	if 'rowspan' in gridcfg:
	    extraRow = gridcfg['rowspan']-1
	    self.lastRow = self.lastRow + extraRow
	    self.lastUsedRow = self.lastUsedRow + extraRow

        return gridcfg



    ###################################################################
    ### Tkinter methods.

    def grabFocus(self):
        """ Method to grab the focus"""
        if self.modal:
            self.root.grab_set()

    def releaseFocus(self):
        """ Method to release the focus of the form"""
        if self.modal:
            self.root.grab_release()

    def destroy(self):
        """ Method to destroy the forms"""
        if self.ownsRoot:
            self.root.destroy()
            #print 'DESTROYING FORM ROOT'
        else:
            #print 'DESTROYING FORM FRAME'
            self.mf.destroy()

            
    def lift(self, event=None):
        """
        Method to lift the inputform
        """
	self.root.winfo_toplevel().lift()

    def withdraw(self, event=None):
        """
        Method to withdraw the inputform
        """
        if self.ownsRoot:
            self.root.winfo_toplevel().withdraw()
        else:
            self.root.forget()
            self.mf.forget()

    def deiconify(self, event=None):
        """
        Method to deiconify the inputform.
        """
        if self.ownsRoot:
            self.root.winfo_toplevel().deiconify()
        else:
            #print 'DEICONIFY'
            self.root.pack(expand=0, fill='both')
            self.mf.pack(fill='both', expand=0)
            #self.root.bind('<Leave>', self.Cancel_cb)
        

def evalString(str):
    if len(str)==0:
        return
    try:
        function = eval("%s"%str)
    except:
        #try:
            obj = compile(str, '<string>', 'exec')
            exec(obj)
            function = eval(obj.co_names[0])
        #except:
        #    raise ValueError
    return function


from tkinter.filedialog import FileDialog
import os

class ModuleDialog(FileDialog):

    title = "Module Selection Dialog"
    def ok_command(self):
        file = self.get_selection()
        if not os.path.isfile(file):
            self.master.bell()
        else:
            self.quit(file)

class CommandDialog(FileDialog):

    title = "Command Selection Dialog"

    def ok_command(self):
        file = self.get_selection()
        if not os.path.isfile(file):
            self.master.bell()
        else:
            self.quit(file)



if __name__=='__main__':
    import os
    root = tkinter.Tk()
##      lb = ListChooser(root, title='Choose a letter',
##                       entries=[ ('a','comment a'),
##                                 ('b','comment b'),
##                                 ('c','comment c'),
##                                 ('d','comment a')] )
##      lb.insert(1, 'Hello', 'value inserted')
##      lb.select(1)
##      result = lb.get()
##      print result
##      res = lb.go()
##      print res

    # should allow to pass a function to check validity
    f = InputForm(root, title='input 2 strings a float and an int',
                  entries = [ {'label':'string 1'},
                              {'label':'string 2', 'defaultValue':'Hello'},
                              {'label':'float', 'type':float, 'required':1},
                              {'label':'int', 'type':int, 'defaultValue':7},
                              {'label':'check', 'widgetType':'Checkbutton'} ] )

    val = f.go()
    md = ModuleDialog(root)
    modulefile = md.go(key='modulekey', pattern = '*Module.py')
    if modulefile:
        print("important string is:", os.path.split(modulefile)[-1])
    cd = CommandDialog(root)
    commandfile = cd.go(key='commandkey', pattern = '*Commands.py')
    if commandfile:
        print("important string is:", os.path.split(commandfile)[-1])
                                                                               
