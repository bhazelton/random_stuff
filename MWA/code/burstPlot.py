import numpy as np
import matplotlib as mpl
mpl.use('WXAgg')
import matplotlib.backends.backend_wxagg as wxmpl
import logging
import wx
import wx.gizmos
import os

MODE_MULTIPLOT=0
MODE_SINGLEPLOT=1

MULTI_NROWS=4
MULTI_NCOLS=2

ID_OPEN=wx.NewId()
ID_SAVE=wx.NewId()
ID_EXIT=wx.NewId()


class MyFormatter(mpl.ticker.ScalarFormatter):
    def __call__(self, x, pos=None):
        N = len(self.locs)
        if pos==(N-1): return ''        # turn off first
        #elif pos==(0): return ''  # turn off last
        else: return int(x)


def readBurst(filename):
    data=np.fromfile(file=filename,dtype=np.float32)
    if len(data) != 4096:
        logging.error("File is not a valid burst average: %s"%filename)
        return None
    out=np.reshape(10*np.log10(data),(16,256))
    return out


class BurstFrame(wx.Frame):
    title="MWA Burst-mode data"
    def __init__(self,parent,id,title):
        wx.Frame.__init__(self,parent,id,title)

        self.create_menu()
        self.create_main_panel()
        self.statusbar=self.CreateStatusBar()

        self.vbox.Fit(self)
        

        self.mode=MODE_MULTIPLOT
        self.singleplot_bf=None

        #self.filelist=['All_1_Rx1_20090806.06920_avg','All_4_Rx1_20090806.06541_avg','NotDipole14_Rx1_20090806.41207_avg','burst.test']
        #self.file_select.Set(self.filelist)
        self.file_select.Clear()   
        self.curr_file=None

        self.do_button_state()

        self.data=None
        self.figtitle=None
        self.get_data(infile=self.curr_file)


        self.do_plot()


    def do_plot(self):
        if self.mode==MODE_MULTIPLOT:
            self.do_multiplot()    
        else:
            self.do_singleplot()
        self.autoscale_all()

    def do_multiplot(self):

        self.fig.clf()

        sax=None
        self.axes={}
        for bf in range(1,9):
            if sax is None:
                self.axes[bf]=self.fig.add_subplot(MULTI_NROWS,MULTI_NCOLS,bf)
                sax=self.axes[bf]
            else:
                self.axes[bf]=self.fig.add_subplot(MULTI_NROWS,MULTI_NCOLS,bf,sharex=sax,sharey=sax)
            
            #self.axes[bf].xaxis.set_visible(False)
            if bf not in (7,8):
                for label in self.axes[bf].xaxis.get_ticklabels():
                    label.set_visible(False)
            else:
                self.axes[bf].xaxis.set_label_text("Freq (MHz)",size=8)
                for label in self.axes[bf].xaxis.get_ticklabels():
                    label.set_fontsize(8)

            if bf%2 == 0:
                for label in self.axes[bf].yaxis.get_ticklabels():
                    label.set_visible(False)
            else:
                self.axes[bf].yaxis.set_label_text("Power (dB)",size=8)
                for label in self.axes[bf].yaxis.get_ticklabels():
                    label.set_fontsize(8)


            self.axes[bf].xaxis.set_major_formatter(MyFormatter())

            self.axes[bf].yaxis.set_major_formatter(MyFormatter())
            self.axes[bf].set_autoscale_on(False)
            
            
        #self.mpl_toolbar.Hide()
        #self.mpl_toolbar.Show()
        


        
        self.fig.subplots_adjust(left=0.05,right=0.95,wspace=0.0,hspace=0.0)
        self.vbox.Layout()
        self.draw_multi_spec()

    def do_singleplot(self):
        self.fig.clf()
        self.axes={}
        self.axes[self.singleplot_bf]=self.fig.add_subplot(1,1,1)
        #self.mpl_toolbar.Show()
        self.vbox.Layout()

        self.draw_single_spec()
        self.autoscale_all()


    def draw_single_spec(self):
        for bf,ax in self.axes.iteritems():
            del ax.lines[:]
        self.draw_bf(self.axes[self.singleplot_bf],self.singleplot_bf)


    def create_main_panel(self):
        self.mainpanel=wx.Panel(self,size=(640,480))
        self.dpi=100
        self.fig=mpl.figure.Figure((6.0,3.0),dpi=self.dpi)
        self.canvas=wxmpl.FigureCanvasWxAgg(self.mainpanel,wx.ID_ANY,self.fig)

        self.file_select=wx.ListBox(self.mainpanel,wx.ID_ANY)
        self.file_select.Bind(wx.EVT_LISTBOX,self.on_file_select_focus)
        self.file_select.Bind(wx.EVT_LISTBOX,self.on_file_select_focus)
        

        self.refresh_button=wx.ToggleButton(self.mainpanel,wx.ID_ANY,"Auto-Refresh")

        self.autoscale_button=wx.Button(self.mainpanel,wx.ID_ANY,"Autoscale")

        self.next_file_button=wx.Button(self.mainpanel,wx.ID_ANY,"Next Spec.")
        self.prev_file_button=wx.Button(self.mainpanel,wx.ID_ANY,"Prev Spec.")

        self.Bind(wx.EVT_TOGGLEBUTTON,self.on_refresh_button,self.refresh_button)
        self.Bind(wx.EVT_BUTTON,self.on_autoscale_button,self.autoscale_button)
        self.Bind(wx.EVT_BUTTON,self.on_next_file,self.next_file_button)
        self.Bind(wx.EVT_BUTTON,self.on_prev_file,self.prev_file_button)
        
        self.canvas.Bind(wx.EVT_LEFT_DCLICK,self.OnDClick)

        #self.canvas.mpl_connect('button_press_event',self.OnDClick)
        


        toolbar=wx.ToolBar(self.mainpanel,wx.TB_HORIZONTAL | wx.EXPAND)
        toolbar.SetToolBitmapSize((24,24))
        toolbar.AddTool(ID_OPEN,wx.ArtProvider.GetBitmap(wx.ART_FILE_OPEN,wx.ART_TOOLBAR,(24,24)))
        toolbar.AddTool(ID_SAVE,wx.ArtProvider.GetBitmap(wx.ART_FILE_SAVE,wx.ART_TOOLBAR,(24,24)))
        toolbar.Realize()
        
        self.mpl_toolbar=wxmpl.NavigationToolbar2WxAgg(self.canvas)
        

        self.vbox=wx.BoxSizer(wx.VERTICAL)
        self.vbox.Add(toolbar,0,wx.EXPAND,0)
        self.vbox.Add(self.canvas,1,wx.LEFT | wx.TOP | wx.GROW)
        self.vbox.Add(self.mpl_toolbar,flag=wx.ALL| wx.TB_HORIZONTAL)

        self.vbox.Add(self.file_select,0,wx.ALL | wx.EXPAND)

        self.hbox=wx.BoxSizer(wx.HORIZONTAL)
        self.hbox.AddSpacer(10)
        self.hbox.Add(self.refresh_button)
        self.hbox.AddSpacer(20)
        self.hbox.Add(self.autoscale_button)
        self.hbox.AddSpacer(20)
        self.hbox.Add(self.prev_file_button)
        self.hbox.Add(self.next_file_button)
        
        self.vbox.Add(self.hbox, 0, flag = wx.ALIGN_LEFT | wx.TOP)
        self.mainpanel.SetSizer(self.vbox)
        
        
        self.vbox.Fit(self)
        self.vbox.Layout()
        
        #self.mpl_toolbar.Hide()
    
        self.refresh_timer=wx.Timer(self)
        self.Bind(wx.EVT_TIMER,self.on_refresh,self.refresh_timer)
        self.refresh_button.SetValue(True)
        self.refresh_timer.Start(1000)



    def OnDClick(self, e):
        x,y=e.GetPositionTuple()


        winsize=self.fig.get_window_extent().get_points()
        y=winsize[1][1]-y

        inbf=None
        for bf,ax in self.axes.iteritems():
            r=ax.get_window_extent().get_points()
            if x>=r[0][0] and x<=r[1][0] and y>=r[0][1] and y<=r[1][1]:
                inbf=bf
                break
        
        if inbf is not None:
            if self.mode==MODE_MULTIPLOT:
                self.mode=MODE_SINGLEPLOT
                self.singleplot_bf=inbf
            elif self.mode==MODE_SINGLEPLOT:
                self.mode=MODE_MULTIPLOT
            self.do_plot()


    def OnExit(self,e):
        self.Close(True)
        

    def OnPick(self,e):
        pass


    def on_autoscale_button(self,e):
        self.autoscale_all()

    def autoscale_all(self):

        for bf,ax in self.axes.iteritems():
            ax.set_autoscale_on(True)
        self.draw_spec()
        for bf,ax in self.axes.iteritems():
            ax.set_autoscale_on(False)

        
    def get_data(self,infile=None):
        try:
            self.data=readBurst(infile)
        except:
            self.data=None

    def draw_spec(self):
        if self.mode==MODE_MULTIPLOT:
            self.draw_multi_spec()
        else:
            self.draw_single_spec()
        if self.figtitle is not None:
            self.figtitle.set_visible(False)
            del self.figtitle
        self.figtitle=self.fig.suptitle(self.curr_file)
        self.canvas.draw()

    def draw_bf(self,ax,bf):
        freqs=1.28*np.array(range(256))
        p1,=ax.plot(freqs,self.data[2*(bf-1),:],label="X",color='green')
        p2,=ax.plot(freqs,self.data[2*(bf-1)+1,:],label="Y",color='blue')
        ax.yaxis.set_major_locator(mpl.ticker.MaxNLocator(4))
        ax.xaxis.set_major_locator(mpl.ticker.MaxNLocator(6))
        ax.text(0.05,0.95,"Slot %i"%bf,linespacing=2,va="top",transform=ax.transAxes)

    def draw_multi_spec(self):

        if self.data is None:
            return

        
        for bf,ax in self.axes.iteritems():
            del ax.lines[:]
            self.draw_bf(ax,bf)
            #plot.set_xlabel("Freq (MHz)")
            #plot.set_ylabel("Power")
            
        
        l1=ax.lines[0]
        l2=ax.lines[1]
        self.fig.legend((l1,l2),("X","Y"),loc=8,ncol=2,mode="expand",borderaxespad=0,bbox_to_anchor=(0,0,1,0.1))
            

    def do_button_state(self):

        nstate=True
        pstate=True
        currnum=self.file_select.GetSelection()
        nitems=self.file_select.GetCount()
        if currnum>=nitems-1:
            nstate=False
        if currnum<=0:
            pstate=False
        
        if nstate:
            self.next_file_button.Enable()
        else:
            self.next_file_button.Disable()
        if pstate:
            self.prev_file_button.Enable()
        else:
            self.prev_file_button.Disable()

    def on_next_file(self,event):
        currnum=self.file_select.GetSelection()
        nitems=self.file_select.GetCount()
        if currnum <= nitems-2:
            self.file_select.Select(currnum+1)
            self.curr_file=self.file_select.GetString(currnum+1)

        self.do_button_state()
        self.get_data(self.curr_file)
        self.draw_spec()

    def on_prev_file(self,event):

        currnum=self.file_select.GetSelection()
        if currnum >0:
            self.file_select.Select(currnum-1)
            self.curr_file=self.file_select.GetString(currnum-1)

        self.do_button_state()
        self.get_data(self.curr_file)
        self.draw_spec()
  
    def on_refresh(self,event):

        old_data=self.data
        self.get_data(self.curr_file)
        if old_data is None or self.data is None:
            self.draw_spec()
            self.autoscale_all()
        elif (self.data != old_data).all():
            self.draw_spec()

    def on_refresh_button(self,e):
        time=self.refresh_timer
        if self.refresh_button.GetValue() is False:
            time.Stop()
        else:
            time.Start()
            

    def OnOpen(self,e):
        file_choices="Burst Average (*_avg)|*_avg"
        dlg = wx.FileDialog(
            self, 
            message="Select averages to open",
            #defaultDir=os.getcwd(),
            defaultFile='',
            wildcard=file_choices,
            style=wx.FD_OPEN | wx.FD_MULTIPLE | wx.FD_CHANGE_DIR)

        if dlg.ShowModal() == wx.ID_OK:
            self.file_select.Clear()
            self.file_select.SetItems(dlg.GetPaths())
            self.file_select.Select(0)
            self.do_button_state()
            self.curr_file=self.file_select.GetString(0)
            self.get_data(self.curr_file)
            self.draw_spec()
            self.autoscale_all()

    def OnSave(self,e):
        file_choices = "PNG (*.png)|*.png|EPS (*.eps)|*.eps"
        
        dlg = wx.FileDialog(
            self, 
            message="Save plot as...",
            defaultDir=os.getcwd(),
            defaultFile="plot.png",
            wildcard=file_choices,
            style=wx.SAVE | wx.FD_OVERWRITE_PROMPT)
        
        if dlg.ShowModal() == wx.ID_OK:
            path = dlg.GetPath()
            self.canvas.print_figure(path, dpi=self.dpi)
            self.flash_status_message("Saved to %s" % path)
    
    def flash_status_message(self, msg, flash_len_ms=1500):
        self.statusbar.SetStatusText(msg)
        self.timeroff = wx.Timer(self)
        self.Bind(
            wx.EVT_TIMER, 
            self.on_flash_status_off, 
            self.timeroff)
        self.timeroff.Start(flash_len_ms, oneShot=True)

    def on_flash_status_off(self, event):
        self.statusbar.SetStatusText('')


    def on_file_select_focus(self,event):
        do_rescale=False
        if self.data is None:
            do_rescale=True
        item=self.file_select.GetSelection()
        if item >=0:
            self.curr_file=self.file_select.GetString(item)
        self.get_data(self.curr_file)
        self.draw_spec()
        if do_rescale:
            self.autoscale_all()

    def create_menu(self):
        # Create the menubar
        file_menu=wx.Menu()
        file_menu.Append(ID_OPEN,"&Open\tCtrl-o","Open a document")
        file_menu.Append(ID_SAVE,"&Save to file\tCtrl-s","Save plot to file")
        file_menu.AppendSeparator()
        file_menu.Append(ID_EXIT,"E&xit","Exit the plotter")
        
        self.menuBar=wx.MenuBar()
        self.menuBar.Append(file_menu,"&File")
        self.SetMenuBar(self.menuBar)

        wx.EVT_MENU(self, ID_EXIT, self.OnExit)
        wx.EVT_MENU(self, ID_OPEN, self.OnOpen)
        wx.EVT_MENU(self, ID_SAVE, self.OnSave)





if __name__=='__main__':
    app=wx.PySimpleApp()
    app.frame=BurstFrame(None,-1,"Burst Plotter")
    app.frame.Show()
    app.MainLoop()
