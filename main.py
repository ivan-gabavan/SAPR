import wx
import wx.grid
import numpy as np
from numpy.linalg import linalg
from wx.lib.plot import PlotCanvas, PlotGraphics, PolyLine, PolyMarker

arr_rods = []
delta = []
sigma = 0.0
E = 0.0
leftSup = False
rightSup = False

class rod:
    def __init__ (self,L,A,q):
        self.L = L
        self.A = A
        self.q = q
    L = 0
    A = 0
    q = 0
class Proc:
    Rods = []
    Forces = []
    def __init__(self):
        global arr_rods
        for i in range(len(arr_rods)//4):
            self.Forces.append(float(arr_rods[i*4]))
            self.Rods.append(rod(float(arr_rods[i*4+1]), float(arr_rods[i*4+2]), float(arr_rods[i*4+3])))
        self.Forces.append(float(arr_rods[len(arr_rods)-1]))
    def ACalculate(self):
        global E,leftSup, rightSup
        rods = self.Rods
        rodsCount = len(rods)
        A = np.zeros((rodsCount+1, rodsCount+1))
        for i in range(rodsCount):
            KK = (E * rods[i].A) / rods[i].L;
            A[i][i] += KK;
            A[i+1][i+1] += KK;
            A[i+1][i] -= KK;
            A[i][i+1] -= KK;
        if (leftSup):
            A[0][0] = 1;
            A[0][1] = 0;
            A[1][0] = 0;
        if (rightSup):
            A[rodsCount][rodsCount] = 1;
            A[rodsCount][rodsCount-1] = 0;
            A[rodsCount-1][rodsCount] = 0;
        return A;

    def bCalculate(self):
        global E, leftSup, rightSup
        rods = self.Rods
        forces = self.Forces
        rodsCount = len(rods)
        b = np.zeros(rodsCount+1)
        for i in range(rodsCount):
            b[i] = forces[i]

        for i in range(rodsCount):
            if (rods[i].q != 0):
                b[i] += rods[i].q * rods[i].L / 2;
                b[i+1] += rods[i].q * rods[i].L / 2;

        if leftSup:
            b[0] = 0
        if rightSup:
            b[rodsCount] = 0
        return b

    def N(self,i,x):
        rod = self.Rods[i]
        global E, delta
        return (E*rod.A) * ( delta[i+1] - delta[i] ) / rod.L + (rod.q * rod.L) /2 * ( 1 - 2*x/rod.L)
    def U(self,i,x):
        rod = self.Rods[i]
        global E, delta
        return (delta[i] + (x / rod.L) * (delta[i + 1] - delta[i]) + (rod.q * rod.L * rod.L * x * (1 - x / rod.L)) / (
                2 * E * rod.A * rod.L))
    def proc(self):
        A = self.ACalculate()
        #A = getTranspose(A);
        b = self.bCalculate()
        global delta
        delta = linalg.solve(A, b);


class MyNGraph(wx.MDIChildFrame):
    def N(self, x, pr):
        global arr_rods;
        arr = arr_rods
        cur_rod_num = 0
        cur_l = 0.
        while (not ((cur_l <= x) and (cur_l + float(arr_rods[(cur_rod_num) * 4 + 1]) > x))):
            if (cur_rod_num + 2 != len(delta) // 4):
                cur_l += float(arr_rods[cur_rod_num * 4 + 1])
                cur_rod_num += 1
        return pr.N(cur_rod_num, x - cur_l)

    def drawN(self, pr):
        global arr_rods;
        total_l = 0
        for i in range(len(arr_rods) // 4):
            total_l += float(arr_rods[i * 4 + 1])
        data1 = np.arange(0, total_l, 0.001)
        data1.shape = (int(total_l / 0.002), 2)
        # data1[:, 1] = self.N(data1[:, 0], pr)
        for i in range(len(data1)):
            data1[i][1] = self.N(data1[i][0], pr)
        markers1 = PolyMarker(data1, legend='Green Markers', colour='green')  # , marker='circle', size=1)
        return PlotGraphics([markers1], "График N(L)", "L", "N(L)")

    def __init__(self, pr, parent, id, name):
        wx.MDIChildFrame.__init__(self, parent, id,
                                  name)

        # Добавляем панель, чтобы она корректно отображалась на всех платформах.
        panel = wx.Panel(self, wx.ID_ANY)

        # Создаём некоторые сайзеры
        mainSizer = wx.BoxSizer(wx.VERTICAL)
        checkSizer = wx.BoxSizer(wx.HORIZONTAL)

        # Создаём виджеты
        self.canvas = PlotCanvas(panel)
        self.canvas.Draw(self.drawN(pr))
        toggleGrid = wx.CheckBox(panel, label="Show Grid")
        toggleGrid.Bind(wx.EVT_CHECKBOX, self.onToggleGrid)
        toggleLegend = wx.CheckBox(panel, label="Show Legend")
        toggleLegend.Bind(wx.EVT_CHECKBOX, self.onToggleLegend)

        # Размещаем виджеты
        mainSizer.Add(self.canvas, 1, wx.EXPAND)
        checkSizer.Add(toggleGrid, 0, wx.ALL, 5)
        checkSizer.Add(toggleLegend, 0, wx.ALL, 5)
        mainSizer.Add(checkSizer)
        panel.SetSizer(mainSizer)

    def onToggleGrid(self, event):
        """"""
        self.canvas.SetEnableGrid(event.IsChecked())

    def onToggleLegend(self, event):
        """"""
        self.canvas.SetEnableLegend(event.IsChecked())


class MyNGraph(wx.MDIChildFrame):
    def N(self, x, pr):
        global arr_rods;
        arr = arr_rods
        cur_rod_num = 0
        cur_l = 0.
        while (not ((cur_l <= x) and (cur_l + float(arr_rods[(cur_rod_num) * 4 + 1]) >= x))):
            if (cur_rod_num + 2 != len(delta) // 4):
                cur_l += float(arr_rods[cur_rod_num * 4 + 1])
                cur_rod_num += 1
        return pr.N(cur_rod_num, x - cur_l)

    def drawN(self, pr):
        global arr_rods;
        total_l = 0
        for i in range(len(arr_rods) // 4):
            total_l += float(arr_rods[i * 4 + 1])
        data1 = np.arange(0, total_l, 0.001)
        data1.shape = (int(total_l / 0.002), 2)
        # data1[:, 1] = self.N(data1[:, 0], pr)
        for i in range(len(data1)):
            data1[i][1] = self.N(data1[i][0], pr)
        markers1 = PolyMarker(data1, legend='Green Markers', colour='green')  # , marker='circle', size=1)
        return PlotGraphics([markers1], "График N(L)", "L", "N(L)")

    def __init__(self, pr, parent, id, name):
        wx.MDIChildFrame.__init__(self, parent, id,
                                  name)

        # Добавляем панель, чтобы она корректно отображалась на всех платформах.
        panel = wx.Panel(self, wx.ID_ANY)

        # Создаём некоторые сайзеры
        mainSizer = wx.BoxSizer(wx.VERTICAL)
        checkSizer = wx.BoxSizer(wx.HORIZONTAL)

        # Создаём виджеты
        self.canvas = PlotCanvas(panel)
        self.canvas.Draw(self.drawN(pr))
        toggleGrid = wx.CheckBox(panel, label="Show Grid")
        toggleGrid.Bind(wx.EVT_CHECKBOX, self.onToggleGrid)
        toggleLegend = wx.CheckBox(panel, label="Show Legend")
        toggleLegend.Bind(wx.EVT_CHECKBOX, self.onToggleLegend)

        # Размещаем виджеты
        mainSizer.Add(self.canvas, 1, wx.EXPAND)
        checkSizer.Add(toggleGrid, 0, wx.ALL, 5)
        checkSizer.Add(toggleLegend, 0, wx.ALL, 5)
        mainSizer.Add(checkSizer)
        panel.SetSizer(mainSizer)

    def onToggleGrid(self, event):
        """"""
        self.canvas.SetEnableGrid(event.IsChecked())

    def onToggleLegend(self, event):
        """"""
        self.canvas.SetEnableLegend(event.IsChecked())


class MySigGraph(wx.MDIChildFrame):
    def Sig(self, x, pr):
        global arr_rods;
        arr = arr_rods
        cur_rod_num = 0
        cur_l = 0.
        while (not ((cur_l <= x) and (cur_l + float(arr_rods[(cur_rod_num) * 4 + 1]) > x))):
            if (cur_rod_num + 2 != len(delta) // 4):
                cur_l += float(arr_rods[cur_rod_num * 4 + 1])
                cur_rod_num += 1
        return pr.N(cur_rod_num, x - cur_l)/float(arr_rods[cur_rod_num*4+2])

    def drawSig(self, pr):
        global arr_rods;
        total_l = 0
        for i in range(len(arr_rods) // 4):
            total_l += float(arr_rods[i * 4 + 1])
        data1 = np.arange(0, total_l, 0.001)
        data1.shape = (int(total_l / 0.002), 2)
        # data1[:, 1] = self.N(data1[:, 0], pr)
        for i in range(len(data1)):
            data1[i][1] = self.Sig(data1[i][0], pr)
        markers1 = PolyMarker(data1, legend='Green Markers', colour='blue')  # , marker='circle', size=1)
        return PlotGraphics([markers1], "График sigma(L)", "L", "sigma(L)")

    def __init__(self, pr, parent, id, name):
        wx.MDIChildFrame.__init__(self, parent, id,
                                  name)

        # Добавляем панель, чтобы она корректно отображалась на всех платформах.
        panel = wx.Panel(self, wx.ID_ANY)

        # Создаём некоторые сайзеры
        mainSizer = wx.BoxSizer(wx.VERTICAL)
        checkSizer = wx.BoxSizer(wx.HORIZONTAL)

        # Создаём виджеты
        self.canvas = PlotCanvas(panel)
        self.canvas.Draw(self.drawSig(pr))
        toggleGrid = wx.CheckBox(panel, label="Show Grid")
        toggleGrid.Bind(wx.EVT_CHECKBOX, self.onToggleGrid)
        toggleLegend = wx.CheckBox(panel, label="Show Legend")
        toggleLegend.Bind(wx.EVT_CHECKBOX, self.onToggleLegend)

        # Размещаем виджеты
        mainSizer.Add(self.canvas, 1, wx.EXPAND)
        checkSizer.Add(toggleGrid, 0, wx.ALL, 5)
        checkSizer.Add(toggleLegend, 0, wx.ALL, 5)
        mainSizer.Add(checkSizer)
        panel.SetSizer(mainSizer)

    def onToggleGrid(self, event):
        """"""
        self.canvas.SetEnableGrid(event.IsChecked())

    def onToggleLegend(self, event):
        """"""
        self.canvas.SetEnableLegend(event.IsChecked())
class MyUGraph(wx.MDIChildFrame):
    def U(self, x, pr):
        global arr_rods;
        arr = arr_rods
        cur_rod_num = 0
        cur_l = 0.
        while( not((cur_l <= x) and (cur_l + float(arr_rods[(cur_rod_num)*4+1])>x))):
            if (cur_rod_num+2 != len(delta)//4):
                cur_l += float(arr_rods[cur_rod_num*4+1])
                cur_rod_num+=1
        return pr.U(cur_rod_num, x - cur_l)

    def drawU(self,pr):
        global arr_rods;
        total_l = 0
        for i in range(len(arr_rods) // 4):
            total_l += float(arr_rods[i*4+1])
        data1 = np.arange(0, total_l, 0.001)
        data1.shape = (int(total_l/0.002), 2)
        #data1[:, 1] = self.N(data1[:, 0], pr)
        for i in range(len(data1)):
            data1[i][1] = self.U(data1[i][0], pr)
        markers1 = PolyMarker(data1, legend='Green Markers', colour='red')#, marker='circle', size=1)
        return PlotGraphics([markers1], "График U(L)", "L", "U(L)")

    def __init__(self, pr,parent,id,name):
        wx.MDIChildFrame.__init__(self, parent, id,
                          name)

        # Добавляем панель, чтобы она корректно отображалась на всех платформах.
        panel = wx.Panel(self, wx.ID_ANY)

        # Создаём некоторые сайзеры
        mainSizer = wx.BoxSizer(wx.VERTICAL)
        checkSizer = wx.BoxSizer(wx.HORIZONTAL)

        # Создаём виджеты
        self.canvas = PlotCanvas(panel)
        self.canvas.Draw(self.drawU(pr))
        toggleGrid = wx.CheckBox(panel, label="Show Grid")
        toggleGrid.Bind(wx.EVT_CHECKBOX, self.onToggleGrid)
        toggleLegend = wx.CheckBox(panel, label="Show Legend")
        toggleLegend.Bind(wx.EVT_CHECKBOX, self.onToggleLegend)

        # Размещаем виджеты
        mainSizer.Add(self.canvas, 1, wx.EXPAND)
        checkSizer.Add(toggleGrid, 0, wx.ALL, 5)
        checkSizer.Add(toggleLegend, 0, wx.ALL, 5)
        mainSizer.Add(checkSizer)
        panel.SetSizer(mainSizer)

    def onToggleGrid(self, event):
        """"""
        self.canvas.SetEnableGrid(event.IsChecked())

    def onToggleLegend(self, event):
        """"""
        self.canvas.SetEnableLegend(event.IsChecked())
class window2(wx.MDIParentFrame):
    def __init__(self, parent):
        super().__init__(parent, pos=(0, 0), size=(700, 400), name = 'Постпроцессор')
        self.pr = Proc()

        self.pr.proc()

        #self.Bind(wx.EVT_PAINT, self.OnPaint)
        self.Maximize(True)
        menubar = wx.MenuBar()
        self.SetMenuBar(menubar)


        win = wx.MDIChildFrame(self, -1, "Меню поспроцессора", size=(800, 500))
        win.fSizer = wx.BoxSizer(wx.VERTICAL)
        win.Show()
        panel = PostProcMenu(win)
        win.fSizer.Add(panel, 1, wx.EXPAND)
        win.SetSizer(win.fSizer)
        win.Show()

        #win2 = wx.MDIChildFrame(self, -1, "Child Window", size=(200, 150))
        #win2.Show()
        #self.Sizer = wx.BoxSizer(wx.VERTICAL)
        #self.SetSizer(self.Sizer)

        #info_text = wx.StaticText(self, label="Введите:")

        #self.fSizer.Add(info_text, wx.EXPAND | wx.LEFT, 5)

        #E_sizer = wx.BoxSizer(wx.HORIZONTAL)
        #E_mpa_text = wx.StaticText(self, label="E (МПа):")
        #E_mpa = wx.TextCtrl(self)
        #E_sizer.Add(E_mpa_text, wx.EXPAND | wx.LEFT, 5)
        #E_sizer.Add(E_mpa, wx.EXPAND | wx.LEFT, 5)
        #self.Sizer.Add(E_sizer,  wx.LEFT, 5)

    def NewNTable(self, row_num, L_start, L_finish, step):
        win = wx.MDIChildFrame(self, -1, "Таблица значений N для стержня " + str(int(row_num) + 1), size=(200, 150))
        win.fSizer = wx.BoxSizer(wx.VERTICAL)
        win.Show()
        panel = N_table(win, row_num, L_start, L_finish, step, self.pr)
        win.fSizer.Add(panel, 1, wx.EXPAND)
        win.SetSizer(win.fSizer)
        win.Show()
    def NewUTable(self, row_num, L_start, L_finish, step):
        win = wx.MDIChildFrame(self, -1, "Таблица значений U для стержня " + str(int(row_num)+1), size=(200, 150))
        win.fSizer = wx.BoxSizer(wx.VERTICAL)
        win.Show()
        panel = U_table(win, row_num, L_start, L_finish, step, self.pr)
        win.fSizer.Add(panel, 1, wx.EXPAND)
        win.SetSizer(win.fSizer)
        win.Show()
    def NewSigTable(self, row_num, L_start, L_finish, step):
        win = wx.MDIChildFrame(self, -1, "Таблица значений sigma для стержня " + str(int(row_num)+1), size=(200, 150))
        win.fSizer = wx.BoxSizer(wx.VERTICAL)
        win.Show()
        panel = Sig_table(win, row_num, L_start, L_finish, step, self.pr)
        win.fSizer.Add(panel, 1, wx.EXPAND)
        win.SetSizer(win.fSizer)
        win.Show()

    def NewNGraf(self, e):
        win = MyNGraph(self.pr, self, id=-1, name="График N для стержня ")
        win.Show()
    def NewSigGraf(self,e):
        win = MySigGraph(self.pr,self, id =-1, name="График N для стержня ")
        win.Show()
    def NewUGraf(self,e):
        win = MyUGraph(self.pr,self, id =-1, name="График N для стержня ")
        win.Show()

class N_table(wx.Window):
    def __init__(self, parent, row_num, L_start, L_finish, step, pr):
        """Constructor"""
        wx.ScrolledWindow.__init__(self, parent)
        self.number_of_buttons = 0
        self.frame = parent


        self.mainSizer = wx.BoxSizer(wx.VERTICAL)

        self.tab = wx.grid.Grid(self)
        self.tab.CreateGrid((L_finish-L_start)//step + 1, 2)
        self.write(row_num, L_start, L_finish, step, pr)

        #self.tab.SetCellValue(0,1,'111')
        if ( ((L_finish-L_start)//step + 5)< 20):
            self.tab.SetMinSize((259, 20* ((L_finish-L_start)//step ) + 5) )
            self.GetParent().SetMaxSize((272, 20* ((L_finish-L_start)//step+2) + 5))
        else:
            self.tab.SetMinSize((259, 20 * 22 + 5))
            self.GetParent().SetMaxSize((272, 22 * 20 + 15))
        #self.GetParent().Maximize(True)

        #self.GetParent().SetSizeHints(-1,270, 20* (L_finish-L_start)//step + 5)


        self.SetLayout()
        self.Centre()


    def write(self, row_num, L_start, L_finish, step, pr):
        self.tab.SetColLabelValue(0, "L")
        self.tab.SetColLabelValue(1, "N")

        i = 0
        for j in np.arange(L_start, L_finish, step):
            self.tab.SetCellValue(i, 0, str(j))
            self.tab.SetCellValue(i, 1, str(pr.N(int(row_num), j)))
            i += 1
    def SetLayout(self):
        sizer = wx.GridBagSizer(5, 5)

        sizer.Add(self.tab, pos=(0, 0), span=(0, 0), flag=wx.TOP | wx.EXPAND, border=0)
        self.SetSizer(sizer)
        self.SetBackgroundColour("#D8BFD8")
class U_table( N_table):
    def write(self, row_num, L_start, L_finish, step, pr):
        self.tab.SetColLabelValue(0, "L")
        self.tab.SetColLabelValue(1, "U")
        global sigma
        i = 0
        for j in np.arange(L_start, L_finish, step):
            self.tab.SetCellValue(i, 0, str(j))
            o=''
            if (pr.N(int(row_num), j) / float(arr_rods[int(row_num) * 4 + 2]) > sigma):
                o = '*'
            self.tab.SetCellValue(i, 1, str(pr.U(int(row_num), j)) + o)
            i+=1
class Sig_table( N_table):
    def write(self, row_num, L_start, L_finish, step, pr):
        self.tab.SetColLabelValue(0, "L")
        self.tab.SetColLabelValue(1, "Sig")

        i = 0
        for j in np.arange(L_start, L_finish, step):
            self.tab.SetCellValue(i, 0, str(j))
            global arr_rods
            self.tab.SetCellValue(i, 1, str(pr.N(int(row_num), j)/float(arr_rods[int(row_num)*4+2])))
            i+=1
class PostProcMenu(wx.ScrolledWindow):
    def __init__(self, parent):
        """Constructor"""
        wx.ScrolledWindow.__init__(self, parent)
        self.number_of_buttons = 0
        self.frame = parent

        self.scrollStep = 10
        self.SetScrollbars(self.scrollStep,
                           self.scrollStep,
                           self.GetMaxWidth()/ self.scrollStep,
                           self.GetMaxWidth() / self.scrollStep)

        self.mainSizer = wx.BoxSizer(wx.VERTICAL)
        staticbox = wx.StaticBox (self, wx.NewId(), label="Таблицы")
        hatSizer = wx.StaticBoxSizer( staticbox, wx.VERTICAL)
        self.mainSizer.Add(hatSizer, wx.LEFT, 5)

        n_table_sizer = wx.BoxSizer(wx.HORIZONTAL)
        n_table_btn = wx.Button(self, label="Coздать таблицу N(x): ")
        choices = []
        global arr_rods
        for i in range(len(arr_rods)//4):
            choices.append(str(i+1))
        rod_num = wx.Choice(self,choices= choices, name = 'Номер стержня')
        rod_num.Bind(wx.EVT_CHOICE, self.on_choice)
        rod_num_text = wx.StaticText ( self, label ="№ cтержня")
        L_start = wx.TextCtrl(self)
        L_start.SetValue('L начала')
        L_finish = wx.TextCtrl(self)
        L_finish.SetValue('L конца')
        L_step = wx.TextCtrl(self)
        L_step.SetValue('Шаг')
        n_table_btn.Bind(wx.EVT_BUTTON, self.NewNTable)
        n_table_sizer.Add(n_table_btn, 0, wx.CENTER | wx.ALL, 5)
        n_table_sizer.Add(rod_num_text, 0, wx.CENTER | wx.ALL, 5)
        n_table_sizer.Add(rod_num, 0, wx.CENTER | wx.ALL, 5)
        n_table_sizer.Add(L_start, 0, wx.CENTER | wx.ALL, 5)
        n_table_sizer.Add(L_finish, 0, wx.CENTER | wx.ALL, 5)
        n_table_sizer.Add(L_step, 0, wx.CENTER | wx.ALL, 5)
        hatSizer.Add(n_table_sizer, wx.EXPAND | wx.LEFT, 5)

        s_table_sizer = wx.BoxSizer(wx.HORIZONTAL)
        s_table_btn = wx.Button(self, label="Coздать таблицу sig(x): ")
        choices = []

        for i in range(len(arr_rods) // 4):
            choices.append(str(i + 1))
        s_rod_num = wx.Choice(self, choices=choices, name='Номер стержня')
        s_rod_num.Bind(wx.EVT_CHOICE, self.on_choice1)
        s_rod_num_text = wx.StaticText(self, label="№ cтержня")
        s_L_start = wx.TextCtrl(self)
        s_L_start.SetValue('L начала')
        s_L_finish = wx.TextCtrl(self)
        s_L_finish.SetValue('L конца')
        s_L_step = wx.TextCtrl(self)
        s_L_step.SetValue('Шаг')
        s_table_btn.Bind(wx.EVT_BUTTON, self.NewSigTable)
        s_table_sizer.Add(s_table_btn, 0, wx.CENTER | wx.ALL, 5)
        s_table_sizer.Add(s_rod_num_text, 0, wx.CENTER | wx.ALL, 5)
        s_table_sizer.Add(s_rod_num, 0, wx.CENTER | wx.ALL, 5)
        s_table_sizer.Add(s_L_start, 0, wx.CENTER | wx.ALL, 5)
        s_table_sizer.Add(s_L_finish, 0, wx.CENTER | wx.ALL, 5)
        s_table_sizer.Add(s_L_step, 0, wx.CENTER | wx.ALL, 5)
        hatSizer.Add(s_table_sizer, wx.EXPAND | wx.LEFT, 5)

        u_table_sizer = wx.BoxSizer(wx.HORIZONTAL)
        u_table_btn = wx.Button(self, label="Coздать таблицу U(x): ")
        choices = []

        for i in range(len(arr_rods) // 4):
            choices.append(str(i + 1))
        u_rod_num = wx.Choice(self, choices=choices, name='Номер стержня')
        u_rod_num.Bind(wx.EVT_CHOICE, self.on_choice2)
        u_rod_num_text = wx.StaticText(self, label="№ cтержня")
        u_L_start = wx.TextCtrl(self)
        u_L_start.SetValue('L начала')
        u_L_finish = wx.TextCtrl(self)
        u_L_finish.SetValue('L конца')
        u_L_step = wx.TextCtrl(self)
        u_L_step.SetValue('Шаг')
        u_table_btn.Bind(wx.EVT_BUTTON, self.NewUTable)
        u_table_sizer.Add(u_table_btn, 0, wx.CENTER | wx.ALL, 5)
        u_table_sizer.Add(u_rod_num_text, 0, wx.CENTER | wx.ALL, 5)
        u_table_sizer.Add(u_rod_num, 0, wx.CENTER | wx.ALL, 5)
        u_table_sizer.Add(u_L_start, 0, wx.CENTER | wx.ALL, 5)
        u_table_sizer.Add(u_L_finish, 0, wx.CENTER | wx.ALL, 5)
        u_table_sizer.Add(u_L_step, 0, wx.CENTER | wx.ALL, 5)
        hatSizer.Add(u_table_sizer, wx.EXPAND | wx.LEFT, 5)

        staticbox2 = wx.StaticBox(self, wx.NewId(), label="Графики")
        hatSizer2 = wx.StaticBoxSizer(staticbox2, wx.HORIZONTAL)
        self.mainSizer.Add(hatSizer2, wx.LEFT, 5)

        N_graph_btn = wx.Button(self, label="Coздать график N(L): ")
        N_graph_btn.Bind(wx.EVT_BUTTON, self.GetParent().GetParent().NewNGraf)
        hatSizer2.Add(N_graph_btn,  wx.LEFT, 5)
        Sig_graph_btn = wx.Button(self, label="Coздать график sig(L): ")
        Sig_graph_btn.Bind(wx.EVT_BUTTON, self.GetParent().GetParent().NewSigGraf)
        hatSizer2.Add(Sig_graph_btn,  wx.LEFT, 5)
        U_graph_btn = wx.Button(self, label="Coздать график U(x): ")
        U_graph_btn.Bind(wx.EVT_BUTTON, self.GetParent().GetParent().NewUGraf)
        hatSizer2.Add(U_graph_btn,  wx.LEFT, 5)
        staticbox3 = wx.StaticBox(self, wx.NewId(), label="Анализатор")
        hatSizer3 = wx.StaticBoxSizer(staticbox3, wx.HORIZONTAL)
        self.mainSizer.Add(hatSizer3, wx.LEFT, 5)

        L_text = wx.StaticText(self, label="L = ")
        L = wx.TextCtrl(self)
        L.Bind(wx.EVT_KILL_FOCUS, self.L_is_entered)
        hatSizer3.Add(L_text,  0, wx.CENTER | wx.ALL, 5)
        hatSizer3.Add(L,  0, wx.CENTER | wx.ALL, 5)
        N_text = wx.StaticText(self, label="N(L) = ")
        N = wx.TextCtrl(self)
        hatSizer3.Add(N_text,  0, wx.CENTER | wx.ALL, 5)
        hatSizer3.Add(N,  0, wx.CENTER | wx.ALL, 5)
        sig_text = wx.StaticText(self, label="sig(L) = ")
        sig = wx.TextCtrl(self)
        hatSizer3.Add(sig_text,  0, wx.CENTER | wx.ALL, 5)
        hatSizer3.Add(sig,  0, wx.CENTER | wx.ALL, 5)
        U_text = wx.StaticText(self, label="U(L) = ")
        U = wx.TextCtrl(self)
        hatSizer3.Add(U_text,  0, wx.CENTER | wx.ALL, 5)
        hatSizer3.Add(U,  0, wx.CENTER | wx.ALL, 5)


        self.SetSizer(self.mainSizer)
    def L_is_entered(self, e):
        sizer = e.GetEventObject().GetParent().Sizer.GetChildren()[2].Sizer.GetChildren()
        cur_l = 0
        cur_rod_num = 0
        global arr_rods
        total_l = 0
        for i in range(len(arr_rods) // 4):
            total_l += float(arr_rods[i * 4 + 1])
        try:
            x = float(sizer[1].Window.GetValue())
        except ValueError:
            sizer[1].Window.SetValue('')
            e.Skip()
            return
        if ( x < 0 ) or ( x > total_l):
            sizer[1].Window.SetValue('')
            e.Skip()
        while (not ((cur_l <= x) and (cur_l + float(arr_rods[(cur_rod_num) * 4 + 1]) > x))):
            if (cur_rod_num + 2 != len(delta) // 4):
                cur_l += float(arr_rods[cur_rod_num * 4 + 1])
                cur_rod_num += 1
        sizer[3].Window.SetValue (str(self.GetParent().GetParent().pr.N(cur_rod_num, x - cur_l)))
        sizer[5].Window.SetValue(str(self.GetParent().GetParent().pr.N(cur_rod_num, x - cur_l)/float(arr_rods[int(cur_rod_num)*4+2])))
        sizer[7].Window.SetValue(str(self.GetParent().GetParent().pr.U(cur_rod_num, x - cur_l)))
        e.Skip()

    def on_choice(self, e):
        sizer = e.GetEventObject().GetParent().Sizer.GetChildren()[0].Sizer.GetChildren()[0].Sizer.GetChildren()
        global arr_rods
        x=sizer[3]
        x=sizer[4]
        sizer[3].Window.SetValue('0')
        x=e.GetEventObject().GetSelection()
        sizer[4].Window.SetValue(str(arr_rods[x*4+1]))
    def on_choice1(self, e):
        sizer = e.GetEventObject().GetParent().Sizer.GetChildren()[0].Sizer.GetChildren()[1].Sizer.GetChildren()
        global arr_rods
        x=sizer[3]
        sizer[3].Window.SetValue('0')
        sizer[4].Window.SetValue(str(arr_rods[int((e.GetEventObject().GetSelection()))*4+1]))
    def on_choice2(self, e):
        sizer = e.GetEventObject().GetParent().Sizer.GetChildren()[0].Sizer.GetChildren()[2].Sizer.GetChildren()
        global arr_rods
        x=sizer[3]
        sizer[3].Window.SetValue('0')
        sizer[4].Window.SetValue(str(arr_rods[int((e.GetEventObject().GetSelection()))*4+1]))
    def NewNTable(self,e):
        sizer = e.GetEventObject().GetParent().Sizer.GetChildren()[0].Sizer.GetChildren()[0].Sizer.GetChildren()
        self.GetParent().GetParent().NewNTable(float(sizer[2].Window.GetSelection()), float(sizer[3].Window.GetValue()), float(sizer[4].Window.GetValue()), float(sizer[5].Window.GetValue()))
    def NewSigTable(self,e):
        sizer = e.GetEventObject().GetParent().Sizer.GetChildren()[0].Sizer.GetChildren()[1].Sizer.GetChildren()
        self.GetParent().GetParent().NewSigTable(float(sizer[2].Window.GetSelection()), float(sizer[3].Window.GetValue()), float(sizer[4].Window.GetValue()), float(sizer[5].Window.GetValue()))
    def NewUTable(self,e):
        sizer = e.GetEventObject().GetParent().Sizer.GetChildren()[0].Sizer.GetChildren()[2].Sizer.GetChildren()
        self.GetParent().GetParent().NewUTable(float(sizer[2].Window.GetSelection()), float(sizer[3].Window.GetValue()), float(sizer[4].Window.GetValue()), float(sizer[5].Window.GetValue()))


class MyPanel(wx.ScrolledWindow):

    E = 0
    sig = 0

    def __init__(self, parent):
        """Constructor"""
        wx.ScrolledWindow.__init__(self, parent)
        self.number_of_buttons = 0
        self.frame = parent

        self.scrollStep = 10
        self.SetScrollbars(self.scrollStep,
                           self.scrollStep,
                           self.GetMaxWidth()/ self.scrollStep,
                           self.GetMaxWidth() / self.scrollStep)

        self.mainSizer = wx.BoxSizer(wx.VERTICAL)


        staticbox = wx.StaticBox (self, wx.NewId(), label="Общая информация о конструкции")
        hatSizer = wx.StaticBoxSizer( staticbox, wx.VERTICAL)
        controlSizer = wx.BoxSizer( wx.HORIZONTAL)
        self.widgetSizer = wx.BoxSizer(wx.VERTICAL)

        info_text = wx.StaticText ( self, label ="Введите:")

        hatSizer.Add(info_text, wx.EXPAND | wx.LEFT, 5 )

        E_sizer = wx.BoxSizer(wx.HORIZONTAL)
        E_mpa_text = wx.StaticText ( self, label = "E (МПа):")
        E_mpa = wx.TextCtrl (self)
        E_mpa.Bind( wx.EVT_KILL_FOCUS, self.E_is_entered )
        E_sizer.Add( E_mpa_text, wx.EXPAND | wx.LEFT, 5 )
        E_sizer.Add(E_mpa, wx.EXPAND | wx.LEFT, 5)

        hatSizer.Add(E_sizer, wx.EXPAND | wx.ALL, 5 )

        sigma_sizer = wx.BoxSizer(wx.HORIZONTAL)
        sigma_text = wx.StaticText ( self, label = "sig (МПа:")
        sigma = wx.TextCtrl(self)
        sigma.Bind(wx.EVT_KILL_FOCUS, self.sig_is_entered)
        sigma_sizer.Add(sigma_text, wx.EXPAND | wx.LEFT, 5)
        sigma_sizer.Add(sigma, wx.EXPAND | wx.LEFT, 5)

        hatSizer.Add(sigma_sizer, wx.EXPAND | wx.ALL, 5)

        checkbox_sizer = wx.BoxSizer(wx.HORIZONTAL)
        support1 = wx.CheckBox ( self, label = "Левая опора")
        support2 = wx.CheckBox ( self, label = "Правая опора")
        checkbox_sizer.Add(support1, wx.EXPAND | wx.LEFT, 5)
        checkbox_sizer.Add(support2, wx.EXPAND | wx.LEFT, 5)

        hatSizer.Add(checkbox_sizer, wx.EXPAND | wx.ALL, 5)

        self.addButton = wx.Button(self, label="Добавить стержень")
        self.addButton.Bind(wx.EVT_BUTTON, self.onAddWidget)
        controlSizer.Add(self.addButton, 0, wx.CENTER | wx.ALL, 5)

        self.removeButton = wx.Button(self, label="Удалить стержень")
        self.removeButton.Bind(wx.EVT_BUTTON, self.onRemoveWidget)
        controlSizer.Add(self.removeButton, 0, wx.CENTER | wx.ALL, 5)

        self.calcButton = wx.Button(self, label="Рассчитать")
        self.calcButton.Bind(wx.EVT_BUTTON, self.onCalculate)
        controlSizer.Add(self.calcButton, 0, wx.CENTER | wx.ALL, 5)

        self.saveButton = wx.Button(self, label="Сохранить :")
        save_path = wx.TextCtrl(self)
        self.saveButton.Bind(wx.EVT_BUTTON, self.save_to_file)
        controlSizer.Add(self.saveButton, 0, wx.CENTER | wx.ALL, 5)
        controlSizer.Add(save_path, 0, wx.CENTER | wx.ALL, 5)

        self.loadButton = wx.Button(self, label="Загрузить :")
        load_path = wx.TextCtrl(self)
        self.loadButton.Bind(wx.EVT_BUTTON, self.open_file)
        controlSizer.Add(self.loadButton, 0, wx.CENTER | wx.ALL, 5)
        controlSizer.Add(load_path, 0, wx.CENTER | wx.ALL, 5)


        node1_sizer = wx.BoxSizer(wx.HORIZONTAL)
        node1_text = wx.StaticText(self, label="Узел 1, сосредоточенная нагрузка(H):")
        node1 = wx.TextCtrl(self)
        node1_sizer.Add(node1_text, wx.EXPAND | wx.LEFT, 5)
        node1_sizer.Add(node1, wx.EXPAND | wx.LEFT, 5)

        self.mainSizer.Add(hatSizer, 0, wx.LEFT)
        self.mainSizer.Add(controlSizer, 0, wx.LEFT)
        self.mainSizer.Add(node1_sizer, 0, wx.LEFT | wx.DOWN, 10)
        self.mainSizer.Add(self.widgetSizer, 0, wx.LEFT, 10)
        self.SetSizer(self.mainSizer)

    def onCalculate(self, event):
        secondWindow = window2(parent=None)
        secondWindow.Show()
        self.GetParent().GetParent().Destroy()
    def onAddWidget(self, event):
        self.number_of_buttons += 1
        segment_sizer = wx.BoxSizer(wx.VERTICAL)
        kernel_sizer = wx.BoxSizer(wx.HORIZONTAL)

        new_kernel_text = wx.StaticText( self, label = "Стержень %s" % self.number_of_buttons)
        kernel_sizer.Add(new_kernel_text, 0, wx.LEFT)
        lenght_text = wx.StaticText(self, label="Длинна стержня L (m): ")
        lenght = wx.TextCtrl(self)
        lenght.Bind(wx.EVT_KILL_FOCUS, self.lenght_is_entered)
       #lenght.Bind(wx.EVT_SET_FOCUS, self.onTextFocus)
        kernel_sizer.Add(lenght_text, 0, wx.LEFT, 30)
        kernel_sizer.Add(lenght, 0, wx.LEFT)
        square_text = wx.StaticText(self, label="Площадь поперечного сечения A (m^2): ")
        square = wx.TextCtrl(self)
        square.Bind(wx.EVT_KILL_FOCUS, self.square_is_entered)
        kernel_sizer.Add(square_text, 0, wx.LEFT, 20)
        kernel_sizer.Add(square, 0, wx.LEFT)
        q_text = wx.StaticText(self, label="Распределённая нагрузка q (kH/m): ")
        q = wx.TextCtrl(self)
        kernel_sizer.Add(q_text, 0, wx.LEFT, 20)
        kernel_sizer.Add(q, 0, wx.LEFT)

        node_sizer = wx.BoxSizer(wx.HORIZONTAL)
        node_text = wx.StaticText(self, label = "Узел %s, сосредоточенная нагрузка(H):" % (self.number_of_buttons + 1))
        node = wx.TextCtrl(self)
        node_sizer.Add(node_text, wx.EXPAND | wx.LEFT)
        node_sizer.Add(node, wx.EXPAND | wx.LEFT, 20)
        node.Bind(wx.EVT_KILL_FOCUS, self.drow_pic)

        segment_sizer.Add(kernel_sizer, 0, wx.DOWN, 5)
        segment_sizer.Add(node_sizer, 0, wx.DOWN, 15)

        self.widgetSizer.Add(segment_sizer, 0, wx.EXPAND | wx.ALL, 5)
        self.mainSizer.Layout()
        #self.widgetSizer.Layout()
        #self.frame.Fit()

    def onRemoveWidget(self, event):

        if self.widgetSizer.GetChildren():
            self.widgetSizer.Hide(self.number_of_buttons - 1)
            self.widgetSizer.Remove(self.number_of_buttons - 1)
            self.number_of_buttons -= 1
            #self.frame.GetParent().GetParent().fSizer.Layout()
            #self.frame.Fit()
        self.drow_pic(event)

    def E_is_entered(self, event):
        if event.GetEventObject().GetValue() == '':
            event.Skip()
            return
        if (int(event.GetEventObject().GetValue()) <= 0):
            event.GetEventObject().Clear()
            event.Skip()
            return
        global E
        E = int(event.GetEventObject().GetValue())
        event.Skip()
    def sig_is_entered(self, event):
        if event.GetEventObject().GetValue() == '':
            event.Skip()
            return
        if (int(event.GetEventObject().GetValue()) <= 0):
            event.GetEventObject().Clear()
            event.Skip()
            return
        global sigma
        sigma = int(event.GetEventObject().GetValue())
        event.Skip()

    def lenght_is_entered(self, event):
        try:
            l = int(event.GetEventObject().GetValue())
        except ValueError:
            event.GetEventObject().Clear()
            event.Skip()
            return
        if ( int(event.GetEventObject().GetValue()) <= 0):
            event.GetEventObject().Clear()
        event.Skip()
    def square_is_entered(self, event):
        if ( int(event.GetEventObject().GetValue()) <= 0):
            event.GetEventObject().Clear()
        event.Skip()

    def open_file(self, event):
        while (self.widgetSizer.GetChildren()):
            self.onRemoveWidget(event)
        open_path = self.mainSizer.GetChildren()[1].Sizer.GetChildren()[6].Window.GetValue()
        str = []
        with open(open_path, 'r') as file:
            str = file.readline()
        arr = str.split()
        self.mainSizer.GetChildren()[0].Sizer.GetChildren()[1].Sizer.GetChildren()[1].Window.SetValue(arr[0])
        global E, sigma
        E = float(arr[0])
        sigma = float(arr[1])
        self.mainSizer.GetChildren()[0].Sizer.GetChildren()[2].Sizer.GetChildren()[1].Window.SetValue(arr[1])
        self.mainSizer.GetChildren()[0].Sizer.GetChildren()[3].Sizer.GetChildren()[0].Window.SetValue(arr[2] == 'True')
        self.mainSizer.GetChildren()[0].Sizer.GetChildren()[3].Sizer.GetChildren()[1].Window.SetValue(arr[3] == 'True')
        self.mainSizer.GetChildren()[2].Sizer.GetChildren()[1].Window.SetValue(arr[4])
        for i in range((len(arr)- 5)//4):
            self.onAddWidget(event)
            segment = self.widgetSizer.GetChildren()[i]
            segment.Sizer.GetChildren()[0].Sizer.GetChildren()[2].Window.SetValue(arr[i*4 + 5])
            segment.Sizer.GetChildren()[0].Sizer.GetChildren()[4].Window.SetValue(arr[i*4+6])
            segment.Sizer.GetChildren()[0].Sizer.GetChildren()[6].Window.SetValue(arr[i*4+7])

            segment.Sizer.GetChildren()[1].Sizer.GetChildren()[1].Window.SetValue(arr[i*4+8])
        self.drow_pic(event)




    def save_to_file(self, event):
        save_path = self.mainSizer.GetChildren()[1].Sizer.GetChildren()[4].Window.GetValue()
        arr=[]
        arr.append(self.mainSizer.GetChildren()[0].Sizer.GetChildren()[1].Sizer.GetChildren()[1].Window.GetValue())
        arr.append(' ')
        arr.append(self.mainSizer.GetChildren()[0].Sizer.GetChildren()[2].Sizer.GetChildren()[1].Window.GetValue())
        arr.append(' ')
        arr.append(str(self.mainSizer.GetChildren()[0].Sizer.GetChildren()[3].Sizer.GetChildren()[0].Window.GetValue()))
        arr.append(' ')
        arr.append(str(self.mainSizer.GetChildren()[0].Sizer.GetChildren()[3].Sizer.GetChildren()[1].Window.GetValue()))
        arr.append(' ')



        F1 = self.mainSizer.GetChildren()[2].Sizer.GetChildren()[1].Window.GetValue()
        try:
            f1 = int(F1)
        except ValueError:
            self.mainSizer.GetChildren()[2].Sizer.GetChildren()[1].Window.Clear()
        arr.append(F1)
        arr.append(' ')
        for segment in self.widgetSizer.GetChildren():
            L = segment.Sizer.GetChildren()[0].Sizer.GetChildren()[2].Window.GetValue()
            A = segment.Sizer.GetChildren()[0].Sizer.GetChildren()[4].Window.GetValue()
            Q = segment.Sizer.GetChildren()[0].Sizer.GetChildren()[6].Window.GetValue()

            F2 = segment.Sizer.GetChildren()[1].Sizer.GetChildren()[1].Window.GetValue()
            try:
                l = int(L)
                a = int(A)
                q = int(Q)
                f2 = int(F2)
            except ValueError:
                segment.Sizer.GetChildren()[0].Sizer.GetChildren()[2].Window.Clear()
                segment.Sizer.GetChildren()[0].Sizer.GetChildren()[4].Window.Clear()
                segment.Sizer.GetChildren()[0].Sizer.GetChildren()[6].Window.Clear()
                segment.Sizer.GetChildren()[1].Sizer.GetChildren()[1].Window.Clear()
                return

            arr.append(L)
            arr.append(' ')
            arr.append(A)
            arr.append(' ')
            arr.append(Q)
            arr.append(' ')
            arr.append(F2)
            arr.append(' ')

        with open (save_path, 'w') as file:
            file.writelines(arr)



    def drow_pic(self, event):
        list = []
        F1 = self.mainSizer.GetChildren()[2].Sizer.GetChildren()[1].Window.GetValue()
        try:
            f1 = int(F1)
        except ValueError:
            self.mainSizer.GetChildren()[2].Sizer.GetChildren()[1].Window.Clear()
        list.append(F1)
        for segment in self.widgetSizer.GetChildren():
            L = segment.Sizer.GetChildren()[0].Sizer.GetChildren()[2].Window.GetValue()
            A = segment.Sizer.GetChildren()[0].Sizer.GetChildren()[4].Window.GetValue()
            Q = segment.Sizer.GetChildren()[0].Sizer.GetChildren()[6].Window.GetValue()

            F2 = segment.Sizer.GetChildren()[1].Sizer.GetChildren()[1].Window.GetValue()
            try:
                l = int(L)
                a = int(A)
                q = int(Q)
                f2 = int(F2)
            except ValueError:
                segment.Sizer.GetChildren()[0].Sizer.GetChildren()[2].Window.Clear()
                segment.Sizer.GetChildren()[0].Sizer.GetChildren()[4].Window.Clear()
                segment.Sizer.GetChildren()[0].Sizer.GetChildren()[6].Window.Clear()
                segment.Sizer.GetChildren()[1].Sizer.GetChildren()[1].Window.Clear()
                event.Skip()
                return

            list.append(L)
            list.append(A)
            list.append(Q)
            list.append(F2)
        x = self.GetParent().GetChildren()
        y = x[1]
        global arr_rods, leftSup, rightSup
        leftSup = self.mainSizer.GetChildren()[0].Sizer.GetChildren()[3].Sizer.GetChildren()[0].Window.GetValue()
        rightSup = self.mainSizer.GetChildren()[0].Sizer.GetChildren()[3].Sizer.GetChildren()[1].Window.GetValue()
        arr_rods = list.copy()
        y.print_pic(self.mainSizer.GetChildren()[0].Sizer.GetChildren()[3].Sizer.GetChildren()[0].Window.GetValue(),
                    self.mainSizer.GetChildren()[0].Sizer.GetChildren()[3].Sizer.GetChildren()[1].Window.GetValue(),
                    list)
        event.Skip()
        return





class MyFrame(wx.Frame ):

    def __init__(self):
        """Constructor"""
        wx.Frame.__init__(self, parent=None, title="SAPR")
        self.Maximize(True)
        self.fSizer = wx.BoxSizer(wx.HORIZONTAL)

        splitter = wx.SplitterWindow(self)
        splitter.Sizer = wx.BoxSizer(wx.HORIZONTAL)
        panel = MyPanel(splitter)
        panel2 = Canvas(splitter)
        splitter.SplitHorizontally(panel, panel2)
        splitter.SetMinimumPaneSize(600)

        self.fSizer.Add(splitter, 1, wx.EXPAND)
        '''self.fSizer = wx.BoxSizer(wx.VERTICAL)'''



        self.SetSizer(self.fSizer)
        self.Show()


class Canvas(wx.ScrolledWindow):
    """
    Canvas stores and renders all nodes and node connections.
    It also handles all user interaction.
    """
    XZoom = 20
    YZoom = 20
    def __init__(self, *args, **kw):
        super(Canvas, self).__init__(*args, **kw)
        self.scrollStep = kw.get("scrollStep", 10)
        self.canvasDimensions = kw.get("canvasDimensions", [2800, 400])
        self.SetScrollbars(self.scrollStep,
                           self.scrollStep,
                           self.canvasDimensions[0] / self.scrollStep,
                           self.canvasDimensions[1] / self.scrollStep)

        self._dcBuffer = wx.Bitmap(*self.canvasDimensions)
        self.Render()
        self.Bind(wx.EVT_KEY_DOWN, self.OnKeyPress)
        self.Bind(wx.EVT_PAINT,
                  lambda evt: wx.BufferedPaintDC(self, self._dcBuffer, wx.BUFFER_VIRTUAL_AREA)
                  )

    def Render(self):
        """Render all nodes and their connection in depth order."""
        cdc = wx.ClientDC(self)
        self.PrepareDC(cdc)
        dc = wx.BufferedDC(cdc, self._dcBuffer)
        dc.Clear()
        gc = wx.GraphicsContext.Create(dc)

        gc.SetPen(wx.Pen('#000000', 2, wx.SOLID))
        #gc.DrawRoundedRectangle(12, 34, 56, 78, 10)
        #gc.DrawRoundedRectangle(112, 134, 156, 178, 10)
        #gc.DrawText("(text)", 22, 44)


    def OnKeyPress(self, evt):
        e = evt.GetKeyCode()
        print("OnKeyDown event %s" % (evt))
        if e == 61:
            self.XZoom *= 1.2

        if e == 45:
            self.XZoom /= 1.2
        if e == 48:
            self.YZoom *= 1.2
        if e == 57:
            self.YZoom /= 1.2
        if e == 32:
            self.XZoom = 20
            self.YZoom = 20
        x = self.GetParent().GetChildren()
        y = x[0]
        y.drow_pic(evt)


    def print_pic(self, ls, rs, list):
        cdc = wx.ClientDC(self)
        self.PrepareDC(cdc)
        dc = wx.BufferedDC(cdc, self._dcBuffer)
        dc.Clear()
        gc = wx.GraphicsContext.Create(dc)

        gc.SetPen(wx.Pen('#000000', 2, wx.SOLID))
        x = 20
        xmax = 0
        for i in range(len(list)//4):
            xmax += int(list[1 + i*4])
        xmax *= self.XZoom
        xmax += 20
        if ls:
            gc.StrokeLine(20, 0, 20, 450)
            for i in range(20, 450, 20):
                gc.StrokeLine(20, i, 0, i-20)

        for i in range(len(list)//4):
            for j in range(i):
                x += int(list[1 + j*4]) * self.XZoom

            ed = int(list[1 + i*4]) / 10
            if (ed/self.XZoom > 10):
                ed = 10*self.XZoom

            if ( int(list[i*4])!= 0):
                if ( int(list[i*4]) > 0 ):
                    gc.StrokeLine(x, 202, x + 2*ed * self.XZoom, 202)
                    gc.StrokeLine(x, 198, x + 2*ed * self.XZoom, 198)
                    gc.StrokeLine(x + ed* self.XZoom, 193, x + 2*ed * self.XZoom, 200)
                    gc.StrokeLine(x + ed* self.XZoom , 207, x + 2*ed * self.XZoom, 200)
                if ( int(list[i*4]) < 0 ):
                    gc.StrokeLine(x - 2*ed* self.XZoom, 202, x, 202)
                    gc.StrokeLine(x - 2*ed* self.XZoom, 198, x , 198)
                    gc.StrokeLine(x - 1*ed * self.XZoom, 193, x - 2*ed * self.XZoom,200 )
                    gc.StrokeLine(x - 1*ed * self.XZoom, 207, x - 2*ed * self.XZoom, 200)
            if ((i+1)*4 == len(list)-1):
                if ( int(list[len(list)-1]) > 0 ):
                    gc.StrokeLine(x + int(int(list[1 + i*4]) * self.XZoom), 202, x + int(int(list[1 + i*4]) * self.XZoom) + 2*ed * self.XZoom, 202)
                    gc.StrokeLine(x + int(int(list[1 + i*4]) * self.XZoom), 198, x + int(int(list[1 + i*4]) * self.XZoom) + 2*ed * self.XZoom, 198)
                    gc.StrokeLine(x + int(int(list[1 + i*4]) * self.XZoom) + ed* self.XZoom, 193, x + int(int(list[1 + i*4]) * self.XZoom) + 2*ed * self.XZoom, 200)
                    gc.StrokeLine(x + int(int(list[1 + i*4]) * self.XZoom) + ed* self.XZoom , 207, x + int(int(list[1 + i*4]) * self.XZoom) + 2*ed * self.XZoom, 200)
                if ( int(list[len(list)-1]) < 0 ):
                    gc.StrokeLine(x + int(int(list[1 + i*4]) * self.XZoom) - 2*ed* self.XZoom, 202, x + int(int(list[1 + i*4]) * self.XZoom), 202)
                    gc.StrokeLine(x + int(int(list[1 + i*4]) * self.XZoom) - 2*ed* self.XZoom, 198, x + int(int(list[1 + i*4]) * self.XZoom) , 198)
                    gc.StrokeLine(x + int(int(list[1 + i*4]) * self.XZoom) - 1*ed * self.XZoom, 193, x + int(int(list[1 + i*4]) * self.XZoom) - 2*ed * self.XZoom,200 )
                    gc.StrokeLine(x + int(int(list[1 + i*4]) * self.XZoom) - 1*ed * self.XZoom, 207, x + int(int(list[1 + i*4]) * self.XZoom) - 2*ed * self.XZoom, 200)

            gc.StrokeLine(x, 200, int(x + int(list[1 + i * 4]) * self.XZoom), 200)
            if (int(list [i*4+3]) != 0):
                gc.StrokeLine( x, 200 , int(x+int(list[1 + i*4]) * self.XZoom), 200)
                for j in range(int(x + 20), int(x+int(list[1 + i*4]) * self.XZoom), 15):
                    if (int(list[i * 4 + 3]) < 0):
                        gc.StrokeLine(j, 195, j - 10, 200)
                        gc.StrokeLine(j, 205, j - 10, 200)
                    if (int(list [i*4+3]) > 0):
                        gc.StrokeLine(j-10, 195, j, 200)
                        gc.StrokeLine(j-10, 205, j, 200)

            gc.DrawRoundedRectangle(x, 200 - int(list[2 + i*4])*self.YZoom//2, int(int(list[1 + i*4]) * self.XZoom), self.YZoom*int(list[2 + i*4]), 10)
            #gc.DrawRoundedRectangle(x, 200 - int(list[2 + i*4])//2*self.YZoom/ed, int(int(list[1 + i*4]) * self.XZoom), self.YZoom*int(list[2 + i*4])/ed, 10)
            if rs:
                gc.StrokeLine(xmax, 0, xmax, 450)
                for ii in range(20, 450, 20):
                    gc.StrokeLine(xmax, ii, xmax+20, ii-20)
            for jj in range(i):
                x -= int(list[1 + jj*4]) * self.XZoom
        return

if __name__ == "__main__":
    app = wx.App(False)
    frame = MyFrame()

    app.MainLoop()
