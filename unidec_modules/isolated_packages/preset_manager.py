import wx
import os
import numpy as np


def make_preset_menu(toppath=None, exclude_dir=None, topi=1500, ext="", exclude_ext=None, exclude_dir_list=[]):
    if toppath is None:
        toppath = "C:\\Python\\UniDec3\\unidec_bin\\Presets\\"
    custommenu = wx.Menu()
    opmenu = custommenu
    masterd = []
    i = 0

    if exclude_dir is None:
        exclude_dir = "SOME STRING THAT COULD NEVER EXIST IN A NORMAL FOLDER. I realize this is weird but it works..."
    if exclude_ext is None:
        exclude_ext = "SOME STRING THAT COULD NEVER EXIST IN A NORMAL Extension. "
    if ext is None:
        ext = ""

    never = "UniDec_Figures_and_Files"

    try:
        for p in os.listdir(toppath):
            path = os.path.join(toppath, p)
            if os.path.isfile(path) and ext in os.path.splitext(path)[1] and exclude_ext not in os.path.splitext(path)[
                1]:
                mitem = opmenu.Append(topi + i, p)
                masterd.append([topi + i, path, mitem])
                i += 1
                # print(p)
            if os.path.isdir(path) and exclude_dir not in path and never not in path \
                    and os.path.split(path)[1] not in exclude_dir_list:
                opmenu2 = wx.Menu()
                opmenu.AppendSubMenu(opmenu2, p)
                for p in os.listdir(path):
                    path2 = os.path.join(path, p)
                    if os.path.isfile(path2) and ext in os.path.splitext(path2)[1] and exclude_ext not in \
                            os.path.splitext(path2)[1]:
                        mitem = opmenu2.Append(topi + i, p)
                        masterd.append([topi + i, path2, mitem])
                        i += 1
                        # print(p)
                    if os.path.isdir(
                            path2) and exclude_dir not in path2 and never not in path2\
                            and os.path.split(path2)[1] not in exclude_dir_list:
                        opmenu3 = wx.Menu()
                        opmenu2.AppendSubMenu(opmenu3, p)
                        for p in os.listdir(path2):
                            path3 = os.path.join(path2, p)
                            if os.path.isfile(path3) and ext in os.path.splitext(path3)[1] and exclude_ext not in \
                                    os.path.splitext(path3)[1]:
                                mitem = opmenu3.Append(topi + i, p)
                                masterd.append([topi + i, path3, mitem])
                                i += 1
                                # print(p)
                            if os.path.isdir(
                                    path3) and exclude_dir not in path3 and never not in path3 \
                                    and os.path.split(path3)[1] not in exclude_dir_list:
                                opmenu4 = wx.Menu()
                                opmenu3.AppendSubMenu(opmenu4, p)
                                for p in os.listdir(path3):
                                    path4 = os.path.join(path3, p)
                                    if os.path.isfile(path4) and ext in os.path.splitext(path4)[
                                        1] and exclude_ext not in os.path.splitext(path4)[1]:
                                        mitem = opmenu4.Append(topi + i, p)
                                        masterd.append([topi + i, path4, mitem])
                                        i += 1
                                        # print(p)
                                    if os.path.isdir(path4):
                                        print("Error: Reached Max Recursion Depth in Presents Folder")
    except Exception as e:
        print("Problem with preset manager: ", e)
    masterd = np.array(masterd)
    # print(masterd)
    return custommenu, masterd


class TestWindow(wx.Frame):
    def __init__(self, parent=None):
        wx.Frame.__init__(self, parent, title="Test")
        menu_bar = wx.MenuBar()
        filemenu = wx.Menu()

        custommenu, self.masterd = make_preset_menu(toppath="C:\\Python\\UniDec3\\unidec_bin\\Example Data",
                                                    exclude_dir="_unidecfiles")
        filemenu.AppendSubMenu(custommenu, "Custom")
        for i, path, item in self.masterd:
            self.Bind(wx.EVT_MENU, self.on_defaults, item)

        menu_bar.Append(filemenu, "&File")
        self.SetMenuBar(menu_bar)
        self.Centre()
        # self.MakeModal(True)

        self.Show(True)

    def on_defaults(self, e):
        try:
            nid = e.GetId()
            ids = self.masterd[:, 0].astype(np.float)
            pos = np.argmin(np.abs(ids - nid))
            path = self.masterd[pos, 1]
            print("Opening Path:", path)
        except Exception as e:
            print(e)


if __name__ == "__main__":
    app = wx.App(False)
    frame = TestWindow()
    app.MainLoop()
