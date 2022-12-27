# provides file dialog functions
# uses "module as a singleton" approach to remember the default directory

import wx, os
import wx.lib.agw.multidirdialog as MDD
import unidec.modules.unidectools as ud
from pathlib import Path

# setup a default path
default_dir = os.path.abspath(os.path.join(os.getcwd(), os.path.pardir, 'data/'))
if not os.path.exists(default_dir):
    default_dir = os.getcwd()
default_dir = ""


# opens a dialog to save a file
def save_file_dialog(message="Save File", file_types="*.*", default_file=""):
    global default_dir

    dlg = wx.FileDialog(None, message, default_dir, default_file, file_types, wx.FD_SAVE | wx.FD_OVERWRITE_PROMPT)
    path = None
    if dlg.ShowModal() == wx.ID_OK:
        path = ud.smartdecode(dlg.GetPath())
        default_dir = ud.smartdecode(os.path.dirname(path))

    dlg.Destroy()
    return path


# opens a dialog to choose a single file
def open_file_dialog(message="Open File", file_types="*.*", default=None):
    global default_dir
    if default is None:
        dlg = wx.FileDialog(None, message, default_dir, "", file_types, wx.FD_OPEN | wx.FD_FILE_MUST_EXIST)
    else:
        dlg = wx.FileDialog(None, message, "", default, file_types, wx.FD_OPEN | wx.FD_FILE_MUST_EXIST)
    path = None
    if dlg.ShowModal() == wx.ID_OK:
        path = ud.smartdecode(dlg.GetPath())
        default_dir = ud.smartdecode(os.path.dirname(path))

    dlg.Destroy()
    return path


# opens a dialog to choose multiple files
def open_multiple_files_dialog(message="Open Files", file_type="*.*"):
    default_dir = ""

    dlg = wx.FileDialog(None, message, default_dir, "", file_type, wx.FD_OPEN | wx.FD_MULTIPLE | wx.FD_FILE_MUST_EXIST)
    file_names = None
    if dlg.ShowModal() == wx.ID_OK:
        file_names = dlg.GetPaths()
        default_dir = dlg.GetDirectory()
        #if path is not None:
        #    file_names = [ud.smartdecode(os.path.join(os.path.dirname(path.encode('utf-8')), f.encode('utf-8'))) for f
        #                  in dlg.GetFilenames()]
        print("Files:", file_names)
    dlg.Destroy()
    return file_names


# opens a dialog to choose a directory
def open_dir_dialog(message="Select a Folder"):
    global default_dir
    dlg = wx.DirDialog(None, message, default_dir)
    path = None
    if dlg.ShowModal() == wx.ID_OK:
        path = ud.smartdecode(dlg.GetPath())
        default_dir = path

    dlg.Destroy()
    return path


def open_multiple_dir_dialog(message, default):
    if default is None:
        global default_dir
        default = default_dir
    dlg = MDD.MultiDirDialog(None, message=message, defaultPath=default,
                             agwStyle=MDD.DD_MULTIPLE | MDD.DD_DIR_MUST_EXIST)
    #dlg = wx.DirDialog(None, message, default, wx.DD_MULTIPLE)
    dirs = None
    if dlg.ShowModal() == wx.ID_OK:
        dirs = ud.smartdecode(dlg.GetPaths())
        dirs2 = []
        for d in dirs:
            p = Path(d)
            drive = p.parts[0]
            if ")" in drive:
                spot = drive.find(":")
                drive_letter = drive[spot-1]
                d=drive_letter+":\\"
                d2 = os.path.join(*p.parts[1:])
                d = os.path.join(d, d2)
            dirs2.append(d)
        dirs = dirs2
    dlg.Destroy()
    return dirs


def open_single_dir_dialog(message, default):
    global default_dir
    dlg = wx.DirDialog(None, message, default)
    # dlg=MDD.MultiDirDialog(None,message=message,defaultPath=default,agwStyle=MDD.DD_NEW_DIR_BUTTON)
    dir = None
    if dlg.ShowModal() == wx.ID_OK:
        dir = ud.smartdecode(dlg.GetPath())
        default_dir = dir
    dlg.Destroy()
    return dir

# TODO: Simple Figure Dialog App to feed into parsers
