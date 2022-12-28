import wx
import os
from natsort import natsorted

class TreeCtrl(wx.TreeCtrl):
    def __init__(self, parent, id, pos, size, style):
        wx.TreeCtrl.__init__(self, parent, id, pos, size, style)

#---------------------------------------------------------------------------

class TreeCtrlPanel(wx.Panel):
    def __init__(self, parent, link):
        # Use the WANTS_CHARS style so the panel doesn't eat the Return key.
        wx.Panel.__init__(self, parent, -1, style=wx.WANTS_CHARS)
        self.Bind(wx.EVT_SIZE, self.OnSize)

        self.link = link

        self.tree = TreeCtrl(self, wx.NewIdRef(),
                             wx.DefaultPosition,
                             (300, 300),
                             wx.TR_DEFAULT_STYLE | wx.TR_MULTIPLE)

        isz = (16, 16)
        il = wx.ImageList(isz[0], isz[1])
        self.fldridx = il.Add(wx.ArtProvider.GetBitmap(wx.ART_FOLDER, wx.ART_OTHER, isz))
        self.fldropenidx = il.Add(wx.ArtProvider.GetBitmap(wx.ART_FILE_OPEN, wx.ART_OTHER, isz))
        self.fileidx = il.Add(wx.ArtProvider.GetBitmap(wx.ART_NORMAL_FILE, wx.ART_OTHER, isz))

        self.tree.SetImageList(il)
        self.il = il

        self.path = None
        self.raw = None

        self.Bind(wx.EVT_TREE_SEL_CHANGED, self.on_selected_changed, self.tree)
        self.Bind(wx.EVT_TREE_ITEM_ACTIVATED, self.on_activate, self.tree)
        self.Bind(wx.EVT_TREE_ITEM_EXPANDED, self.on_item_expanded, self.tree)

        self.populate_tree()


    # instantiates the tree, uses path in gui, which has a default value
    def populate_tree(self):
        path = self.link.folder_path.GetValue()
        if os.path.isdir(path):
            if path != self.path:
                self.path = path
                self.tree.DeleteAllItems()
                root = self.add_root(os.path.basename(path), path)
                self.add_children(root, path, 2)
                self.tree.Expand(root)


    def add_root(self, name, data=None):
        """
        Add root folder to tree
            returns root object
        """
        root = self.tree.AddRoot(name)
        self.tree.SetItemData(root, data)
        self.tree.SetItemImage(root, self.fldridx, wx.TreeItemIcon_Normal)
        self.tree.SetItemImage(root, self.fldropenidx, wx.TreeItemIcon_Expanded)
        return root


    def add_children(self, parent, path, depth_limit=1):
        """
        Recursively adds children up to the depth specified
        """
        try:
            for folder in natsorted(os.listdir(path)):
                new_path = os.path.join(path, folder)
                if not os.path.isfile(new_path) and \
                   not folder.startswith('.')  and \
                   not folder.startswith('_'):
                    if os.path.isdir(new_path):
                        child = self.tree.AppendItem(parent, os.path.basename(new_path))
                        self.tree.SetItemData(child, new_path)
                        self.tree.SetItemImage(child, self.fldridx, wx.TreeItemIcon_Normal)
                        self.tree.SetItemImage(child, self.fldropenidx, wx.TreeItemIcon_Expanded)
                        if depth_limit > 0:
                            self.add_children(child, new_path, depth_limit - 1)
        except OSError:
            pass


    def raw_file_info(self, path):
        if not '.raw' in path:
            self.link.desc.Clear()
        else:
            if path != self.raw:
                self.raw = path

                self.link.desc.Clear()

                add_lines = ''

                des_path = os.path.join(path, '_HEADER.TXT')
                if os.path.isfile(des_path):
                    # get description from header
                    f = open(des_path, 'r')#, encoding="utf-8")
                    lines = f.readlines()
                    f.close()
                    for line in lines:
                        if line.startswith('$$ Sample Description:'):
                            add_lines += 'Description from Header File\n'
                            add_lines += '----------------------------------\n'
                            description = line.split('$$ Sample Description:')[1].strip()
                            add_lines += description + '\n\n'

                param_path = os.path.join(path, '_extern.inf')
                if os.path.isfile(param_path):
                    # read in parameter file
                    f = open(param_path, 'r')
                    lines = f.readlines()
                    f.close()
                    add_lines += 'Instrument Parameter File\n'
                    add_lines += '-------------------------------\n'
                    for line in lines:
                        add_lines += line# + '\n'

                if add_lines != '':
                    print("add_lines is blank")

                #add_lines = unicode(add_lines, errors='ignore')
                self.link.desc.SetValue(add_lines)


    def OnSize(self, event):
        w, h = self.GetClientSize()
        self.tree.SetSize(0, 0, w, h)


    def on_selected_changed(self, event):
        self.item = event.GetItem()
        if self.item:
            self.raw_file_info(self.tree.GetItemData(self.item))
            try:
                self.raw_file_info(self.tree.GetItemData(self.item))
            except Exception as e:
                print("Error in file header decoding: ", e)
        event.Skip()


    def on_activate(self, event):
        if self.item:
            if self.tree.IsExpanded(self.item):
                self.tree.CollapseAllChildren(self.item)
            else:
                # remove any existing children:
                self.tree.DeleteChildren(self.item)
                self.add_children(self.item, self.tree.GetItemData(self.item), 1)
                self.tree.Expand(self.item)


    def on_item_expanded(self, event):
        self.item = event.GetItem()
        if self.item:
            # remove any existing children:
            self.tree.DeleteChildren(self.item)
            self.add_children(self.item, self.tree.GetItemData(self.item), 1)
