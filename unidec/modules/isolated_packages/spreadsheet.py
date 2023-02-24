import wx.grid as gridlib
import wx
import unidec.modules.miscwindows as misc
import pandas as pd
import os
import wx.lib.mixins.inspection
import wx.html
import webbrowser
from unidec.modules.matchtools import file_to_df


# Taken from https://stackoverflow.com/questions/28509629/
# work-with-ctrl-c-and-ctrl-v-to-copy-and-paste-into-a-wx-grid-in-wxpython
class MyGrid(wx.grid.Grid):
    def __init__(self, parent, panel):
        wx.grid.Grid.__init__(self, panel, wx.ID_ANY, wx.DefaultPosition, wx.DefaultSize, 0)
        self.Bind(wx.EVT_KEY_DOWN, self.on_key)
        self.Bind(wx.grid.EVT_GRID_CELL_CHANGING, self.on_change)
        self.Bind(wx.grid.EVT_GRID_LABEL_RIGHT_CLICK, self.on_label_right_click)
        self.Bind(wx.grid.EVT_GRID_CELL_RIGHT_CLICK, self.on_cell_right_click)
        self.Bind(wx.grid.EVT_GRID_CELL_LEFT_DCLICK, self.open_link)
        self.selected_rows = []
        self.selected_cols = []
        self.history = []
        self.parent = parent

    def set_col_headers(self, headers):
        ncol = self.GetNumberCols()
        if ncol < len(headers):
            self.InsertCols(ncol, len(headers) - ncol)
        for i, h in enumerate(headers):
            self.SetColLabelValue(i, h)

    def get_col_headers(self):
        return [self.GetColLabelValue(col) for col in range(self.GetNumberCols())]

    def get_table(self):
        for row in range(self.GetNumberRows()):
            result = {}
            for col, header in enumerate(self.get_col_headers()):
                result[header] = self.GetCellValue(row, col)
            yield result

    def get_list(self):
        data = []
        for row in range(self.GetNumberRows()):
            row_list = []
            for col, header in enumerate(self.get_col_headers()):
                result = self.GetCellValue(row, col)
                row_list.append(result)
            data.append(row_list)
        return data

    def remove_empty(self, col_header=None):
        delete_rows = []
        for row in range(self.GetNumberRows()):
            all_bad = True
            for col, header in enumerate(self.get_col_headers()):
                if header == col_header:
                    result = self.GetCellValue(row, col)
                    if result == "":
                        delete_rows.append(row)
                elif col_header is None:
                    result = self.GetCellValue(row, col)
                    if result != "":
                        all_bad = False
            if all_bad and col_header is None:
                delete_rows.append(row)

        for row in delete_rows[::-1]:
            self.DeleteRows(row)

    def fill_empty(self, fillvalue="0"):
        for row in range(self.GetNumberRows()):
            for col, header in enumerate(self.get_col_headers()):
                result = self.GetCellValue(row, col)
                if result == "":
                    self.SetCellValue(row, col, fillvalue)

    def get_df(self):
        self.fill_empty()
        data = self.get_list()
        df = pd.DataFrame(data, columns=self.get_col_headers())
        return df

    def set_df(self, df):
        headers = df.columns.tolist()
        self.set_col_headers(headers)
        nrow = self.GetNumberRows()
        if nrow < len(df):
            self.InsertRows(nrow, len(df) - nrow)
        for row in range(len(df)):
            for col, header in enumerate(self.get_col_headers()):
                value = str(df.iloc[row, col])
                self.SetCellValue(row, col, value)
                if ".html" in value or ".htm" in value or "http" in value:
                    self.SetCellTextColour(row, col, wx.BLUE)
                    self.SetCellFont(row, col,
                                     wx.Font(10, wx.FONTFAMILY_DEFAULT, wx.FONTSTYLE_NORMAL, wx.FONTWEIGHT_BOLD,
                                             wx.TEXT_ATTR_FONT_UNDERLINE))
                    self.SetReadOnly(row, col, True)
        self.autoscale_column_width()

    def clear_all(self):
        for row in range(self.GetNumberRows()):
            for col, header in enumerate(self.get_col_headers()):
                self.SetCellValue(row, col, str(""))

    def add_rows(self, event=None):
        for row in self.selected_rows:
            self.InsertRows(row)
        self.add_history({"type": "add_rows", "rows": self.selected_rows})

    def set_row_values(self, row, **kwargs):
        found_row = False
        for col, header in enumerate(self.get_col_headers()):
            for k in kwargs.keys():
                if k == header:
                    if self.GetCellValue(row, col) == "":
                        self.SetCellValue(row, col, kwargs[k])
                        found_row = True
        return found_row

    def add_line(self, **kwargs):
        found_row = False
        row_count = 0
        for row in range(self.GetNumberRows()):
            if not found_row:
                row_count = row
                found_row = self.set_row_values(row, **kwargs)

        if not found_row:
            row_count += 1
            self.InsertRows(row_count)
            self.set_row_values(row_count, **kwargs)

        return row_count

    def delete_rows(self, event):
        self.cut(event)
        rows = []
        for row in reversed(self.selected_rows):
            rows.append((
                row,
                {  # More attributes can be added
                    "label": self.GetRowLabelValue(row),
                    "size": self.GetRowSize(row)
                }
            ))
            self.DeleteRows(row)
        self.add_history({"type": "delete_rows", "rows": rows})

    def add_cols(self, event):
        for col in self.selected_cols:
            self.InsertCols(col + 1)
            try:
                self.parent.on_add_col(event, col + 1)
            except Exception as e:
                print(e)
        self.add_history({"type": "add_cols", "cols": self.selected_cols})

    def rename_col(self, event):
        dlg = misc.SingleInputDialog(self)
        dlg.initialize_interface(title="Rename Column", message="Column Name:", defaultvalue="")
        dlg.ShowModal()
        try:
            name = str(dlg.value)
        except (ValueError, TypeError, AttributeError):
            return
        for col in self.selected_cols:
            self.SetColLabelValue(col, name)
            pass

    def delete_cols(self, event):
        self.delete(event)
        cols = []
        for col in reversed(self.selected_cols):
            cols.append((
                col,
                {  # More attributes can be added
                    "label": self.GetColLabelValue(col),
                    "size": self.GetColSize(col)
                }
            ))
            self.DeleteCols(col)
        self.add_history({"type": "delete_cols", "cols": cols})

    def on_cell_right_click(self, event):
        menus = [(wx.NewId(), "Cut", self.cut),
                 (wx.NewId(), "Copy", self.copy),
                 (wx.NewId(), "Paste", self.paste)]
        popup_menu = wx.Menu()
        for menu in menus:
            if menu is None:
                popup_menu.AppendSeparator()
                continue
            popup_menu.Append(menu[0], menu[1])
            self.Bind(wx.EVT_MENU, menu[2], id=menu[0])

        self.PopupMenu(popup_menu, event.GetPosition())
        popup_menu.Destroy()
        return

    def on_label_right_click(self, event):
        menus = [(wx.NewId(), "Cut", self.cut),
                 (wx.NewId(), "Copy", self.copy),
                 (wx.NewId(), "Paste", self.paste),
                 None]

        # Select if right clicked row or column is not in selection
        if event.GetRow() > -1:
            if not self.IsInSelection(row=event.GetRow(), col=1):
                self.SelectRow(event.GetRow())
            self.selected_rows = self.GetSelectedRows()
            menus += [(wx.NewId(), "Add row", self.add_rows)]
            menus += [(wx.NewId(), "Delete row", self.delete_rows)]
        elif event.GetCol() > -1:
            if not self.IsInSelection(row=1, col=event.GetCol()):
                self.SelectCol(event.GetCol())
            self.selected_cols = self.GetSelectedCols()
            menus += [(wx.NewId(), "Rename column", self.rename_col)]
            menus += [(wx.NewId(), "Add column", self.add_cols)]
            menus += [None]
            menus += [(wx.NewId(), "Delete column", self.delete_cols)]

        else:
            return

        popup_menu = wx.Menu()
        for menu in menus:
            if menu is None:
                popup_menu.AppendSeparator()
                continue
            popup_menu.Append(menu[0], menu[1])
            self.Bind(wx.EVT_MENU, menu[2], id=menu[0])

        self.PopupMenu(popup_menu, event.GetPosition())
        popup_menu.Destroy()
        return

    def on_change(self, event):
        cell = event.GetEventObject()
        row = cell.GetGridCursorRow()
        col = cell.GetGridCursorCol()
        attribute = {"value": self.GetCellValue(row, col)}
        self.add_history({"type": "change", "cells": [(row, col, attribute)]})

    def add_history(self, change):
        self.history.append(change)

    def undo(self):
        if not len(self.history):
            return

        action = self.history.pop()
        if action["type"] == "change" or action["type"] == "delete":
            for row, col, attribute in action["cells"]:
                self.SetCellValue(row, col, attribute["value"])
                if action["type"] == "delete":
                    self.SetCellAlignment(row, col, *attribute["alignment"])  # *attribute["alignment"] > horiz, vert

        elif action["type"] == "delete_rows":
            for row, attribute in reversed(action["rows"]):
                self.InsertRows(row)
                self.SetRowLabelValue(row, attribute["label"])
                self.SetRowSize(row, attribute["size"])

        elif action["type"] == "delete_cols":
            for col, attribute in reversed(action["cols"]):
                self.InsertCols(col)
                self.SetColLabelValue(col, attribute["label"])
                self.SetColSize(col, attribute["size"])

        elif action["type"] == "add_rows":
            for row in reversed(action["rows"]):
                self.DeleteRows(row)

        elif action["type"] == "add_cols":
            for col in reversed(action["cols"]):
                self.DeleteCols(col)
        else:
            return

    def on_key(self, event):
        """
        Handles all key events.
        """
        # print(event.GetKeyCode())
        # Ctrl+C or Ctrl+Insert
        if event.ControlDown() and event.GetKeyCode() in [67, 322]:
            self.copy(event)

        # Ctrl+V
        elif event.ControlDown() and event.GetKeyCode() == 86:
            self.paste(event)

        # DEL
        elif event.GetKeyCode() == 127:
            self.delete(event)

        # Ctrl+A
        elif event.ControlDown() and event.GetKeyCode() == 65:
            self.SelectAll()

        # Ctrl+Z
        elif event.ControlDown() and event.GetKeyCode() == 90:
            self.undo()

        # Ctrl+X
        elif event.ControlDown() and event.GetKeyCode() == 88:
            # Call delete method
            self.cut(event)

        # Ctrl+V or Shift + Insert
        elif (event.ControlDown() and event.GetKeyCode() == 67) \
                or (event.ShiftDown() and event.GetKeyCode() == 322):
            self.paste(event)

        else:
            event.Skip()

    def get_selection(self):
        """
        Returns selected range's start_row, start_col, end_row, end_col
        If there is no selection, returns selected cell's start_row=end_row, start_col=end_col
        """
        if not len(self.GetSelectionBlockTopLeft()):
            selected_columns = self.GetSelectedCols()
            selected_rows = self.GetSelectedRows()
            if selected_columns:
                start_col = selected_columns[0]
                end_col = selected_columns[-1]
                start_row = 0
                end_row = self.GetNumberRows() - 1
            elif selected_rows:
                start_row = selected_rows[0]
                end_row = selected_rows[-1]
                start_col = 0
                end_col = self.GetNumberCols() - 1
            else:
                start_row = end_row = self.GetGridCursorRow()
                start_col = end_col = self.GetGridCursorCol()
        elif len(self.GetSelectionBlockTopLeft()) > 1:
            wx.MessageBox("Multiple selections are not supported", "Warning")
            return []
        else:
            start_row, start_col = self.GetSelectionBlockTopLeft()[0]
            end_row, end_col = self.GetSelectionBlockBottomRight()[0]

        return [start_row, start_col, end_row, end_col]

    def get_selected_cells(self):
        # returns a list of selected cells
        selection = self.get_selection()
        if not selection:
            return

        start_row, start_col, end_row, end_col = selection
        for row in range(start_row, end_row + 1):
            for col in range(start_col, end_col + 1):
                yield [row, col]

    def get_selected_rows(self):
        # returns a list of selected rows
        selection = self.get_selection()
        if not selection:
            return

        start_row, start_col, end_row, end_col = selection
        for row in range(start_row, end_row + 1):
            yield row

    def copy(self, event):
        """
        Copies range of selected cells to clipboard.
        """

        selection = self.get_selection()
        if not selection:
            return []
        start_row, start_col, end_row, end_col = selection

        data = u''

        rows = range(start_row, end_row + 1)
        for row in rows:
            columns = range(start_col, end_col + 1)
            for idx, column in enumerate(columns, 1):
                if idx == len(columns):
                    # if we are at the last cell of the row, add new line instead
                    data += self.GetCellValue(row, column) + "\n"
                else:
                    data += self.GetCellValue(row, column) + "\t"

        text_data_object = wx.TextDataObject()
        text_data_object.SetText(data)

        if wx.TheClipboard.Open():
            wx.TheClipboard.SetData(text_data_object)
            wx.TheClipboard.Close()
        else:
            wx.MessageBox("Can't open the clipboard", "Warning")

    def paste(self, event):
        if not wx.TheClipboard.Open():
            wx.MessageBox("Can't open the clipboard", "Warning")
            return False

        clipboard = wx.TextDataObject()
        wx.TheClipboard.GetData(clipboard)
        wx.TheClipboard.Close()
        data = clipboard.GetText()
        if data[-1] == "\n":
            data = data[:-1]

        try:
            cells = self.get_selected_cells()
            cell = next(cells)
        except StopIteration:
            return False

        start_row = end_row = cell[0]
        start_col = end_col = cell[1]
        max_row = self.GetNumberRows()
        max_col = self.GetNumberCols()

        history = []
        out_of_range = False

        for row, line in enumerate(data.split("\n")):
            target_row = start_row + row
            if not (0 <= target_row < max_row):
                out_of_range = True
                break

            if target_row > end_row:
                end_row = target_row

            for col, value in enumerate(line.split("\t")):
                target_col = start_col + col
                if not (0 <= target_col < max_col):
                    out_of_range = True
                    break

                if target_col > end_col:
                    end_col = target_col

                # save previous value of the cell for undo
                history.append([target_row, target_col, {"value": self.GetCellValue(target_row, target_col)}])

                self.SetCellValue(target_row, target_col, value)

        self.SelectBlock(start_row, start_col, end_row, end_col)  # select pasted range
        if out_of_range:
            wx.MessageBox("Pasted data is out of Grid range", "Warning")

        self.add_history({"type": "change", "cells": history})

    def delete(self, event):
        cells = []
        for row, col in self.get_selected_cells():
            attributes = {
                "value": self.GetCellValue(row, col),
                "alignment": self.GetCellAlignment(row, col)
            }
            cells.append((row, col, attributes))
            self.SetCellValue(row, col, "")

        self.add_history({"type": "delete", "cells": cells})

    def cut(self, event):
        self.copy(event)
        self.delete(event)

    def open_link(self, event):
        row = event.GetRow()
        col = event.GetCol()
        value = self.GetCellValue(row, col)
        if ".html" in value or "http" in value:
            print("Opening link:", value)
            webbrowser.open(value)
        elif "Open In UniDec" in value:
            print("Opening UniDec")
            self.parent.open_unidec(row)
        else:
            event.Skip()

    def save_file(self, filename=None):
        """
        Saves the grid to a file.
        :param filename: filename to save to
        :return: None
        """
        if filename is not None:
            df = self.get_df()
            if filename.endswith(".csv"):
                df.to_csv(filename, index=False)
            elif filename.endswith(".xlsx"):
                df.to_excel(filename, index=False)
            elif filename.endswith(".txt"):
                df.to_csv(filename, sep="\t", index=False)
            else:
                df.to_csv(filename, index=False)
            print("Saved Spreadsheet to: ", filename)

    def autoscale_column_width(self):
        # self.SetDefaultRenderer(CutomGridCellAutoWrapStringRenderer())
        self.AutoSizeColumns()
        for col in range(self.GetNumberCols()):
            if self.GetColSize(col) > 250:
                self.SetColSize(col, 250)
        self.AutoSizeRows()


class SpreadsheetPanel:
    def __init__(self, parent, panel, nrows=0, ncolumns=0):
        self.ss = MyGrid(parent, panel)
        self.nrows = nrows
        self.ncolumns = ncolumns
        self.ss.CreateGrid(nrows, ncolumns)
        self.parent = parent

    def set_col_labels(self, labels):
        for i, l in enumerate(labels):
            self.ss.SetColLabelValue(i, l)

    def add_col(self, label=""):
        self.ss.AppendCols(1)


class SpreadsheetFrame(wx.Frame):
    """"""

    def __init__(self, nrows, ncolumns, title=""):
        """Constructor"""
        wx.Frame.__init__(self, parent=None, title=title)
        panel = wx.Panel(self)

        self.ss = SpreadsheetPanel(self, panel, nrows, ncolumns).ss

        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(self.ss, 1, wx.EXPAND)
        panel.SetSizer(sizer)
        self.Show()


from wx.lib import wordwrap


class CutomGridCellAutoWrapStringRenderer(wx.grid.GridCellRenderer):
    def __init__(self):
        wx.grid.GridCellRenderer.__init__(self)

    def Draw(self, grid, attr, dc, rect, row, col, isSelected):
        text = grid.GetCellValue(row, col)
        dc.SetFont(attr.GetFont())
        text = wordwrap.wordwrap(text, grid.GetColSize(col), dc, breakLongWords=True)
        hAlign, vAlign = attr.GetAlignment()
        if isSelected:
            bg = grid.GetSelectionBackground()
            fg = grid.GetSelectionForeground()
        else:
            bg = attr.GetBackgroundColour()
            fg = attr.GetTextColour()
        dc.SetTextBackground(bg)
        dc.SetTextForeground(fg)
        dc.SetBrush(wx.Brush(bg, wx.SOLID))
        dc.SetPen(wx.TRANSPARENT_PEN)
        dc.DrawRectangle(rect)
        grid.DrawTextRectangle(dc, text, rect, hAlign, vAlign)

    def GetBestSize(self, grid, attr, dc, row, col):
        text = grid.GetCellValue(row, col)
        dc.SetFont(attr.GetFont())
        text = wordwrap.wordwrap(text, grid.GetColSize(col), dc, breakLongWords=False)
        w, h, lineHeight = dc.GetMultiLineTextExtent(text)
        return wx.Size(w, h)

    def Clone(self):
        return CutomGridCellAutoWrapStringRenderer()


'''
class HtmlRenderer(wx.grid.GridCellRenderer):
    def __init__(self, html_str):
        wx.grid.GridCellRenderer.__init__(self)
        self.html_str = html_str
        self.html = wx.html.HtmlDCRenderer()

    def Draw(self, grid, attr, dc, rect, row, col, isSelected):
        self.html.SetDC(dc)
        self.html.SetSize(100, 20)
        self.html.SetHtmlText(self.html_str)

        dc.SetBackgroundMode(wx.SOLID)
        if isSelected:
            dc.SetBrush(wx.Brush(wx.BLUE, wx.SOLID))
            dc.SetPen(wx.Pen(wx.BLUE, 1, wx.SOLID))
        else:
            dc.SetBrush(wx.Brush(wx.WHITE, wx.SOLID))
            dc.SetPen(wx.Pen(wx.WHITE, 1, wx.SOLID))
        dc.DrawRectangle(rect)


        width, height = html.GetTotalWidth(), html.GetTotalHeight()
        if width > rect.width - 2:
            width = rect.width - 2
        if height > rect.height - 2:
            height = rect.height - 2

        self.html.Render(5, -2)
        # dc.Blit(rect.x + 1, rect.y + 1, width, height, html, 0, 0, wx.COPY, True)

    def on_open_link(self, event):
        print("Test")
   
# renderer = HtmlRenderer("<a href=\"http:\\www.google.com\">test</a>")
# frame.ss.SetCellRenderer(0, 0, renderer)
# frame.ss.SetReadOnly(0, 0)
# frame.ss.Bind(wx.html.EVT_HTML_CELL_CLICKED, renderer.on_open_link)     
'''


class MyApp(wx.App, wx.lib.mixins.inspection.InspectionMixin):
    def OnInit(self):
        self.Init()
        self.frame = SpreadsheetFrame(12, 8)
        self.frame.Show()
        self.SetTopWindow(self.frame)
        return True


if __name__ == "__main__":
    # app = wx.App()
    app = MyApp(redirect=False)
    frame = app.frame
    frame.ss.SetCellRenderer(0, 0, CutomGridCellAutoWrapStringRenderer())
    frame.ss.SetCellValue(0, 0, "http://www.google.com/test/test/test/test")
    frame.ss.SetReadOnly(0, 0)
    # frame = SpreadsheetFrame(12, 8)
    # path = "C:\\Data\\Luis Genentech\\Merged Glycan List.csv"
    # df = pd.read_csv(path)
    # frame.ss.set_df(df)
    # print(df)

    app.MainLoop()
