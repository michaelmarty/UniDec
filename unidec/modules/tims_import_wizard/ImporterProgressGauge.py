# threaded operations used by ionMS

import wx

from pubsub import pub


dlg = None
counter = 1
maximum = 1
pub_message = None


# this is just used in the main thread
def progress_dialog(parent_window, title, message, num_files, pub_mess=None):
    global dlg, counter, maximum, pub_message
    maximum = num_files
    counter = 1
    pub_message = pub_mess

    pub.subscribe(on_callback, 'RAW DATA ADDED TO MODEL')

    # create a file loading progress dialog
    dlg = wx.ProgressDialog(
            title,
            message,
            maximum=maximum,
            parent=parent_window,
            style=wx.PD_APP_MODAL
    )

def on_callback():
    global dlg, counter, maximum, pub_message
    if counter == maximum:
        #message = "Finished Import"
        #(keepGoing, skip) = dlg.Update(counter, message)
        #print "Done"
        #wx.MilliSleep(200)
        dlg.Destroy()
        #wx.MilliSleep(200) # if this isn't here, then pythonw.exe will crash (python.exe is fine though), was sometimes crashing with just 50 ms
        # actually sometimes just crashes regardless, maybe this isn't the cause

        # check on dlg's parent_window (not sure what attribute)
        # to decide what kind of message to send
        # is a bit hacky, but..
        #if pub_message is None:
        #    pub.sendMessage('DATA ADDED')
        #else:
        #    pub.sendMessage(pub_message)


    else:
        #print "Keep Going"
        message = "Importing file %s of %s, please wait..." % (counter + 1, maximum)
        (keepGoing, skip) = dlg.Update(counter, message)
        counter += 1
