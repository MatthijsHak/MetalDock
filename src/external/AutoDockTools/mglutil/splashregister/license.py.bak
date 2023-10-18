# $Header: /mnt/raid/services/cvs/python/packages/share1.5/mglutil/splashregister/license.py,v 1.15 2010/08/19 20:14:31 sargis Exp $
# $Id: license.py,v 1.15 2010/08/19 20:14:31 sargis Exp $
#

import Tkinter

tk_root = Tkinter.Tk()
tk_root.title("Commercial Usage")
txt = """
 The software component for computing molecular surfaces (MSMS) 
 is not free for commercial usage. If you plan to use MSMS for commercial 
 research please contact sanner@scripps.edu

 Some software components such at the volume rendering and 
 isocontouring were developed at UT Austin.

 If you publish scientific results generated using this software 
 please cite the appropriate software components.
 A list of papers is provided under Help -> Citation Information 
 menu in PMV and ADT. 
"""
Tkinter.Label(tk_root, text=txt, justify=Tkinter.LEFT).pack()
Tkinter.Button(tk_root, text="OK", command=tk_root.quit).pack()

# The following is to position the window in the center of the screen:
# A hack to get the window size. Temporarily hide the
# window to avoid update_idletasks() drawing the window in the wrong
# position.
tk_root.withdraw()
tk_root.update_idletasks()  # Update "requested size" from geometry manager

x = (tk_root.winfo_screenwidth() - tk_root.winfo_reqwidth()) / 2
y = (tk_root.winfo_screenheight() - tk_root.winfo_reqheight()) / 2
tk_root.geometry("+%d+%d" % (x, y))

# This seems to draw the window frame immediately, so only call deiconify()
# after setting correct window position
tk_root.deiconify()

#tk_root.eval('tk::PlaceWindow %s center' % tk_root.winfo_pathname(tk_root.winfo_id())) # this  also seems to work (tested on Linux)
tk_root.mainloop()
