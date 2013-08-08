#!/usr/bin/env python
# template user1.py
# To add some user widgets in yao interfaces
# 1) modify user*.glade:
# - implement button/widgets in user*.glade & edit for callbacks.
# 2) in user*.py:
# - change label
# - add callback functions to your buttons, etc.
# first panel is user1
# second panel is user2

import gobject
import gtk
import gtk.glade
import os
import string
import sys


class user1:

   def __init__(self, path='.', parent=None, py2yo=None):

      self.py2yo=py2yo
      self.gladefile = 'user1.glade'
      self.glade = gtk.glade.XML(os.path.join(path,self.gladefile), root='top')
      self.top = self.glade.get_widget('top')
      self.glade.signal_autoconnect(self)

      if parent:
         parent.foreach(parent.remove)
         parent.add(self.top)

      parent.show()

   def label(self):
      return 'my button'
