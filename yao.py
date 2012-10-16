#!/usr/bin/env python
#
# yao.py
#
# python functions:
#   - glade callback functions
#   - calls to yorick
#
# This file is part of the yao package, an adaptive optics
# simulation tool.
#
# $Id: yao.py,v 1.9 2010/04/15 02:36:53 frigaut Exp $
#
# Copyright (c) 2002-2007, Francois Rigaut
#
# This program is free software; you can redistribute it and/or  modify it
# under the terms of the GNU General Public License  as  published  by the
# Free Software Foundation; either version 2 of the License,  or  (at your
# option) any later version.
#
# This program is distributed in the hope  that  it  will  be  useful, but
# WITHOUT  ANY   WARRANTY;   without   even   the   implied   warranty  of
# MERCHANTABILITY or  FITNESS  FOR  A  PARTICULAR  PURPOSE.   See  the GNU
# General Public License for more details (to receive a  copy  of  the GNU
# General Public License, write to the Free Software Foundation, Inc., 675
# Mass Ave, Cambridge, MA 02139, USA).
#
# To Do:
# - implement gtk dialog to warn user of danger of re-using the same file name
#   when editing a current par file and saving
# 
# $Log: yao.py,v $
# Revision 1.9  2010/04/15 02:36:53  frigaut
#
#
# final commit to upgrade this repo to yao 4.5.1
#
# Revision 1.8  2007/12/20 13:34:53  frigaut
# - various bug fixes
# - better handlng of default parfile path
# - better handling of options menu (WFS and DM)
#
# Revision 1.7  2007/12/19 19:44:19  frigaut
# - solved a number of bugs and inconsistencies between regular yao call and
#   GUI calls.
# - fixed misregistration for curvature systems
# - change: misregistration entry from the GUI is now in pupil diameter unit,
#   not in subaperture unit!
# - changed default efd in c188-bench.par
#
# Revision 1.6  2007/12/19 15:45:32  frigaut
# - implemented yao.conf which defines the YAO_SAVEPATH directory where
# all temporary files and result files will be saved
# - modified yao.i and aoutil.i to save in YAO_SAVEPATH
# - bumped version to 4.2.0
# - slight changes to GUI (edit conf file)
#
# Revision 1.5  2007/12/19 13:18:59  frigaut
# - explicit message when screens are not present/found
# - more messages in statusbar
# - added statusbar1 (that can hide/show) for strehl status header
#
# Revision 1.4  2007/12/18 19:03:20  frigaut
# - reworked Y_PYTHON and search for yao.py
# - added Y_GLADE and path to yao.glade
# - now removes CVS directories in install of examples and doc
#
# Revision 1.3  2007/12/17 20:21:04  frigaut
# - renamed yaogtk -> yao (and updated Makefile accordingly)
# - gotten rid of usleep() calls in yorick -> python communication. Instead,
# am using a pyk_flush, which send a flush request to python every seconds.
# This is still a hack to turn around the BLOCK bug in python, but at least
# it does not use usleep (cleaner hack?).
# - added debug python <> yorick entry in GUI help menu (set/unset pyk_debug)
#
# Revision 1.2  2007/12/13 00:58:31  frigaut
# added license and header
#
#

import gtk
import gtk.glade
import sys
import gobject
import os, fcntl, errno
from time import *

class yao:
   
   def destroy(self, wdg, data=None):
      self.py2yo('yaopy_quit')
#      gtk.main_quit()
      
   def __init__(self,path2glade,dpi):
      self.path2glade = path2glade
      self.usercmd = 'STOP'
      
      self.glade = gtk.glade.XML(self.path2glade+'/yao.glade') 
      self.window = self.glade.get_widget('window1')
      if (self.window):
         self.window.connect('destroy', self.destroy)

      self.glade.signal_autoconnect(self)

      self.editor = self.glade.get_widget('window2')
      if (self.editor):
         self.editor.connect('delete_event', self.on_editor_close_activate)

      self.editor2 = self.glade.get_widget('window4')
      if (self.editor2):
         self.editor2.connect('delete_event', self.on_editor2_close_activate)

      self.statusbar = self.glade.get_widget('statusbar')
      self.statusbar1 = self.glade.get_widget('statusbar1')
      self.progressbar = self.glade.get_widget('progressbar')

      self.aoinit_popup_menu = self.glade.get_widget('aoinit_popup_menu')
      self.aoinit_popup_button = self.glade.get_widget('aoinit_popup_button')
      self.aoinit_popup_button.connect_object("event", self.on_aoinit_popup_button_clicked, self.aoinit_popup_menu)

      self.aoloop_popup_menu = self.glade.get_widget('aoloop_popup_menu')
      self.aoloop_popup_button = self.glade.get_widget('aoloop_popup_button')
      self.aoloop_popup_button.connect_object("event", self.on_aoloop_popup_button_clicked, self.aoloop_popup_menu)

      # set stdin non blocking, this will prevent readline to block
      fd = sys.stdin.fileno()
      flags = fcntl.fcntl(fd, fcntl.F_GETFL)
      fcntl.fcntl(fd, fcntl.F_SETFL, flags | os.O_NONBLOCK)
      
      # add stdin to the event loop (yorick input pipe by spawn)
      gobject.io_add_watch(sys.stdin,gobject.IO_IN | gobject.IO_HUP,self.yo2py)

      self.wfs_panel_set_sensitivity(0,0)
      self.dm_panel_set_sensitivity(0)
      # self.glade.get_widget('wfs_and_dms').hide()
      # declare and set default (overwritten by gui_update)
      self.yuserdir = "./"
      self.yaopardir = "/"
      self.pyk_debug = 0
      
      # set size of graphic areas:
      dsx = int(635.*dpi/75)+4
      dsy = int(650.*dpi/75)+25
      self.glade.get_widget('drawingarea1').set_size_request(dsx,dsy)
      # self.drawingarea_size_allocate(dpi)
      
      # run
      gtk.main()


   dispflag = 1
   init = 0
   buffer  = gtk.TextBuffer()
   buffer2 = gtk.TextBuffer()

   def on_about_activate(self,wdg):
      dialog = self.glade.get_widget('aboutdialog')
      dialog.run()
      dialog.hide()

   def on_debug_toggled(self,wdg):
      if (wdg.get_active()):
         self.pyk_debug=1
         self.py2yo("pyk_set pyk_debug 1")
      else:
         self.pyk_debug=0
         self.py2yo("pyk_set pyk_debug 0")

   def on_show_wfss_and_dms_toggled(self,wdg):
      show_state = self.glade.get_widget('show_wfss_and_dms').get_active()
      if (show_state):
         try:
            s = self.size
         except:
            s = 0
         self.size = self.window.get_size()
         self.glade.get_widget('wfs_and_dms').show()
         if (s):
            self.window.resize(s[0],s[1])
      else:
         s = self.size
         self.size = self.window.get_size()
         self.glade.get_widget('wfs_and_dms').hide()
         self.window.resize(s[0],s[1])
         
   def  on_create_phase_screens_activate(self,wdg):
      self.set_cursor_busy(1)
      self.py2yo("wrap_create_phase_screens")
   
   #
   # EDITORS
   #
   
   def on_modified_editor(self,wdg):
      self.editor.set_title(' * '+self.yaoparfile)   

   def on_editor_save_activate(self,wdg):
      # get buffer content
      textarea = self.glade.get_widget('textarea')
      self.buffer = textarea.get_buffer()
      params=self.buffer.get_text(self.buffer.get_start_iter(),self.buffer.get_end_iter())
      # contruct/set new yaoparfile
      #self.yaoparfile = self.get_yaoparfile_stamped()
      # set main window yaoparfile entry to the new name
      #self.glade.get_widget('yaoparfile').set_text(self.yaoparfile)
      # save file
      f = open(self.yaopardir+'/'+self.yaoparfile,'w')
      f.write(params)
      f.close()
      # refresh editor window title
      self.editor.set_title(self.yaoparfile)
      self.on_aoread_clicked(wdg)

   def on_editor_save_as_activate(self,wdg):
      chooser = gtk.FileChooserDialog(title='Save YAO parfile as...',action=gtk.FILE_CHOOSER_ACTION_SAVE,buttons=(gtk.STOCK_CANCEL,gtk.RESPONSE_CANCEL,gtk.STOCK_SAVE,gtk.RESPONSE_OK))
      res = chooser.run()
      if res == gtk.RESPONSE_CANCEL:
         chooser.destroy()
         return
      name = chooser.get_filename()
      self.yaopardir = os.path.dirname(name)
      self.yaoparfile = os.path.basename(name)
      chooser.destroy()
      # get buffer content
      textarea = self.glade.get_widget('textarea')
      self.buffer = textarea.get_buffer()
      params=self.buffer.get_text(self.buffer.get_start_iter(),self.buffer.get_end_iter())
      # save the content in tmp.par
      f = open(self.yaopardir+'/'+self.yaoparfile,'w')
      f.write(params)
      f.close()
      # get new name
      # self.yaoparfile = self.get_yaoparfile_stamped()
      # set parfile name to yaoparfile in main window
      self.glade.get_widget('yaoparfile').set_text(self.yaoparfile)
      # apply (triger aoread())
      self.on_aoread_clicked(wdg)

   def on_editor_close_activate(self,wdg,*args):
      self.editor.hide()
      return True
   
   def get_yaoparfile_stamped(self):
      tmp = self.yaoparfile.replace('.par','')
      tmp = tmp+'_'+strftime("%Y%B%d_%H%M%S",localtime())+'.par'
      return tmp

   def on_modified_editor2(self,wdg):
      self.editor2.set_title(' * yao.conf')

   def on_editor2_save_activate(self,wdg):
      # get buffer content
      textarea2 = self.glade.get_widget('textarea2')
      self.buffer2 = textarea2.get_buffer()
      params=self.buffer2.get_text(self.buffer2.get_start_iter(),self.buffer2.get_end_iter())
      # save file
      f = open(self.yuserdir+'yao.conf','w')
      f.write(params)
      f.close()
      self.py2yo('read_conf')
      self.editor2.set_title('yao.conf')

   def on_editor2_close_activate(self,wdg,*args):
      self.editor2.hide()
      return True    
   
   #
   # Main Panel Events Handlers
   #
      
   def on_yaoparfile_activate(self,wdg):
      self.glade.get_widget('edit').set_sensitive(1)
      self.glade.get_widget('aoread').set_sensitive(1)
      self.glade.get_widget('aoinit').set_sensitive(0)
      self.glade.get_widget('aoloop').set_sensitive(0)
      self.glade.get_widget('go').set_sensitive(0)
      self.glade.get_widget('pause').set_sensitive(0)
      self.glade.get_widget('step').set_sensitive(0)
      self.glade.get_widget('restart').set_sensitive(0)
      self.glade.get_widget('displays').set_sensitive(0)
      self.glade.get_widget('wfss').set_sensitive(0)
      self.glade.get_widget('dms').set_sensitive(0)
      self.glade.get_widget('seeingframe').set_sensitive(0)
      self.window.set_title(wdg.get_text())
      self.glade.get_widget('aoread').grab_focus()

   def on_yaoparfile_select_clicked(self,wdg):
      chooser = gtk.FileChooserDialog(title='YAO parfile selection',action=gtk.FILE_CHOOSER_ACTION_OPEN,buttons=(gtk.STOCK_CANCEL,gtk.RESPONSE_CANCEL,gtk.STOCK_OPEN,gtk.RESPONSE_OK))
      filter = gtk.FileFilter()
      filter.add_pattern('*.par')
      filter.set_name('YAO parfiles')
      chooser.add_filter(filter)
      chooser.set_current_folder(self.yaopardir)
      res = chooser.run()
      if res == gtk.RESPONSE_OK:
         self.yaoparfile=chooser.get_filename()
         self.yaoparfile=self.yaoparfile.split('/')[-1]
         self.glade.get_widget('yaoparfile').set_text(self.yaoparfile)
         self.yaopardir = chooser.get_current_folder()
         self.window.set_title(self.yaoparfile)
      chooser.destroy()
      self.glade.get_widget('edit').set_sensitive(1)
      self.glade.get_widget('aoread').set_sensitive(1)
      self.glade.get_widget('aoinit').set_sensitive(0)
      self.glade.get_widget('aoloop').set_sensitive(0)
      self.glade.get_widget('go').set_sensitive(0)
      self.glade.get_widget('pause').set_sensitive(0)
      self.glade.get_widget('step').set_sensitive(0)
      self.glade.get_widget('restart').set_sensitive(0)
      self.glade.get_widget('displays').set_sensitive(0)
      self.glade.get_widget('wfss').set_sensitive(0)
      self.glade.get_widget('dms').set_sensitive(0)
      self.glade.get_widget('seeingframe').set_sensitive(0)
      self.glade.get_widget('aoread').grab_focus()

   def on_edit_clicked(self,wdg):
      self.editor.set_title(self.yaoparfile)
      self.editor.show()
      f = open(self.yaopardir+'/'+self.yaoparfile,'r')
      params = f.read()
      f.close()
      self.buffer.set_text(params)
      self.buffer.connect('changed',self.on_modified_editor)
      textarea = self.glade.get_widget('textarea')
      textarea.set_buffer(self.buffer)

   def on_edit2_activate(self,wdg):
      self.editor2.set_title('yao.conf')
      self.editor2.show()
      f = open(self.yuserdir+'yao.conf','r')
      params = f.read()
      f.close()
      self.buffer2.set_text(params)
      self.buffer2.connect('changed',self.on_modified_editor2)
      textarea2 = self.glade.get_widget('textarea2')
      textarea2.set_buffer(self.buffer2)

   def on_aoread_clicked(self,wdg):
      self.set_cursor_busy(1)
      yaoparfile=self.glade.get_widget('yaoparfile').get_text()
      self.py2yo('pyk_set yaopardir "%s"' % self.yaopardir)
      self.py2yo('pyk_set yaoparfile "%s"' % yaoparfile)
      self.py2yo('wrap_aoread')
      self.glade.get_widget('aoinit').set_sensitive(1)
      self.glade.get_widget('aoloop').set_sensitive(0)
      self.glade.get_widget('go').set_sensitive(0)
      self.glade.get_widget('pause').set_sensitive(0)
      self.glade.get_widget('step').set_sensitive(0)
      self.glade.get_widget('restart').set_sensitive(0)
      self.glade.get_widget('displays').set_sensitive(0)
      self.glade.get_widget('seeingframe').set_sensitive(0)
      self.glade.get_widget('aoinit').grab_focus()
      self.glade.get_widget('wfss').set_sensitive(1)
      self.glade.get_widget('dms').set_sensitive(1)
      self.progressbar.set_fraction(0.)
      self.progressbar.set_text('')
      self.init = 0
      self.py2yo('pyk_set initdone 0')

   def on_aoinit_clicked(self,wdg):
      self.set_cursor_busy(1)
      self.set_aoinit_flags()
      self.py2yo('do_aoinit_disp')
      self.glade.get_widget('aoloop').set_sensitive(1)
      self.glade.get_widget('aoloop').grab_focus()
      self.progressbar.set_fraction(0.)
      self.progressbar.set_text('')
      self.init = 1
      self.py2yo('pyk_set initdone 1')
   
   def on_aoinit_popup_button_clicked(self,wdg,event):
      if event.type == gtk.gdk.BUTTON_PRESS:
         wdg.popup(None, None, None, event.button, event.time)
         # Tell calling code that we have handled this event the buck
         return True
      # Tell calling code that we have not handled this event pass it on.
      return False

   def set_aoinit_flags(self):
      disp = self.glade.get_widget('aoinit_disp').get_active()
      clean = self.glade.get_widget('aoinit_clean').get_active()
      forcemat = self.glade.get_widget('aoinit_forcemat').get_active()
      svd = self.glade.get_widget('aoinit_svd').get_active()
      keepdmconfig = self.glade.get_widget('aoinit_keepdmconfig').get_active()
      self.py2yo('set_aoinit_flags %d %d %d %d %d' % (disp,clean,forcemat,svd,keepdmconfig))
   
   def set_aoloop_flags(self):
      disp = self.glade.get_widget('aoloop_disp').get_active()
      savecb = self.glade.get_widget('aoloop_savecb').get_active()
      self.py2yo('set_aoloop_flags %d %d' % (disp,savecb))
   
   def on_aoloop_popup_button_clicked(self,wdg,event):
      if event.type == gtk.gdk.BUTTON_PRESS:
         wdg.popup(None, None, None, event.button, event.time)
         # Tell calling code that we have handled this event the buck
         return True
      # Tell calling code that we have not handled this event pass it on.
      return False
   
   def on_aoloop_clicked(self,wdg):
      self.set_cursor_busy(1)
      self.glade.get_widget('go').set_sensitive(1)
      self.glade.get_widget('pause').set_sensitive(1)
      self.glade.get_widget('step').set_sensitive(1)
      self.glade.get_widget('restart').set_sensitive(1)
      self.glade.get_widget('displays').set_sensitive(1)
      self.glade.get_widget('wfss').set_sensitive(1)
      self.glade.get_widget('dms').set_sensitive(1)
      self.glade.get_widget('seeingframe').set_sensitive(1)
#      self.glade.get_widget('loopgain').set_sensitive(1)
#      self.glade.get_widget('imlambda').set_sensitive(1)
      self.glade.get_widget('go').grab_focus()
      disprate = self.glade.get_widget('disp_rate').get_value()
      self.py2yo('do_aoloop_disp %d' % disprate)

   def on_go_clicked(self,wdg):
      #self.py2yo('toggle_animate 1')
      self.py2yo('cont')

   def on_pause_clicked(self,wdg):
      self.py2yo('stop')
      #self.py2yo('toggle_animate 0')

   def on_step_clicked(self,wdg):
      self.py2yo('go 1')
      
   def on_restart_clicked(self,wdg):
      self.glade.get_widget('go').set_sensitive(1)
      self.glade.get_widget('pause').set_sensitive(1)
      self.glade.get_widget('step').set_sensitive(1)
      self.py2yo('do_aoloop_disp')

   def aoloop_to_end(self):
      self.glade.get_widget('go').set_sensitive(0)
      self.glade.get_widget('pause').set_sensitive(0)
      self.glade.get_widget('step').set_sensitive(0)

   #
   # Display Panel Events Handlers
   #
      
   def disp_panel_set_sensitivity(self,sens):
      self.glade.get_widget('disp_rate').set_sensitive(sens)
      self.glade.get_widget('disp_label').set_sensitive(sens)
      self.glade.get_widget('instavg_label').set_sensitive(sens)
      self.glade.get_widget('instavg_hbox').set_sensitive(sens)
      
   def on_disp_pause_clicked(self,wdg):
      #self.py2yo('toggle_animate 0')
      self.py2yo('fma')
      self.py2yo('funcset dispFlag 10000')
      self.disp_panel_set_sensitivity(0)
      self.glade.get_widget('disp_pause').set_sensitive(0)
      self.glade.get_widget('disp_resume').set_sensitive(1)

   def on_disp_resume_clicked(self,wdg):
      self.py2yo('funcset dispFlag %d' % self.dispflag)
      #self.py2yo('toggle_animate 1')
      self.disp_panel_set_sensitivity(1)
      self.glade.get_widget('disp_pause').set_sensitive(1)
      self.glade.get_widget('disp_resume').set_sensitive(0)

   def on_image_disp_inst_clicked(self,wdg):
      if wdg.get_active():
         self.py2yo('toggle_im_imav 0')

   def on_image_disp_avg_clicked(self,wdg):
      if wdg.get_active():
         self.py2yo('toggle_im_imav 1')

   def on_disp_rate_value_changed(self,wdg):
      self.dispflag = wdg.get_value()
      self.py2yo('funcset dispFlag %d' % self.dispflag)

   #
   # General parameters Panel Events Handlers
   #
      
   def on_seeing_value_changed(self,wdg):
      self.py2yo('change_seeing %f' % wdg.get_value())
      
   def on_loopgain_value_changed(self,wdg):
      if self.init:
         self.py2yo('set_loop_gain %f' % wdg.get_value())

   def on_imlambda_value_changed(self,wdg):
      if self.init:
         self.py2yo('change_target_lambda %f' % wdg.get_value())

   #
   # WFS parameters Panel Events Handlers
   #

   def on_subtract_background_toggled(self,wdg):
      if self.init:
         self.py2yo('wfs_subtract_background %d' % wdg.get_active())

   def on_noise_toggled(self,wdg):
      if self.init:
#      if wdg.get_active():
#         self.py2yo('set_wfs_noise 1')
#      else:
#         self.py2yo('set_wfs_noise 0')
         self.py2yo('set_wfs_noise %d' % wdg.get_active())
      
   def on_correct_up_tt_toggled(self,wdg):
      if self.init:
         self.py2yo('wfs_set_uptt %d' % wdg.get_active())
   
   def on_efd_value_changed(self,wdg):
      if self.init:
         self.py2yo('set_wfs_efd %f' % wdg.get_value())

   def on_pyr_mod_value_changed(self,wdg):
      if self.init:
         self.py2yo('set_wfs_pyr_mod %f' % wdg.get_value())

   def on_gsmag_value_changed(self,wdg):
      if self.init:
         self.py2yo('set_gs_mag %f' % wdg.get_value())
         
   def on_gsalt_value_changed(self,wdg):
      if self.init:
         self.py2yo('set_gs_alt %f' % wdg.get_value())
         
   def on_gsdepth_value_changed(self,wdg):
      if self.init:
         self.py2yo('set_gs_depth %f' % wdg.get_value())
         
   def on_ron_value_changed(self,wdg):
      if self.init:
         self.py2yo('set_wfs_ron %f' % wdg.get_value())

   def on_sh_threshold_value_changed(self,wdg):
      if self.init:
         self.py2yo('set_wfs_threshold %f' % wdg.get_value())
      
   def on_sh_kernel_value_changed(self,wdg):
      if self.init:
         self.py2yo('set_wfs_kernel %f' % wdg.get_value())

   def on_ninteg_cycles_value_changed(self,wdg):
      if self.init:
         self.py2yo('set_wfs_nintegcycles %f' % wdg.get_value())
   
   def y_set_nwfs(self,nwfs):
      self.nwfs = nwfs
      for i in range(nwfs):
         self.glade.get_widget('wfs'+str(i+1)).show()
         self.glade.get_widget('wfs'+str(i+1)).set_active(0)
         
   def on_wfs_select_toggled(self,wdg):
      for i in range(self.nwfs):
         ok = self.glade.get_widget('wfs'+str(i+1)).get_active()
         self.py2yo('set_okwfs %d %d' % (i+1,ok))

   def wfs_panel_set_sensitivity(self,sens,wfstype):
      # sens = 0 or 1
      # wfstype = 1 (sh), 2 (curvature) or 3 (pyramid) (0=undefined)
      self.glade.get_widget('subtract_background').set_sensitive(sens)
      self.glade.get_widget('noise').set_sensitive(sens)
      self.glade.get_widget('correct_up_tt').set_sensitive(sens)
      self.glade.get_widget('efd').set_sensitive(sens)
      self.glade.get_widget('pyr_mod').set_sensitive(sens)
      self.glade.get_widget('gsmag').set_sensitive(sens)
      self.glade.get_widget('gsalt').set_sensitive(sens)
      self.glade.get_widget('gsdepth').set_sensitive(sens)
      self.glade.get_widget('ron').set_sensitive(sens)
      self.glade.get_widget('sh_threshold').set_sensitive(sens)
      self.glade.get_widget('sh_kernel').set_sensitive(sens)
      self.glade.get_widget('ninteg_cycles').set_sensitive(sens)
      if (wfstype!=1):
         self.glade.get_widget('sh_threshold').set_sensitive(0)
         self.glade.get_widget('sh_kernel').set_sensitive(0)
         self.glade.get_widget('gsalt').set_sensitive(0)
         self.glade.get_widget('gsdepth').set_sensitive(0)
         self.glade.get_widget('correct_up_tt').set_sensitive(0)
      if (wfstype!=2):
         self.glade.get_widget('efd').set_sensitive(0)
      if (wfstype!=3):
         self.glade.get_widget('pyr_mod').set_sensitive(0)
         
         
   #
   # DM parameters Panel Events Handlers
   #
         
   def on_dmreset_clicked(self,wdg):
      self.py2yo('dm_reset')

   def on_dmflatten_clicked(self,wdg):
      self.py2yo('dm_flatten')

   def on_dmgain_value_changed(self,wdg):
      self.py2yo('dm_gain %f' % wdg.get_value())
         
   def on_xmisreg_value_changed(self,wdg):
      self.py2yo('dm_xmisreg %f' % wdg.get_value())
         
   def on_ymisreg_value_changed(self,wdg):
      self.py2yo('dm_ymisreg %f' % wdg.get_value())

   def on_sat_voltage_value_changed(self,wdg):
      self.py2yo('dm_satvolt %f' % wdg.get_value())


   def y_set_ndm(self,ndm):
      self.ndm = ndm
      for i in range(ndm):
         self.glade.get_widget('dm'+str(i+1)).show()
         self.glade.get_widget('dm'+str(i+1)).set_active(0)

   def on_dm_select_toggled(self,wdg):
      for i in range(self.ndm):
         ok = self.glade.get_widget('dm'+str(i+1)).get_active()
         self.py2yo('set_okdm %d %d' % (i+1,ok))
         
   def dm_panel_set_sensitivity(self,sens):
      # sens = 0 or 1
      self.glade.get_widget('dmreset').set_sensitive(sens)
      self.glade.get_widget('dmflatten').set_sensitive(sens)
      self.glade.get_widget('extrapolated').set_sensitive(sens)
      self.glade.get_widget('dmgain').set_sensitive(sens)
      self.glade.get_widget('xmisreg').set_sensitive(sens)
      self.glade.get_widget('ymisreg').set_sensitive(sens)
      self.glade.get_widget('sat_voltage').set_sensitive(sens)
                  
   #
   # Yorick to Python Wrapper Functions
   #

   def y_parm_update(self,name,val):
      self.glade.get_widget(name).set_value(val)

   def y_text_parm_update(self,name,txt):
      self.glade.get_widget(name).set_text(txt)

   def y_set_checkbutton(self,name,val):
      self.glade.get_widget(name).set_active(val)
      
   def pyk_error(self,msg):
      dialog = gtk.MessageDialog(type=gtk.MESSAGE_ERROR,buttons=gtk.BUTTONS_OK,message_format=msg)
      dialog.run()
      dialog.destroy()
      
   def pyk_info(self,msg):
      dialog = gtk.MessageDialog(type=gtk.MESSAGE_INFO,buttons=gtk.BUTTONS_OK,message_format=msg)
      dialog.run()
      dialog.destroy()

   def pyk_info_w_markup(self,msg):
      dialog = gtk.MessageDialog(type=gtk.MESSAGE_INFO,buttons=gtk.BUTTONS_OK)
      dialog.set_markup(msg)
#      dialog.set_size_request(600,-1)
      dialog.run()
      dialog.destroy()

   def pyk_warning(self,msg):
      dialog = gtk.MessageDialog(type=gtk.MESSAGE_WARNING,buttons=gtk.BUTTONS_OK,message_format=msg)
      dialog.run()
      dialog.destroy()
      
   def on_quit1_activate(self,*args):
      # tell the ascam image server to quit (or not?)
      self.py2yo('yaopy_quit')
#      raise SystemExit
   
   def on_window1_map_event(self,wdg,*args):
      drawingarea = self.glade.get_widget('drawingarea1')
      mwid = drawingarea.window.xid;
      self.py2yo('yao_win_init %d' % mwid)
      # update parameters from yorick:
      self.py2yo('gui_update')
   
#   def on_drawingarea1_map_event(self,wdg,*args):   
#      # only reparent once the widget is mapped.
#      # note this only need to happen once
#      mwid = wdg.window.xid;
#      self.py2yo('yao_win_init %d' % mwid)
#      # update parameters from yorick:
#      self.py2yo('gui_update')
#
      return False

   #
   # minimal wrapper for yorick/python communication
   #

   def yo2py_flush(self):
      sys.stdin.flush()
   
   def py2yo(self,msg):
      # sends string command to yorick's eval
      sys.stdout.write(msg+'\n')
      sys.stdout.flush()
      
   def yo2py(self,cb_condition,*args):
      if cb_condition == gobject.IO_HUP:
         raise SystemExit, "lost pipe to yorick"
      # handles string command from yorick
      # note: individual message needs to end with /n for proper ungarbling
      while 1:
         try:
            msg = sys.stdin.readline()
            msg = "self."+msg
            #  self.py2yo('\"%s\"' % msg)
            try:
               exec(msg)
            except Exception, e:
               sys.stderr.write('yo2py eval: '+str(e)+'\n')
         except IOError, e:
            if e.errno == errno.EAGAIN:
               # the pipe's empty, good
               break
            # else bomb out
            raise SystemExit, "yo2py unexpected IOError:" + str(e)
         except Exception, ee:
            raise SystemExit, "yo2py unexpected Exception:" + str(ee)
      # carefull with the ident here
      return True


   def set_cursor_busy(self,state):
      if state:
         self.window.window.set_cursor(gtk.gdk.Cursor(gtk.gdk.WATCH))
      else:
         self.window.window.set_cursor(gtk.gdk.Cursor(gtk.gdk.LEFT_PTR))
         

if len(sys.argv) != 3:
   print 'Usage: yao.py path_to_glade dpi'
   raise SystemExit

path2glade = str(sys.argv[1])
dpi = int(sys.argv[2])
top = yao(path2glade,dpi)
