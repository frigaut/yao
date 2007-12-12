#!/usr/bin/env python
# yao.py
# To Do:
# - implement gtk dialog to warn user of danger of re-using the same file name
#   when editing a current par file and saving

import gtk
import gtk.glade
import sys
import gobject
import os, fcntl, errno
from time import *

class yao:
   
   def destroy(self, wdg, data=None):
      self.py2yo('quit')
      gtk.main_quit()
      
   def __init__(self,yaotop):
      self.yaotop = yaotop
      self.usercmd = 'STOP'
      
      # callbacks and glade UI
      dic = {
         'on_about_activate': self.on_about_activate,
         'on_quit1_activate' : self.on_quit1_activate,
         'on_show_wfss_and_dms_toggled' : self.on_show_wfss_and_dms_toggled,
         'on_drawingarea1_map_event' : self.on_drawingarea1_map_event,
         'on_edit_clicked' : self.on_edit_clicked,
         'on_create_phase_screens_activate': self.on_create_phase_screens_activate,
         # MAIN
         'on_aoread_clicked' : self.on_aoread_clicked,
         'on_yaoparfile_activate': self.on_yaoparfile_activate,
         'on_yaoparfile_changed': self.on_yaoparfile_activate,
         'on_yaoparfile_select_clicked' : self.on_yaoparfile_select_clicked,
         'on_aoinit_clicked' : self.on_aoinit_clicked,
         'on_aoloop_clicked' : self.on_aoloop_clicked,
         'on_go_clicked' : self.on_go_clicked,
         'on_pause_clicked' : self.on_pause_clicked,
         'on_step_clicked' : self.on_step_clicked,
         'on_restart_clicked' : self.on_restart_clicked,
         # DISPLAYS
         'on_disp_pause_clicked' : self.on_disp_pause_clicked,
         'on_disp_resume_clicked' : self.on_disp_resume_clicked,
         'on_image_disp_inst_clicked' : self.on_image_disp_inst_clicked,
         'on_image_disp_avg_clicked' : self.on_image_disp_avg_clicked,
         'on_disp_rate_value_changed' : self.on_disp_rate_value_changed,
         # GEN PAR
         'on_loopgain_value_changed' : self.on_loopgain_value_changed,
         'on_imlambda_value_changed' : self.on_imlambda_value_changed,
         'on_seeing_value_changed' : self.on_seeing_value_changed,
         # DM
         'on_dmreset_clicked': self.on_dmreset_clicked,
         'on_dmflatten_clicked': self.on_dmflatten_clicked,
         'on_dmgain_value_changed' : self.on_dmgain_value_changed,
         'on_xmisreg_value_changed' : self.on_xmisreg_value_changed,
         'on_ymisreg_value_changed' : self.on_ymisreg_value_changed,
         'on_sat_voltage_value_changed' : self.on_sat_voltage_value_changed,
         'on_dm_select_toggled': self.on_dm_select_toggled,
         # WFS
         'on_subtract_background_toggled' : self.on_subtract_background_toggled,
         'on_noise_toggled': self.on_noise_toggled,
         'on_correct_up_tt_toggled' : self.on_correct_up_tt_toggled,
         'on_efd_value_changed' : self.on_efd_value_changed,         
         'on_gsmag_value_changed' : self.on_gsmag_value_changed,
         'on_gsalt_value_changed' : self.on_gsalt_value_changed,
         'on_gsdepth_value_changed' : self.on_gsdepth_value_changed,
         'on_ron_value_changed': self.on_ron_value_changed,
         'on_sh_threshold_value_changed' : self.on_sh_threshold_value_changed,
         'on_sh_kernel_value_changed' : self.on_sh_kernel_value_changed,
         'on_ninteg_cycles_value_changed' : self.on_ninteg_cycles_value_changed,
         'on_wfs_select_toggled': self.on_wfs_select_toggled,
         # EDITOR
         'on_editor_save_activate' : self.on_editor_save_activate,
         'on_editor_save_as_activate' : self.on_editor_save_as_activate,
         'on_editor_close_activate' : self.on_editor_close_activate,
#         'on_save_warning_dialog_response' : self.on_save_warning_dialog_response,
         }
      
#      self.yaotop = os.environ['YAOTOP']
      self.glade = gtk.glade.XML(os.path.join(self.yaotop,'glade/yao.glade')) 
      self.window = self.glade.get_widget('window1')
      if (self.window):
         self.window.connect('destroy', self.destroy)
      self.glade.signal_autoconnect(dic)

      self.editor = self.glade.get_widget('window2')
      if (self.editor):
         self.editor.connect('delete_event', self.on_editor_close_activate)

      self.statusbar = self.glade.get_widget('statusbar')
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
      gobject.io_add_watch(sys.stdin,gobject.IO_IN|gobject.IO_HUP,self.yo2py,None)

      self.wfs_panel_set_sensitivity(0,0)
      #self.glade.get_widget('wfs_and_dms').hide()
      
      # run
      gtk.main()


   dispflag = 1
   init = 0
   buffer = gtk.TextBuffer()
   
   def on_about_activate(self,wdg):
      dialog = self.glade.get_widget('aboutdialog')
      dialog.run()
      dialog.hide()

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
   # EDITOR
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
      f = open(self.yaoparfile,'w')
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
      self.yaoparfile = chooser.get_filename()
      chooser.destroy()
      # get buffer content
      textarea = self.glade.get_widget('textarea')
      self.buffer = textarea.get_buffer()
      params=self.buffer.get_text(self.buffer.get_start_iter(),self.buffer.get_end_iter())
      # save the content in tmp.par
      f = open(self.yaoparfile,'w')
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
      res = chooser.run()
      if res == gtk.RESPONSE_OK:
         self.yaoparfile=chooser.get_filename()
         self.yaoparfile=self.yaoparfile.split('/')[-1]
         self.glade.get_widget('yaoparfile').set_text(self.yaoparfile)
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
      self.window.set_title(self.yaoparfile)
      self.glade.get_widget('aoread').grab_focus()

   def on_edit_clicked(self,wdg):
      self.editor.set_title(self.yaoparfile)
      self.editor.show()
      f = open(self.yaoparfile,'r')
      params = f.read()
      f.close()
      self.buffer.set_text(params)
      self.buffer.connect('changed',self.on_modified_editor)
      textarea = self.glade.get_widget('textarea')
      textarea.set_buffer(self.buffer)

   def on_aoread_clicked(self,wdg):
      self.set_cursor_busy(1)
      yaoparfile=self.glade.get_widget('yaoparfile').get_text()
      self.py2yo('pyk_set yaoparfile "%s"' % yaoparfile)
      self.py2yo('wrap_aoread \"'+yaoparfile+'\"')
      self.glade.get_widget('aoinit').set_sensitive(1)
      self.glade.get_widget('aoloop').set_sensitive(0)
      self.glade.get_widget('go').set_sensitive(0)
      self.glade.get_widget('pause').set_sensitive(0)
      self.glade.get_widget('step').set_sensitive(0)
      self.glade.get_widget('restart').set_sensitive(0)
      self.glade.get_widget('displays').set_sensitive(0)
      self.glade.get_widget('wfss').set_sensitive(0)
      self.glade.get_widget('dms').set_sensitive(0)
      self.glade.get_widget('seeingframe').set_sensitive(0)
      self.glade.get_widget('aoinit').grab_focus()
      self.progressbar.set_fraction(0.)
      self.progressbar.set_text('')
      self.init = 0

   def on_aoinit_clicked(self,wdg):
      self.set_cursor_busy(1)
      self.set_aoinit_flags()
      self.py2yo('do_aoinit_disp')
      self.glade.get_widget('aoloop').set_sensitive(1)
      self.glade.get_widget('aoloop').grab_focus()
      self.progressbar.set_fraction(0.)
      self.progressbar.set_text('')
      self.init = 1
   
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
      self.py2yo('do_aoloop_disp')

   def on_go_clicked(self,wdg):
      self.py2yo('toggle_animate 1')
      self.py2yo('go')

   def on_pause_clicked(self,wdg):
      self.py2yo('stop')
      self.py2yo('toggle_animate 0')

   def on_step_clicked(self,wdg):
      self.py2yo('go 1')
      
   def on_restart_clicked(self,wdg):
      self.py2yo('do_aoloop_disp')
      
   #
   # Display Panel Events Handlers
   #
      
   def on_disp_pause_clicked(self,wdg):
      self.py2yo('toggle_animate 0')
      self.py2yo('fma')
      self.py2yo('funcset dispFlag 10000')
      
   def on_disp_resume_clicked(self,wdg):
      self.py2yo('funcset dispFlag %d' % self.dispflag)
      self.py2yo('toggle_animate 1')

   def on_image_disp_inst_clicked(self,wdg):
      if wdg.get_active():
         self.py2yo('toggle_im_imav 1')

   def on_image_disp_avg_clicked(self,wdg):
      if wdg.get_active():
         self.py2yo('toggle_im_imav 2')

   def on_disp_rate_value_changed(self,wdg):
      self.dispflag = wdg.get_value()
      self.py2yo('funcset dispFlag %d' % self.dispflag)

   #
   # General parameters Panel Events Handlers
   #
      
   def on_seeing_value_changed(self,wdg):
      if self.init:
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
      # wfstype = 1 (sh) or 2 (curvature) (0=undefined)
      self.glade.get_widget('subtract_background').set_sensitive(sens)
      self.glade.get_widget('noise').set_sensitive(sens)
      self.glade.get_widget('correct_up_tt').set_sensitive(sens)
      self.glade.get_widget('efd').set_sensitive(sens)
      self.glade.get_widget('gsmag').set_sensitive(sens)
      self.glade.get_widget('gsalt').set_sensitive(sens)
      self.glade.get_widget('gsdepth').set_sensitive(sens)
      self.glade.get_widget('ron').set_sensitive(sens)
      self.glade.get_widget('sh_threshold').set_sensitive(sens)
      self.glade.get_widget('sh_kernel').set_sensitive(sens)
      self.glade.get_widget('ninteg_cycles').set_sensitive(sens)
      if (wfstype==1):
         self.glade.get_widget('efd').set_sensitive(0)
      if (wfstype==2):
         self.glade.get_widget('sh_threshold').set_sensitive(0)
         self.glade.get_widget('sh_kernel').set_sensitive(0)
         self.glade.get_widget('gsalt').set_sensitive(0)
         self.glade.get_widget('gsdepth').set_sensitive(0)
         self.glade.get_widget('correct_up_tt').set_sensitive(0)
         
         
   #
   # DM parameters Panel Events Handlers
   #
         
   def on_dmreset_clicked(self,wdg):
      self.py2yo('dm_reset')

   def on_dmflatten_clicked(self,wdg):
      self.py2yo('dm_flatten')

   def on_dmgain_value_changed(self,wdg):
      if self.init:
         self.py2yo('dm_gain %f' % wdg.get_value())
         
   def on_xmisreg_value_changed(self,wdg):
      if self.init:
         self.py2yo('dm_xmisreg %f' % wdg.get_value())
         
   def on_ymisreg_value_changed(self,wdg):
      if self.init:
         self.py2yo('dm_ymisreg %f' % wdg.get_value())

   def on_sat_voltage_value_changed(self,wdg):
      if self.init:
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
      self.py2yo('quit')
#      raise SystemExit
   
   def on_drawingarea1_map_event(self,wdg,*args):   
      # only reparent once the widget is mapped.
      # note this only need to happen once
      mwid = wdg.window.xid;
      self.py2yo('yao_win_init %d' % mwid)
      # update parameters from yorick:
      self.py2yo('gui_update')

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
      # note: inidividual message needs to end with /n for proper ungarbling
      while 1:
         try:
            msg = sys.stdin.readline()
            msg = "self."+msg
            exec(msg)
         except IOError, e:
            if e.errno == errno.EAGAIN:
               # the pipe's empty, good
               break
            # else bomb out
            raise SystemExit, "yo2py unexpected IOError:" + str(e)
         except Exception, ee:
            raise SystemExit, "yo2py unexpected Exception:" + str(ee)
         return True

   def set_cursor_busy(self,state):
      if state:
         self.window.window.set_cursor(gtk.gdk.Cursor(gtk.gdk.WATCH))
      else:
         self.window.window.set_cursor(gtk.gdk.Cursor(gtk.gdk.LEFT_PTR))
         

if len(sys.argv) != 2:
   print 'Usage: yao.py path_to_yao'
   raise SystemExit

yaotop = str(sys.argv[1])
top = yao(yaotop)
