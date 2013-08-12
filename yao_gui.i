/*
 * yao_gui.i
 *
 * Routines for a quick and dirty GUI to report on progress for yao
 *
 * This file is part of the yao package, an adaptive optics simulation tool.
 *
 * Copyright (c) 2002-2013, Francois Rigaut
 *
 * This program is free software; you can redistribute it and/or  modify it
 * under the terms of the GNU General Public License  as  published  by the
 * Free Software Foundation; either version 2 of the License,  or  (at your
 * option) any later version.
 *
 * This program is distributed in the hope  that  it  will  be  useful, but
 * WITHOUT  ANY   WARRANTY;   without   even   the   implied   warranty  of
 * MERCHANTABILITY or  FITNESS  FOR  A  PARTICULAR  PURPOSE.   See  the GNU
 * General Public License for more details (to receive a  copy  of  the GNU
 * General Public License, write to the Free Software Foundation, Inc., 675
 * Mass Ave, Cambridge, MA 02139, USA).
 *
 */

func progress_bar(percent,&id,init=)
/* DOCUMENT func progress_bar(percent,&id,init=)
   percent: progress -in percent- to be displayed
   id: internal use. Returned by init. To be used in further calls.
   init = [x0,y0,x1,y1]: position to display the progress bar
   NOTE: One has to use the palette "earth.gp" to get a nice blue progress bar
   example:
   
   SEE ALSO:
 */
{
  extern PBids,PBcoords;

  if (!is_void(init)) {
    if (numberof(init) != 4) {
      error,"progress_bar init keyword takes exactly 4 parameters (x0,y0,x1,y1)";
    }
    if (is_void(PBids)) {
      PBids=1;
      PBcoords = init(,-);
    } else {
      grow,PBids,PBids(0)+1;
      grow,PBcoords,init;
    }
    id = PBids(0);
  }
  coo = PBcoords(,id);
  if (!is_void(init)) {
    plg,[coo(2),coo(2),coo(4),coo(4),coo(2)],[coo(1),coo(3),coo(3),coo(1),coo(1)],marks=0;
    return;
  }
  pli,array('\27',[2,2,2]),coo(1),coo(2),coo(1)+(coo(3)-coo(1))*percent/100.,coo(4);
  plg,[coo(2),coo(2),coo(4),coo(4),coo(2)],[coo(1),coo(3),coo(3),coo(1),coo(1)],
    width=2,marks=0;
  limits,0.,1.,0.,1.;
}

func draw_button(color,pos,size)
{
  n = 32;
  if (color=="green")  but = char(255-(dist(n) <(n/2))*5);
  if (color=="yellow") but = char(255-(dist(n) <(n/2))*9);
  if (color=="red")    but = char(255-(dist(n) <(n/2))*4);
  pli,but,pos(1)-size/2.,pos(2)-size/2.,pos(1)+size/2.,pos(2)+size/2.;
  limits,0.,1.,0.,1.;
}

func draw_status(status,&id,init=)
/* DOCUMENT func draw_status(status,&id,init=)
   status = 1 (green), 0 (yellow), -1 (red)
   id = returned after init, use as is in futher calls.
   init = [xpos,ypos,size]
   example:
   first call/init:
   draw_status,1,ids2,init=[0.3,0.6,0.02]
   init the light and place it at positions (0.3,0.6) with size 0.02
   draw red/yellow/green light:
   draw_status,1,ids2
   SEE ALSO:
 */
{
  extern Sids,Scoords;

  if (!is_void(init)) {
    if (numberof(init) != 3) {
      error,"progress_bar init keyword takes exactly 3 parameters (xpos,ypos,size)";
    }
    if (is_void(Sids)) {
      Sids=1;
      Scoords = init(,-);
    } else {
      grow,Sids,Sids(0)+1;
      grow,Scoords,init;
    }
    id = Sids(0);
  }
  coo = Scoords(,id);
  s = span(0.,2*pi,30);
  //  if (!is_void(init)) {
  //    plg,coo(3)*cos(s)/2+coo(2),coo(3)*sin(s)/2+coo(1),width=3,marks=0;
  //    return;
  //  }
  color = ["green","yellow","red"](-status+2);
  draw_button,color,[coo(1),coo(2)],coo(3);
  //  plg,coo(3)*cos(s)/2+coo(2),coo(3)*sin(s)/2+coo(1),width=3;
  limits,0.,1.,0.,1.;
}
  
