# YAO package

Copyright (c) Francois Rigaut (frigaut) 2002-2016.
Initial release June 2002.

## Developers
Francois Rigaut         frigaut@gmail.com  
Marcos van Dam          marcos@flatwavefronts.com  
Ralf Flicker            †July 2009  
Damien Gratadour        damien.gratadour@obspm.fr  
Aurea Garcia-Rissmann   agrissmann@gmail.com

## Synopsis
Yao is a Monte-Carlo simulation tool for Adaptive optics (AO) systems, written as a yorick plugin. It uses a number of custom developed functions to simulate wavefront sensors (WFS), deformable mirrors (DM) and many other aspects of an AO loop. The core functions are written in C, hence are very fast.

## Help pages
http://frigaut.github.com/yao/index.html


## Required Dependencies

- FFTW (http://www.fftw.org)
- yorick-yutils (https://github.com/frigaut/yorick-yutils)
- yorick-imutil (https://github.com/frigaut/yorick-imutil)
- yorick-spydr (https://github.com/frigaut/yorick-spydr)
- yorick-usleep (https://github.com/frigaut/yorick-usleep)

## Optional Dependencies
- yorick-soy (https://github.com/frigaut/yorick-soy)
- ygsl (https://github.com/emmt/ygsl)
- ylapack (available from Eric Thiebaut)