# PY_ANTENNA

by Mikhail Kravchenko
mkravche@vols.utk.edu
-----------------

A library in c++ and python that deals with RF computation.

This file project the classes that define the EM fields. 

This is me consolidating (or attempting to consoldiate) and documenting various code and notes I have lazily accreted in the last few years. I am historically very bad at showing my work perhaps because I am too self concious about my implmentation. I will tend to publish work that I feel is 'good enough' to be showable. I would rather code I push to be readable and usuable by other people. The code I keep for myself needs to be better documented, optimized, and in some cases moved from paper to computer.

This code is not particularly valuable feel free to take and use my examples if it helps you. If you want to give credit, feel free.

Included files:

README.md   &nbsp;&nbsp;&nbsp;&nbsp; (this file)

LICENSE     &nbsp;&nbsp;&nbsp;&nbsp;        (MIT License)

T_junction.bmp  &nbsp;&nbsp;&nbsp;&nbsp;    Input coordinates defined by a raster.

Output3.webm    &nbsp;&nbsp;&nbsp;&nbsp;    The pretty output. I'm honestly doing this for the aestetic.

T_junction.py &nbsp;&nbsp;&nbsp;&nbsp; This code extracts the relevant coordinates from the BMP image. The black pixels are the wall and are detected by finding pixels with B=0
The 'pink' pixels are the source coordinates and are for pixels with G=127. The BMP has BGR pixel orientation.




