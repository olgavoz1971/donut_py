# LICENCE:
#; Copyright (C) 2005, SJT, P. Chanial
# This program is free software; you can redistribute it and/or modify  
# it under the terms of the GNU General Public License as published by  
# the Free Software Foundation; either version 2 of the License, or     
# (at your option) any later version.                                   
# manually converted FROM GDL dist.pro by ovoz
#

import idlwrap
from numpy import sqrt
def dist(m,n=None):
  if n is None: n = m
  if (m <= 0) or (n <= 0):
     print("Array dimensions must be greater than 0.")
     return -1
  x = idlwrap.findgen(m)
  x = idlwrap.operator_(x, "<", (m-x))**2
  result = idlwrap.fltarr(m,n)
  result[0,:] = sqrt(x)
#
  for i in range(1, int(n/2)+1):
     dist = sqrt(x + i**2.0)
     result[i,:]   = dist
     result[n-i,:] = dist
  return result
