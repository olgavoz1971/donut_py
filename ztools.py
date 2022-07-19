from numpy import *

# Toolbox of Zernike polynomials
#--------------------------------------------------
# Compute the  Noll covariance matrix of Zernike polynomials
from math import gamma
def cova_zern1(jmax):
   c = zeros([jmax - 1, jmax - 1])
   for j in arange(2, (jmax+1)):
      for jp in arange(2, (jmax+1)):
#         print(j,jp)
         n,m   = zern_num(j)
         np,mp = zern_num(jp)
#k_zz = 7.2e-3 *(!Dpi)^(8./3.)* gamma(14D0/3D0) ; gives 2.242434, c_11 = 0.44815
         k_zz = 2.242434        
         K = k_zz * (-1) ** ((n + np - 2.0 * m) / 2.0) * sqrt((n + 1) * (np + 1))
#         print(n, np, m, K)
         if (m == mp) and ((j*jp/2.0 != int(j*jp/2.0)) or ((j/2.0 == int(j/2.0)) and \
                                                      (jp/2.0 == int(jp/2.0))) or (m == 0)) :
              c[j-2,jp-2] = K * gamma((n+np-5.0/3.0)/2.0) /(gamma((n-np+17.0/3.0)/2.0) * \
                                gamma((np-n+17.0/3.0)/2.0) * gamma((n+np+23.0/3.0)/2.0))
         else:
           c[j-2,jp-2] = 0.0
   return c

#-------------------------------------------------------------------
def zern_num(num, m=None, n=None, info=None):
#  Description
#  ZERN_NUM computes the azimuthal degree M and the radial order N
#  for the sequential Zernike polynomial NUM
#  INFO  =  result display
   n_params = 1
   _opt = (m, n, info)
   if num < 1:
       print("ERROR: NUM must be an integer greater or equal 1") 
   j = int(num)
   n = int(int( sqrt(8*j-7)-1 )/2)
#   n = int((sqrt(8*j-7)-1)/2)
   if n % 2:                     # odd n
      m = 1 + 2 * int( int(j - 1 - n*(n+1)/2) / 2 )
   else:                            # even n
      m = 2 * int( int(j - n*(n+1)/2) / 2 )
   if (info is not None):   
      print(j, n, m)   
   return n, m

#-------------------------------------------------------------------
from math import factorial
def zernike_estim(mode, grid):
#
# Input: grid = set of points (polar co.) on which the Zernike must be evaluated
# output: vector of Zernike values on the grid
#    R=R+(-1.)^J*Fact(n-J)/(Fact(S)*Fact((n+m)/2-J)$_
#                          *Fact((n-m)/2-J))*grid(*,0)^(n-2*J)

  n,m = zern_num(mode)
  p = mode % 2
  R = 0.0
  for j in range(0, int((n-m)/2)+1):
    s = j   
    R = R + (-1.)**j * factorial(n-j) / (factorial(s) * factorial(int((n+m)/2) - j) \
                          * factorial(int((n-m)/2) - j)) * grid[0,:]**(n-2*j)
#
  if m == 0: 
    zz = sqrt(n+1.0) * R
#
  if m != 0:
    if p == 0:
      zz = sqrt(n+1.0) * sqrt(2.0) * cos( m * grid[1,:] ) * R
    if p > 0:
      zz = sqrt(n+1.0) * sqrt(2.0) * sin( m * grid[1,:] ) * R
#
  return zz


def  svd_invert(matrix, threshold):
#
# Returns the SVD-inverted matrix 
#
  u,ws,vv = linalg.svd(matrix, full_matrices=False)
  v = vv.T
  ww = ws.max()
  n = ws.size
  invw = identity(n)
  ncount = 0
# 
  for i in arange(0, n):
    if ws[i] < ww*threshold:
       invw[i,i] = 0.
       ncount = ncount + 1
    else:
       invw[i,i]=1./ws[i]
#
  print(ncount, 'singular values rejected in inversion')
  inv_matrix = matmul( matmul(v,invw), u.transpose())
  return inv_matrix

# -----------------------
#
# This function calculates the x and y derivative coefficients needed to compute 
# the derivative of the jth Zernike polynomial.
# (d/dx)Zj=SUM_j' gammax_j' Z_j'
# (d/dy)Zj=SUM_j' gammay_j' Z_j'
# gammax and gammay is the output vector gamma=[2,j]
#
# Date : 9 December 1999
# Written by Elise Viard, eviard@eso.org
#
def zern_derivx(j):
   n,m = zern_num(j)
   gam = zeros([j])
   for j2 in arange(1, j+1):
      n2,m2 = zern_num(j2)
      if (m-m2)**2 == 1:
         if (m != 0) and (m2 != 0):
            if ( (j % 2 == 0) and (j2 % 2 == 0) ) or ((j % 2 != 0) and (j2 % 2 != 0)):
               gam[j2-1] = sqrt((n+1)*(n2+1))
            else:
               gam[j2-1] = 0 
         elif ((m == 0) and (j2 % 2 == 0)):
            gam[j2-1] = sqrt(2.*(n+1)*(n2+1))
         elif ((m2 == 0) and (j % 2 == 0)):
            gam[j2-1] = sqrt(2.*(n+1)*(n2+1))
         else:
            gam[j2-1] = 0 
      else:
         gam[j2-1] = 0
   return gam


def zern_derivy(j):
   n,m = zern_num(j)
   gam = zeros([j])
   for j2 in arange(1, j+1):
      n2,m2 = zern_num(j2)
      if (m-m2)**2 == 1:
         if (m != 0) and (m2 != 0):
            if ( (j % 2 == 0) and (j2 % 2 != 0) ) or ((j2 % 2 == 0) and (j % 2 != 0)):
               if m2 == (m+1) and (j % 2 != 0):
                 sig = -1
               elif m2 == (m-1) and (j % 2 == 0):
                 sig = -1
               else:
                 sig = 1
               gam[j2-1] =  sig*sqrt((n+1)*(n2+1)) 
            else:
               gam[j2-1] = 0
         elif ((m == 0) and (j2 % 2 != 0)):
            gam[j2-1] = sqrt(2.*(n+1)*(n2+1))
         elif ((m2 == 0) and (j % 2 != 0)):
            gam[j2-1] = sqrt(2.*(n+1)*(n2+1))
         else:
            gam[j2-1] = 0 
      else:
         gam[j2-1] = 0
   return gam

def zern_deriv(j):
  gam = zeros([j,2])
  gam[:,0] = zern_derivx(j) 
  gam[:,1] = zern_derivy(j)
  return gam


# NOT USED
###ngrid = 64     # tmp, take it from donut.par
###Rpix = 100     # tmp !!!!!
#import idlwrap
#from dist import dist


#def getftzer(Jzer):
## Compute the Fourier Transform of Zernike mode
## 
#  x_ = idlwrap.findgen(2*ngrid) - ngrid
#  x__ = repeat(1.,2*ngrid) 
#  x = idlwrap.operator_(x_, "#", x__)
### x = (findgen(2*ngrid) - ngrid) # replicate(1.,2*ngrid)
#  y = transpose(x)
#  theta = arctan2(y,x)
#  theta[ngrid,ngrid] = 0.
##  
#  n,m = zern_num(Jzer)
##
#  f = idlwrap.shift(dist(2*ngrid), ngrid, ngrid) / (2*ngrid) *Rpix
#  f[ngrid,ngrid] = 1e-3
##  f =  shift(dist(2*ngrid),ngrid,ngrid)/(2*ngrid)*Rpix & f(ngrid,ngrid) = 1e-3
#  ftmod = sqrt(n+1.0) * idlwrap.beselj( 2*pi*f, n+1) / (pi*f)
##  ftmod = sqrt(n+1D0)*beselj(2*!dpi*f,n+1)/(!dpi*f)   
#  if m == 0: zz = ftmod * idlwrap.dcomplex(0, 1.0)**int(n/2)
##  IF (m EQ 0) Then  zz = ftmod*complex(0,1D0)^(n/2)
#  if m != 0:
#    if (Jzer % 2) == 0: 
#      fact = sqrt(2.0)*cos(m * theta)
#    else:
#      fact = sqrt(2.0)*sin(m * theta)
#    zz = ftmod * fact * (-1)**( (n - int(m/2)) ) * idlwrap.dcomplex(0, 1.0)**m
##    iF ((Jzer MOD 2) EQ 0) then fact=sqrt(2.)*Cos(m*theta) $ 
##                           else fact=sqrt(2.)*sin(m*theta)
##    zz = ftmod*fact*(-1)^((n-m/2))*complex(0,1D0)^m
##;    tvscl, zz
##;    tvscl, imaginary(zz), 2*ngrid, 0
##; aber = shift(fft(shift(zz,ngrid,ngrid)), ngrid, ngrid) ; get the mode shape
#    return zz
#
#
##res = getftzer(3)
##print('res =', res)
