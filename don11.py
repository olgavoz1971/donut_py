import idlwrap    # https://r4lv.github.io/idlwrap/index.html
import ztools
from dist import dist
import common
from numpy import *

# Measuring low-order aberrations from defocused images
# May 3, 2006
# A.A.Tokovinin
# ovoz: IDL to python3 translation
#------------------------------------------------

#from collections import namedtuple
from recordtype import recordtype # Mytable !

def read_params(filename):
  di = {}
  with open(filename) as f:
    for line in f:
      if(line[0] == '#'):
        continue
# remove stuff after comment-symbol ';' or '#'
      line = line.split('#')[0].split(';')[0]
      li = line.split(':')
      if len(li) > 1:
        key   = li[0].strip().lower()
# tuning:
        if key == 'lambda': key = 'lambd'
        if key == 'd': key = 'D'
        value = li[1].strip().replace(',','').replace('}','')
# convert to float if possible:
        try:
          value1 = float(value)
          if value1.is_integer(): value1 = int(value1)
        except ValueError:
          value1 = value
        di[key] = value1
  return di

# Pre-compute the parameters and save them in the COMMOM block
def init(filename_params):
   di = read_params(filename_params)
#   Donpar = namedtuple('Donpar', di)
   Donpar = recordtype('Donpar', di)   # Mutable one!
   common.donpar = Donpar(**di)

   donpar = common.donpar
#
   common.ngrid = donpar.ngrid
   D            = donpar.D
   eps          = donpar.eps
   common.lambd = donpar.lambd
   common.ccdpix= donpar.pixel
   common.sflag = 0 # fast seeing blur, set to 1 for slow calc.

#   Find reasonable limits on the grid parameters
   common.asperpix =  206265.*(1e-6 * common.lambd)/D # maximum fine-pixel size
# 
   k = idlwrap.floor( log10( common.ccdpix/common.asperpix ) / log10(2.0) ) + 1
   common.npixperpix = 2**k
   common.fovpix = int(2 * common.ngrid / common.npixperpix)     # CCD field size
   common.asperpix = common.ccdpix / common.npixperpix
   size = 206266.*(1e-6*common.lambd) / common.asperpix
   common.Rpix = common.ngrid / size * D
#
   print('Rebinning factor: ', common.npixperpix)
   print('Grid pixel: ',common.asperpix,' arcsec')
   print('Grid size: ', size,' m')
   print('CCD pixel:  ',common.ccdpix,' arcsec')
   print('CCD field:  ',common.fovpix*common.ccdpix,' arcsec')
   print('CCD format: ', common.fovpix)
   print('Rpix:', common.Rpix)
#
   ngrid = common.ngrid
   Rpix  = common.Rpix
   common.r = idlwrap.shift( dist(2*ngrid), ngrid, ngrid ) # distance from grid center, pixs
   common.inside_mask = logical_and((common.r <= Rpix), (common.r >= Rpix*eps))
   common.pupil = idlwrap.fltarr(2*ngrid, 2*ngrid)    # zero array
   common.pupil[common.inside_mask] = 1
   n = common.pupil[common.inside_mask].size
#
   x_ = idlwrap.findgen(2*ngrid) - ngrid
   x__ = repeat(1.,2*ngrid) 
   x = idlwrap.operator_(x_, "#", x__)
#
   theta = arctan2(transpose(x),x)
   theta[ngrid, ngrid] = 0.
   common.zgrid = idlwrap.fltarr(n, 2)
   common.zgrid[0,:] = common.r[common.inside_mask] / Rpix
   common.zgrid[1,:] = theta[common.inside_mask]

#  -------------------------------------------------------
from rebin import rebin
def getimage(z):
# z is the Zernike vector in microns, starting from Z=2 (tip)
# z[0] is seeing in arcseconds

  donpar   = common.donpar
  ngrid    = common.ngrid
  asperpix = common.asperpix
  fovpix   = common.fovpix
  ccdpix   = common.ccdpix
  
  fact = 2.0 * pi / donpar.lambd
  nzer = z.size
  phase = common.zgrid[0,:] * 0.  # empty array for phase
  for j in range(1, nzer):
     phase += fact * z[j] * ztools.zernike_estim (j+1, common.zgrid )

  tmp = idlwrap.fltarr(ngrid*2, ngrid*2 )
  common.uampl = idlwrap.dcomplex(tmp, tmp)

  common.uampl[common.inside_mask] = idlwrap.dcomplex(cos(phase), sin(phase))
  common.seeing = z[0]

#---------  compute the image ----------------------
  imh = abs( idlwrap.shift( idlwrap.fft( idlwrap.shift( common.uampl, ngrid, ngrid ), inverse=True), ngrid, ngrid ) )**2
  if common.sflag > 0: # exact seeing blur
    common.filter2 = exp( -2.0 * pi**2 * (common.seeing/2.35/asperpix/2/ngrid)**2 * r**2) # unbinned image
    imh = abs( idlwrap.fft( idlwrap.shift( fft(imh), ngrid, ngrid) * common.filter2 ) )
    impix = rebin( imh, fovpix, fovpix ) # rebinning into CCD pixels
  else:
    rr = idlwrap.shift( dist(fovpix), int(fovpix/2), int(fovpix/2) ) 
    common.filter2 = exp(-2. * pi**2 * (common.seeing/2.35/ccdpix/fovpix)**2 * rr**2) # binned image
    impix = rebin( imh, fovpix, fovpix) # rebinning into CCD pixels
    impix = abs( idlwrap.fft( idlwrap.shift( idlwrap.fft(impix), fovpix/2, fovpix/2 ) * common.filter2)) # Seeing blur
# show it !!! jk !!!
#  tvscl, congrid(impix,2*ngrid,2*ngrid)
  return impix/idlwrap.total(impix)  

def newimage(a, jzer):
# a is the amplitude change of aberration (micron), Jzer is the Zernike number 
# (1=seeing, 4=focus etc.)
#
  donpar   = common.donpar
  zgrid    = common.zgrid
  asperpix = common.asperpix
  ngrid    = common.ngrid
  inside_mask = common.inside_mask
  ccdpix   = common.ccdpix
  fovpix   = common.fovpix
  filter2  = common.filter2
  sflag    = common.sflag
  seeing   = common.seeing
  newampl  = common.uampl.copy()

  if jzer > 1:      # Change Zernike coefficient 
    newphase =  2.0 * pi / donpar.lambd * a * ztools.zernike_estim(jzer, zgrid)  
    newampl[inside_mask] *= idlwrap.dcomplex(cos(newphase), sin(newphase))
    filter = filter2.copy()
  else:     # new seeing
    newseeing = seeing + a
    if sflag > 0:
       filter = exp(-2.0 * pi**2 * (newseeing/2.35/asperpix/2/ngrid)**2 * r**2) # unbinned image
    else:
      rr = idlwrap.shift( dist(fovpix), int(fovpix/2), int(fovpix/2)) 
      filter = exp(-2.0 * pi**2 * (newseeing/2.35/ccdpix/fovpix)**2 * rr**2)    # binned image
#
# ---------  compute the image ----------------------
#
  imh = abs(idlwrap.shift( idlwrap.fft( idlwrap.shift( newampl, ngrid, ngrid), inverse=True), ngrid, ngrid ) )**2
  if sflag > 0:    # exact seing blur
    imh   = abs( idlwrap.fft( idlwrap.shift( idlwrap.fft(imh), ngrid, ngrid)*filter))
    impix = rebin(imh, fovpix, fovpix)  # rebinning into CCD pixels
  else:
    impix = rebin(imh, fovpix, fovpix)  # rebinning into CCD pixels
    impix = abs( idlwrap.fft( idlwrap.shift( idlwrap.fft(impix), fovpix/2, fovpix/2) * filter)) # Seeing blur
# !!! draw it JK !!!
#  tvsc#l, congrid(impix,2*ngrid,2*ngrid)
  return impix/idlwrap.total(impix)  
#
# -------------------------------------------------------
#
# --- tmp ----:
import matplotlib.pyplot as plt
from astropy.visualization import SqrtStretch, simple_norm
# ---
def getmom(impix1):
#
# returns the vector of Zernike aberrations in microns
#
  donpar = common.donpar

  n = common.ngrid / common.npixperpix
# ---- original toko code
#  xx_  = idlwrap.findgen(2*n) - n
#  xx__ = repeat(1, 2*n)
#  xx   = idlwrap.operator_(xx_, "#", xx__)
#  yy   = transpose(xx)
# ----
  impix  = impix1.copy()
  thresh = impix.max() * donpar.thresh     # take third value from max instead!
  impix[impix<= thresh] = 0
  imh0   = impix.sum()
  print('imh0 = ', imh0, 'n =',n)
# --- search moments:
  from skimage.measure import moments, moments_central
  image = impix
  M =  moments(image)
  Mc = moments_central(image)
  print('xc =', M[0,1]/M[0,0]-n-0.5, 'yc =', M[1,0]/M[0,0]-n-0.5)
# ------------
#
#  xc  = idlwrap.total( xx * impix ) / imh0
#  yc  = idlwrap.total( yy * impix ) / imh0
#  print('xc, yc =', xc, yc)
#  mxx = idlwrap.total( impix * (xx - xc)**2 ) / imh0
#  myy = idlwrap.total( impix * (yy - yc)**2 ) / imh0
#  mxy = idlwrap.total( impix * (xx - xc) * (yy - yc) ) / imh0
#  print('mxx, myy, mxy =', mxx, myy, mxy)
#  print('               ', Mc[0,2],Mc[2,0], Mc[1,1])
  mxx = Mc[0,2]
  myy = Mc[2,0]
  mxy = Mc[1,1]
# Decenter:
  xc = M[0,1]/M[0,0] - n - 0.5   # ovoz
  yc = M[1,0]/M[0,0] - n - 0.5
#
  scale = common.npixperpix / (common.ngrid / common.Rpix)
#
#  a2 = scale * (xc + 0.5) * pi * 0.5     # toko
#  a3 = scale * (yc + 0.5) * pi * 0.5
  a2 = scale * (xc) * pi * 0.5
  a3 = scale * (yc) * pi * 0.5
  a4 = scale * sqrt( (mxx + myy) * 0.5 ) / 1.102
# MAX does not work with arrays
  a4 = sqrt( idlwrap.operator_(a4**2 - (0.5/2.35)**2, ">",0))
  a5 = scale * mxy * ( mxx * myy )**(-0.25) / 1.45
  a6 = scale * (mxx - myy) * 0.5 * (mxx * myy)**(-0.25)/1.45
  zestim = array([0., a2, a3, a4, a5, a6]) * donpar.lambd / (2.0 * pi) # estimated Zernike aberrations
#  zestim[0] = 0.5    # JK !!!! take it from params now !!! 
  zestim[0] = donpar.seeing    # JK !!!! take it from params now !!! 
  print('zestim = ', zestim)
  print('\n\n\n')
#
# tmp --------
#  imgtmp = impix1
#  norm = simple_norm(imgtmp, 'sqrt', percent=99)
#  plt.imshow(imgtmp, norm=norm, origin='lower', cmap='Greys_r',
#           interpolation='nearest')
#  print('Check-stop\nClose graph-window to continue...')
#  plt.show()
# -------------
  return zestim

#from congrid import congrid
#         tmp[:,256:512] = congrid(model, (256,256))

# -------------------------------------------------------
def find(impix, zres, nzer, display):
#
# Nzer is the highest Zernike number, zres is the  result
#
# ovoz display = display method from GUI
  nzer = max(nzer,6)
  n1 = zres.size
  nzer1 = max(n1, nzer)
  z0 = idlwrap.fltarr(nzer1)
  z0[:n1] = zres.copy()
#
  impixnorm = impix / idlwrap.total(impix)
  impixmax = impix.max()
#
  xi = idlwrap.fltarr(nzer)
  for j in range(1, nzer):
      xi[j] = common.donpar.lambd / (2.0 * pi) * 0.5 / int( int( sqrt( 8*(j+1)-6 )-1 ) / 2 ) # 0.5/n amplitudes
## for j=1,nzer-1 do xi[j] = donpar.lambda/(2.*!pi)*0.5/(long(sqrt(8L*(j+1L)-6L)-1L)/2L) ; 0.5/n amplitudes
  xi[0] = 0.1

  indonut_mask = impixnorm > impixmax * common.donpar.thresh
  im = impixnorm[impixnorm > impixmax * common.donpar.thresh]
  n = im.size
  chi2old = n**2
  chi2norm = idlwrap.total( im**2 )
#
  ncycle = 20
  thresh0 = 0.01    # initial SVD threshold
  norm = impix.max()
  thresh = thresh0  # SVD inversion threshold, initial
#  print('Z  ', idlwrap.indgen(nzer) + 1)    # format='(A4,100I8)'  
  print('Z  ', end='')
  for i in idlwrap.indgen(nzer) + 1: print('%8i' % i, end='')
  print()

  lambd = 1. # for L-M method

  for k in range(0, ncycle):
#    print('zres start cycle = ', zres)
    model = getimage(z0)
    im0   = model[indonut_mask]
    chi2  = sqrt( idlwrap.total( (im0 - im)**2 ) / chi2norm )
    print('\n\n ---------------- Cycle:', k+1, 'RMS = ', chi2*100, 'percent')
    print('um:', end='')
    for z in z0[:nzer]: print('%8.3f' % z, end='')
    print()
#    with  printoptions(precision=3, suppress=False): print(z0[:nzer])  #, format='(A4,100F8.3)'

    thresh = thresh * 0.5 # > 1e-3

    if chi2 < 5e-5:
      print('chi2 < 1e-4', chi2)
      break

    if chi2 <= chi2old:
      zres  = z0.copy() 
      lambd = lambd * 0.1
      if (chi2 >= chi2old * 0.99) and (k > 3):
        print('Stop on chi2 >= chi2old*0.99', chi2, chi2old, chi2old*0.99);
        break
      chi2old = chi2 
    else:
      z0 = zres.copy()
      thresh = thresh0
      lambd = lambd * 10.0
      print('Divergence... Now LM parameter =', lambd)
 
    if k % 2 == 0:
      imat = idlwrap.fltarr(n, nzer)
      print('Computing the interaction matrix...')

      for j in range(0, nzer):
        imat[j,:] = ( (newimage(xi[j], j+1))[indonut_mask] - im0 ) / xi[j]
      tmat = idlwrap.operator_( transpose(imat), "#", imat )
      tmp = ztools.svd_invert(tmat, thresh)
      invmat = idlwrap.operator_(tmp, "#", transpose(imat))

    dif = im - im0
    dz = matmul(dif, invmat)  # idlwrap.operator_# does not work with these shapes
#   dz = invmat # dif  
    z0[:nzer] += 0.7 * dz
    z0[0] = max(z0[0],  0.2)
    z0[0] = min(z0[0], 1.5)

    d1 = dif.min()
    d2 = dif.max()
#
# ---------------   display the image (left: input, right: model)
    if display is not None:
      display(model)

#    ny, nx = impix.shape
#    tmp = zeros((ny, nx*2))
#    tmp[:,:nx] = impix
#    tmp[:,nx:] = model
#    imgtmp = tmp
#    norm = simple_norm(imgtmp, 'sqrt', percent=99)
#    plt.imshow(imgtmp, norm=norm, origin='lower', cmap='Greys_r',
#        interpolation='nearest')
#    print('Check-stop\nClose graph-window to continue...')
#    plt.pause(0.1)
  print('Fitting done!')
#  plt.show()
# -------------
  return zres, model, chi2

# -------------------------------------------------------
def fit(impix, display):
#
#  display is GUI-method
# preliminary and final fitting
  zres = getmom(impix)
#  print(' +++++++++++++++++++++++  getmom: zres =', zres)
# to suit original DONUT example
#  zres[:6] = [0.500000,     -1.60505,     -1.84489,      2.11016,    -0.100961,     0.140742]
#  print(' +++++++++++++++++++++++  getmom: zres =', zres)


  nzer = common.donpar.nzer
#  if common.donpar.static != 0:
#     readz(z0, donpar.static )
#  else:
#     z0 = idlwrap.fltarr(nzer)
#
  z0 = idlwrap.fltarr(nzer)
  z0[:6] = zres
# z0[0:5] = zres  IDL
  if common.donpar.efoc < 0: z0[3:6] *= -1.
# if (efoc lt 0) then z0[3:5] *= -1.
  zres = z0
# 
  zres, immod, chi2 = find(impix, zres, nzer, display)
#
  return zres, immod, chi2

# -------------------------------------------------------
def writepar(filename):
# Write parameters into a file
#
  outfile = open(filename, 'w')  
  for key, value in zip(donpar._fields, donpar):
    outfile.write(key + ' : ' + value)
  outfile.close()
  print('Parameters are saved in', filename)

# ------------------------ 
#pro savez, z, file
#; Save Zernike vector in ASCII file
#
# openw, unit, file,/get_lun
# fmt='(F8.3,1X,I2)' 
# for i=0,n_elements(z)-1 do printf, unit, z[i], i+1, format=fmt
# close, unit 
# free_lun, unit 
#  print, 'Zernike vector is saved in ',file 
#end
#
# ------------------------------------------------
#pro readz, z, file
# Read Zernike vector from ASCII file
#
#  res = findfile(file, count=c)
#  if (c eq 0) then begin
#    print, 'File ',file,' is not found, exiting'
#    return
#  endif
#
# openr, unit, file,/get_lun
# tmp = idlwrap.fltarr(300)
# i = 0 & x = 0.0
# while not eof(unit) do begin
#   readf, unit, x
#   tmp[i] = x
#   i +=1
# endwhile
# close, unit 
# free_lun, unit 
# if (i eq 0) then z=-1 else z = tmp[0:i-1]
# print, i,' Zernike coefficients are retrieved from ',file 
#
#end
#;------------------------------------------------
#pro saveres, resfile, z, chi2, imfile
#
#@donut.common
#
#  openw, unit, resfile,/get_lun, /append
#  printf, unit, imfile, donpar.xc, donpar.yc, flux, chi2, z,  $
#format='(A20,2I6,E10.3,F8.4,100F8.3)'
#  close, unit 
#  free_lun, unit 
#  print, 'Results are saved!'
#
#end
# ------------------------------------------------
def extract(img, xc, yc, nccd):
# image extraction

  ix1 = max(xc-nccd, 0)
  ix2 = min(xc+nccd, img[0,:].size)# - 1
  iy1 = max(yc-nccd, 0)
  iy2 = min(yc+nccd, img[:,0].size)# - 1

  img1 = img[iy1:iy2, ix1:ix2]    # cut out the required part
  img1 = img1 - img1.min()        # subtract background

# Alternative skimage method:   Works better!
  from skimage.measure import regionprops, moments
#
# Find strong background:
  idx = int(round(0.3 * img1.size))
  sorted_idx = argsort(img1.ravel())
  i = sorted_idx[idx]
  threshold_value = (img1.ravel())[i] # 30% quantile of pixel distribution
  print('threshold_value =', threshold_value)
  labeled_foreground = (img1 > threshold_value).astype(int)
  properties = regionprops(labeled_foreground, img1)

# My algorithm (ovoz):
  ix1 = round(properties[0].weighted_centroid[1] - int(nccd/2) +0.5)
  ix2 = ix1 + nccd
  iy1 = round(properties[0].weighted_centroid[0] - int(nccd/2) +0.5)
  iy2 = iy1 + nccd

  if (ix1 < 0) or (iy1 < 0) or ((ix2-ix1) < nccd) or ((iy2-iy1) < nccd):
      print('Image is cut on one side!')
      return(-1)
#
  impix = img1[iy1:iy2, ix1:ix2]
#  impix = img1[9:25,9:25]  # to suit original donut
  impix_ravel = impix.ravel()
  idx = idlwrap.fix(0.1 * nccd**2)
  sorted_idx = argsort(impix_ravel)
  i = sorted_idx[idx]
  backgr = impix_ravel[i] # 10% quantile of pixel distribution
  print('backgr =', backgr)

  impix = impix - backgr
  common.flux = idlwrap.total(impix)
  print('Total flux, ADU:', common.flux)
  impix = impix / common.flux
# max does not work with arrays !!!
  common.sigpix = idlwrap.operator_(impix, ">", 0) * common.flux * common.donpar.eadu \
       + common.donpar.ron**2  # variance in each pixel 

  return impix
# ------ end of the program ----------------------------  
