##Author: Joel Chacón c.
##Soruces: cat.pgm; PGM version-->P2 
##Example for reading PGM files 
##More information bout PGM format: https://www.usna.edu/Users/cs/nchamber/courses/si204/s18/lab/l08/lab.html
import numpy as np
from libpgm import pgmread, pgmwrite


def normalize(image):
  maxv = max(map(max, image)) 
  minv = min(map(min, image)) 
  return (255*(image-minv)/(maxv-minv)).astype(np.int32)
def factorial(n):
  if n <= 1 :
   return 1
  return n*factorial(n-1)
def forward_table(idx, n, row): ##Pag 380 Rajesh
  if n==0:
    return row[idx];
  return forward_table(idx+1, n-1, row)-forward_table(idx, n-1, row);

def backward_table(idx, n, row): ##Pag 380 Rajesh
  if n==0:
    return row[idx];
  return backward_table(idx, n-1, row)-backward_table(idx-1, n-1, row);

def central_table(idx, n, row): ##Pag 381 Rajesh homework!!!!!
  return 0

def forward_difference_by_cols_2(original): ##The most simple example two points
  nrows, ncols= np.shape(original)
  gx = np.full_like(dataimage[0],0) ##make room for gx with the same dimensions that the original image with zero values...
  ##Apply forward difference to rows..
  ##notes: h is always one..
  for i in range(nrows):
   for j in range(ncols-1):
     f_x_h_1 = (original[i][j+1])
     f_x_h = (original[i][j])
     gx[i][j] = (f_x_h_1-f_x_h)
  return gx

def forward_difference_by_cols(original, npoints, fac): ##Pag 498 Rajesh  Y-gradient
  nrows, ncols= np.shape(original)
  gx = np.full_like(dataimage[0],0) ##make room for gx with the same dimensions that the original image with zero values...
  ##Apply forward difference to rows..
  ##notes: h is always one..
  for i in range(nrows):
   for j in range(ncols-npoints):
      for k in range(npoints-1): ##pow --> npoints-1
        gx[i][j] += (1.0/fac[k+1])*forward_table(j, k+1, original[i])
  return gx

def forward_difference_by_rows(original, npoints, fac): ##Pag 498 Rajesh X-gradient
  nrows, ncols= np.shape(original)
  gx = np.full_like(dataimage[0],0) ##make room for gx with the same dimensions that the original image with zero values...
  ##Apply forward difference to rows..
  ##notes: h is always one..
  for i in range(nrows-npoints):
   for j in range(ncols):
      for k in range(npoints-1): ##pow --> npoints-1
        gx[i][j] += (1.0/fac[k+1])*forward_table(i, k+1, original[:, j])
  return gx

fac = np.zeros(10)
for k in range(10):
 fac[k]=factorial(k) ##this is to avoid recomputation of factorial...

#################simple example...
#dataimage = pgmread('shark.pgm') ## information ---> [0]:array, [1], width, [2] height 
#original = np.copy(dataimage[0])
#gy = forward_difference_by_cols_2(original)
#gy = forward_difference_by_cols(original, 5, fac)
#pgmwrite(gy, "out.pgm")  ##export pgm..

##########################MAIN CODE!!!
dataimage = pgmread('shark.pgm') ## information ---> [0]:array, [1], width, [2] height 
original = np.copy(dataimage[0])
for npoints in [2, 3, 5, 7, 10]:
 gy = forward_difference_by_rows(original, npoints, fac)
 gy = normalize(gy)##Important!!!
 pgmwrite(gy, "gy_points_"+str(npoints) +".pgm")  ##export pgm..

 gx = forward_difference_by_cols(original, npoints, fac)
 gx = normalize(gx)##Important!!!
 pgmwrite(gx, "gx_points_"+str(npoints) +".pgm")  ##export pgm..




