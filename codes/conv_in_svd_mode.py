# -*- coding: utf-8 -*-
#/usr/bin/python

## Illustration of a convolution with Singular Value Decomposition
## _______________________________________________________________
from numpy.linalg import svd
import numpy as np

def svd_convolution():
    I = 
  #import sys,Image,numpy
  #imagefile = "original_image.jpg"	
  origImage = Image.open(imagefile)
  colors = origImage.split()
  if(len(colors) != 1):
    print "Only one channel is handled in the script"
  
  (img_width, img_height) = origImage.size

  reshaped_colors = numpy.reshape(colors[0], (img_height, img_width))
  
  # Gaussian kernel
  # Uselessly large
  h_filt = 128
  w_filt = 128
  X,Y = numpy.mgrid[0:h_filt,0:w_filt]
  sigma_h = 2.0
  sigma_w = 4.0
  kernel = 1.0/(sigma_h * sigma_w * 2.0 * numpy.pi) * numpy.exp(-0.5*(((X-(h_filt-1.0)/2.0)/sigma_h)**2 + ((Y-(w_filt-1.0)/2.0)/sigma_w)**2))

  (u,s,vh) = svd(kernel)
  
  # Compute the rank of the matrix
  rank = 0
  cmp_threshold = 1e-10
  for i in range(len(s)):
    if(s[i] >= cmp_threshold):
      rank = rank + 1
  print "Rank of the kernel :", rank
  
  resultImage = numpy.zeros(reshaped_colors.shape)
  imgTemp = numpy.zeros(reshaped_colors.shape)

  for i in range(rank):
    # First pass
    for j in range(reshaped_colors.shape[1]):
      imgTemp[:,j] = numpy.convolve(reshaped_colors[:,j], u[:,i],'same')
    for j in range(reshaped_colors.shape[0]):
      resultImage[j,:] += s[i] * numpy.convolve(imgTemp[j,:], vh[i,:],'same')

  image_file = Image.new("L", (img_width, img_height))
  image_file.putdata(resultImage.reshape(resultImage.size))
  image_file.save("result_image.jpg")
  

if __name__=='__main__':
  from timeit import Timer
  t = Timer("svd_convolution()","from __main__ import svd_convolution")
  nb_repetitions = 1
  print "%.4f seconds" % (t.timeit(nb_repetitions)/nb_repetitions)
