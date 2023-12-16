import cv2 
from matplotlib import pyplot as plt 
  
img = cv2.imread('INSERT_IMAGE_DIR_HERE.jpg',0) 
  
# find frequency of pixels in range 0-255 
histr = cv2.calcHist([img],[0],None,[256],[0,256]) 
  
plt.hist(img.ravel(),256,[2,256]) 
plt.show() 