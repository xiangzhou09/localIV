rm(list=ls())

# extract Y, X, D, and Z

nlsy <- read.table("C:/Users/Xiang/Dropbox/MTE/NLSY data analysis/nlsyfinal.csv",sep=",",header=T)

formula_Y <- wage~exp+d57+d58+d59+d60+d61+d62+d63+expsq+urban14+numsibs+I(numsibs^2)+cafqt+I(cafqt^2)+
  mhgc+I(mhgc^2)+lavlocwage17+I(lavlocwage17^2)+avurate+I(avurate^2)+lwage5+lurate

formula_D <- state~urban14+numsibs+I(numsibs^2)+mhgc+I(mhgc^2)+cafqt+I(cafqt^2)+
  d57+d58+d59+d60+d61+d62+d63+lavlocwage17+I(lavlocwage17^2)+avurate+I(avurate^2)+
  pub4+pub4*cafqt+pub4*mhgc+pub4*numsibs+lwage5_17+lwage5_17*cafqt+lwage5_17*mhgc+
  lwage5_17*numsibs+lurate_17+lurate_17*cafqt+lurate_17*mhgc+lurate_17*numsibs+
  tuit4c+tuit4c*cafqt+tuit4c*mhgc+tuit4c*numsibs

mte_fit <- mte(formula_D, formula_Y, data = nlsy)

