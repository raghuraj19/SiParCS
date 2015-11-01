#load fields
#library(fields)
library(fields)

#load all MAGMA-accelerated functions
thisDir = "/glade/u/home/jpaige/git/fieldsmagma_all/fieldsMAGMA"

dyn.load(paste0(thisDir, '/magmaCholesky/magmaCholesky.so'))
#dyn.load(paste0(thisDir, '/magmaCholesky/magmaCholesky_m.so'))
#dyn.load(paste0(thisDir, '/magmaCholesky/smagmaCholesky.so'))
#dyn.load(paste0(thisDir, '/magmaCholesky/smagmaCholesky_m.so'))
#dyn.load(paste0(thisDir, '/magmaCholesky/magmaCholeskyNew_m.so'))
#dyn.load(paste0(thisDir, '/magmaCholesky/smagmaCholeskyNew_m.so'))

source(paste0(thisDir, '/RwrapperMAGMA.r'), chdir=TRUE)
source(paste0(thisDir, '/accelMKrig/mKrigMAGMA.r'), chdir=TRUE)
source(paste0(thisDir, '/accelMKrig.MLE/mKrigMAGMA.MLE.r'), chdir=TRUE)
