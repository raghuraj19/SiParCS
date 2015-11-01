#The following script computes how long mKrig takes for approximately 10000 observations
#using the CO2 dataset from fields.  Accelerated using 1 GPU. Uses chordal distances.

#load accelerated fields
library('fields')
thisDir = "~/git/fieldsmagma_all/fieldsMAGMA/examples"
source(paste0(thisDir, "/../fieldsMAGMA.r"))
thisDir = "~/git/fieldsmagma_all/fieldsMAGMA/examples"

#this function performs spatial analysis for given lon/lat coordinates (lonlat) and observations (y)
workflow = function(x, y, nGPU) {
  #start with spatial parameter guess of theta=8
  theta <<- 8
  
  print('optimizing over lambda')
  #calculate likelihood for different lambdas using mKrig.MLE and these lambda:
  lambda = 10^seq(-5, 0, length=10)
  out = mKrigMAGMA.MLE(x, y, theta=theta, lambda=lambda, cov.args= list(Covariance="Exponential", Distance="rdist"), lambda.profile=FALSE, nGPUs=nGPU)
  lnLikes = out$summary[,2]
  lambda.MLE <<- out$lambda.MLE
  
  #use spline interpolator to find max likelihood:
  interpGrid = 10^(seq(-5, 0, length=150))
  interpLnLikes = splint(lambda, lnLikes, interpGrid)
  index = which.max(interpLnLikes)
  
  #compute MLE krig:
  out.mle <- mKrigMAGMA(x, y, theta=theta, lambda=lambda.MLE, cov.args= list(Covariance="Exponential", Distance="rdist"), nGPUs=nGPU)
  
  #Now predict surfaces and exact errors using Krig and mKrig objects: 
  print('computing Kriging surfaces...')
  out.p <<- predictSurface(out.mle)
  
  #compute error:
  print('approximating Kriging error')
  minLat = min(x[,2])
  maxLat = max(x[,2])
  size = maxLat - minLat
  gridExpansion = max(c(1 + 1e-07, 10*theta/size))

  #make sure grid is at least 10 times as big as theta to avoid error
  set.seed(1)
  out.sim <<- try(sim.mKrig.approx(out.mle, M=5, gridExpansion=gridExpansion))
  if(class(out.sim) == 'try-error') {
    print('Used extra large gridExpansion')
    out.sim <<- sim.mKrig.approx(out.mle, M=5, gridExpansion=2*gridExpansion)
    print('Finished using extra large grid expansion')
  }
  
  invisible(NULL)
}

#now time the workflow for Krig or mKrig using CO2 dataset from fields:
data(CO2)

#filter dataset by lon/lat coords using "lim": only use a subset of 
#the data that is close enough to(lon=0,lat=0)
lim <- 105
ind <- (-lim < CO2$lon.lat[,1]) & (CO2$lon.lat[,1] < lim)
ind <- (ind & -lim/2 < CO2$lon.lat[,2]) & (CO2$lon.lat[,2] < lim/2)
x = CO2$lon.lat[ind,]
y <- CO2$y[ind]
n <- length(y)

#record timings (note that there may be some GPU startup cost in the accelerated timing):
print(paste0('Using ', n, ' data points with lim: ', lim))
time = system.time(workflow(x, y, 0))[3]
print('finished non-accelerated workflow, beginning accelerated workflow')
accelTime = system.time(workflow(x, y, 1))[3]

print(paste0('non-accelerated time: ', time))
print(paste0('accelerated time: ', accelTime))

#plot prediction surface
png(height=5, width=7, units="in", res=300, file="CO2_prediction.png")
surface(out.p, xlab="Longitude", ylab='Latitude', main="CO2 Concentration")
world(add=TRUE)
points(x, pch=20, cex=.1)
dev.off()

#plot estimated error surface
png(height=5, width=7, units="in", res=300, file="CO2_SE_est.png")
surSE<- apply( out.sim$Ensemble, 1, sd )
print(paste0("out.sim names: ", names(out.sim)))
image.plot(as.surface(out.sim$predictionPoints, surSE), xlab='Longitude', ylab='Latitude', main=sprintf('Kriging Estimated Standard Error for Lambda = %.1f, Theta = %.1f', lambda.MLE, theta))
world(add=TRUE)
points(x, pch=20, cex=.1)
dev.off()

#Smooth estimated error surface using thin plate spline
png(height=5, width=7, units='in', res=300, file="CO2_SE_smooth.png")
smoothSE = Tps(out.sim$predictionPoints, surSE)
smoothSE.surf = predictSurface(smoothSE, nx=100, ny=100)
surface(smoothSE.surf, xlab='Longitude', ylab='Latitude', main=sprintf('Smoothed, Estimated Standard Error for Lambda = %.1f, Theta = %.0f', lambda.MLE, theta))
world(add=TRUE)
points(x, pch=20, cex=.1)
dev.off()

#save data:
save(list=c("n", "time", "accelTime", "lim"), file=paste0(thisDir, "/mKrigWorkflowComparison1_n", n, ".RData"))

