### cd("/glade/u/home/vinayr/")
include("MRAfunctions.jl")
using Distances  # for distance computation (pairwise())

J=2; M=2; r=3;
nper=r;
nAll=nper*J^M;
domain=[-1e-5;1]*nAll;
locsAll=[1:nAll;]/nAll*domain[2];

### covariance function
function co(locs1,locs2)
  x=pairwise(Euclidean(),locs1',locs2')
  return .5*exp(-(x./.6)) + .3*exp(-(x./2).^2) + 0*(abs(x).<1e-16)
end;


### simulate data

srand(100);
predlocsInd=roundfun(nAll*3/5)+1:nAll;
locsInd=1:roundfun(nAll*3/5);
locs=locsAll[locsInd];
predvec=locsAll[predlocsInd];
if nAll<10000
    y=vec(chol(co(locsAll,locsAll))*randn(nAll,1));
  else
    y=DHsim(vec(co([locsAll[1];],locsAll)));
end
y=[1:nAll;]*1.0;
z=y[locsInd]; # +rnorm(n)*.1
n=length(z);



### partitioning
J=J*ones(Int64,M);
r=r*ones(Int64,M+1);
@time (indices,knots,data,bounds,predlocs)=partition1D(J,domain,locs,z,r,predvec);



### calculate likelihood
@time MRA(co,data,knots,indices,varEpsilon=0.2)
