N=32400;

for i=1:10
	x=randn(N,N);
	x=x+x';
	for i=1:N
        	x[i,i] = x[i,i] + 1000;
	end
	y=x;
	y[1,1] = y[1,1] * 1.0;
	tic()
	 ccall((:magmaCholeskyFinal,"mylib.so"),Void,(Ptr{Cdouble},Int32),x,N)
	toc();
end 

