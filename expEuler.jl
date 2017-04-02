function expEuler(f,u0,tspan,steps = 1000)
	A = f[1]	#Matrix A of SplitODE
	g = f[2]	#function g(u,t) of SplitODE
	t = tspan[1]	#initialize t
	tmax = tspan[2]
	u = u0
	nsteps = steps
	h = tmax/nsteps		#step size
	Size = size(A)
	I = eye(Size[1])	#Identity matrix
	phi1 = ((expm(h*A)-I)/A)
	for k in 1:nsteps
  		u += phi1*(A*u + g(t,u))
  		t = k*h
	end
	return u
end
