
using SpecialMatrices, DiffEqBase
A = -full(Strang(11))
A[end,1] = 1
A[1,end] = 1
g(t,u) = 2-u
u0 = zeros(11); u0[6] = 1
nsteps = 1000;
tmax = 10.0;
h = tmax/nsteps;
u = u0;
t = 0
u = zeros(11);u[6]=1
phi1 = ((expm(h*A)-I)/A)
for k in 1:nsteps
  u += phi1*(A*u + g(t,u))
  t = k*h
end
@show u

