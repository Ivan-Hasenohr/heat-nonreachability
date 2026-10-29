function base = cap_quadratic(N)
% Exact Hermite representation of q(x)=x(1-x), with h-scaled derivatives.
x = intval((0:N)')/N;
base = [x(2:N).*(1-x(2:N));(1-2*x)/N];
end
