function target = iv_coupled_constant_control_endpoint(N,T,regions,level)
%IV_COUPLED_CONSTANT_CONTROL_ENDPOINT Enclose the first N endpoint sine modes.
% Zero initial state, Dirichlet boundary conditions on (0,1), and equations
%   y1_t = y1_xx - y1 + 2*y2,
%   y2_t = (3/2)*y2_xx - 2*y1 - y2 + chi_regions*level.
% The control is constant in time and on each of the disjoint regions.
% Row k contains the coefficients of sin(k*pi*x) in both components.
% The returned coefficients DEFINE projected targets; FEM/Euler propagation
% and its error estimates are still used to certify non-reachability.

validateattributes(N,{'double'},{'scalar','real','finite','integer','positive'});
assert(isa(T,'intval') && isequal(size(T),[1,1]) && isreal(T) ...
    && isfinite(sup(T)) && inf(T)>=0, 'T must be a nonnegative scalar intval.');
assert(isa(level,'intval') && isequal(size(level),[1,1]) && isreal(level) ...
    && isfinite(inf(level)) && isfinite(sup(level)), ...
    'level must be a finite scalar intval.');
assert(isa(regions,'intval') && isreal(regions) && size(regions,2)==2 ...
    && size(regions,1)>=1 && all(inf(regions(:,1))>=0) ...
    && all(sup(regions(:,2))<=1) ...
    && all(sup(regions(:,1))<inf(regions(:,2))), ...
    'regions must have one ordered pair of interval endpoints per row.');
if size(regions,1)>1
    assert(all(sup(regions(1:end-1,2))<=inf(regions(2:end,1))), ...
        'Control regions must be ordered and disjoint.');
end

P = intval('pi');
k = intval((1:N)');
lambda = (P*k).^2;
alpha = lambda+1;
beta = intval(3)/2*lambda+1;

% On frequency k, the generator is F = [-alpha,2;-2,-beta].
% Its eigenvalues are tau +/- delta, with delta real for every k >= 1.
delta = sqrt((beta-alpha).^2/4-4);
assert(all(inf(delta)>0));
tau = -(alpha+beta)/2;
difference = (beta-alpha)/2;
slow = exp((tau+delta)*T);
fast = exp((tau-delta)*T);
e12 = (slow-fast)./delta;
e22 = (slow+fast)/2-(difference./(2*delta)).*(slow-fast);

% Sine coefficients of the forcing, including the factor 2 from ||sin||^2.
moments = intval(zeros(N,1));
for j = 1:size(regions,1)
    moments = moments+(cos(P*k*regions(j,1))-cos(P*k*regions(j,2)))./(P*k);
end
forcing = 2*level*moments;
determinant = alpha.*beta+4;

% Integral_0^T exp(t*F)*[0;forcing] dt = F^(-1)*(exp(T*F)-Id)*[0;forcing].
target = [forcing.*(-beta.*e12-2*(e22-1))./determinant, ...
    forcing.*(2*e12-alpha.*(e22-1))./determinant];
end
