function checks = iv_heat_sine_filter_checks(T,coefficients,subdivisions)
%IV_HEAT_SINE_FILTER_CHECKS Comparison and norm obstructions for scalar heat.
% Fixed model: y_t=y_xx+chi_omega*u, zero initial/boundary data on (0,1),
% omega=(1/5,2/5) union (4/5,1), and 0<=u<=1.
% The point target is sum_k coefficients(k)*sin(k*pi*x).
% INTLAB must already be initialized. No FEM mesh or witness is needed.
% Status: 'satisfied', 'nonreachable', or 'inconclusive'. A satisfied
% necessary condition does not establish reachability.

if nargin<3
    subdivisions = 400;
end
validateattributes(subdivisions,{'double'}, ...
    {'scalar','real','finite','integer','>=',10});
assert(isa(T,'intval') && isequal(size(T),[1,1]) && isreal(T) ...
    && isfinite(sup(T)) && inf(T)>0, 'T must be a positive scalar intval.');
assert(isa(coefficients,'intval') && isreal(coefficients) ...
    && size(coefficients,2)==1 && size(coefficients,1)>=1 ...
    && all(isfinite(inf(coefficients))) && all(isfinite(sup(coefficients))), ...
    'Provide finite interval sine coefficients as a column.');

P = intval('pi');
N = 12;
k = intval((1:N)');

% w solves -w''=chi_omega, with zero boundary values.
% The maximal reachable trajectory is v(T)=w-exp(T*Delta)*w.
stationary_coefficients = 2*(cos(P*k/5)-cos(2*P*k/5) ...
    +cos(4*P*k/5)-cos(P*k))./(P*k).^3;
transient_coefficients = stationary_coefficients.*exp(-(P*k).^2*T);

% |w_k|<=8/(pi*k)^3 and |sin(k*pi*x)|<=k*sin(pi*x).
% Thus the omitted transient, divided by sin(pi*x), is bounded by
% 8*exp(-pi^2*(N+1)^2*T)/(pi^3*N), using sum_{k>N}1/k^2<=1/N.
tail = 8*exp(-P^2*intval(N+1)^2*T)/(P^3*intval(N));

%% Cover the whole domain with interval cells, split at control endpoints
knots = [intval(0),intval(1)/5,intval(2)/5,intval(4)/5,intval(1)];
left = intval(zeros(4*subdivisions,1));
right = left;
steady_ratio = left;
for part = 1:4
    rows = (part-1)*subdivisions+(1:subdivisions);
    a = knots(part)+(knots(part+1)-knots(part))*intval((0:subdivisions-1)')/subdivisions;
    b = knots(part)+(knots(part+1)-knots(part))*intval((1:subdivisions)')/subdivisions;
    X = infsup(max(0,inf(a)),min(1,sup(b)));
    left(rows) = a;
    right(rows) = b;
    switch part
        case 1
            steady_ratio(rows) = intval(4)/25*distance_over_sine(X,P);
        case 2
            steady_ratio(rows) = (intval(4)/25*X-(X-intval(1)/5).^2/2)./sin(P*X);
        case 3
            steady_ratio(rows) = (intval(3)/50-X/25)./sin(P*X);
        case 4
            V = infsup(max(0,inf(1-X)),min(1,sup(1-X)));
            steady_ratio(rows) = (X/2-intval(13)/50).*distance_over_sine(V,P);
    end
end
X = infsup(max(0,inf(left)),min(1,sup(right)));
Z = cos(P*X);

% sin(k*pi*x)/sin(pi*x)=U_(k-1)(cos(pi*x)), including endpoint limits.
target_ratio = sine_ratio(coefficients,Z);
transient_ratio = sine_ratio(transient_coefficients,Z);
lower_gap = target_ratio;
upper_gap = steady_ratio-transient_ratio-target_ratio ...
    +midrad(zeros(size(X)),sup(tail)*ones(size(X)));

checks.lower_comparison = summarize_gap(lower_gap);
checks.upper_comparison = summarize_gap(upper_gap);
checks.lower_comparison.violation_point = violation_point(lower_gap,X);
checks.upper_comparison.violation_point = violation_point(upper_gap,X);
checks.transient_tail_over_sine = sup(tail);
checks.spatial_cells = 4*subdivisions;

%% Necessary L2 and H1-seminorm bounds for every admissible trajectory
modes = intval((1:size(coefficients,1))');
target_norm = sqrt(sum(coefficients.^2)/2);
target_gradient_norm = sqrt(sum((P*modes.*coefficients).^2)/2);
forcing_bound = sqrt(intval(2)/5);

% L2 contraction: ||S(t)||<=exp(-pi^2*t).
l2_bound = forcing_bound*(1-exp(-P^2*T))/P^2;
% Cauchy-Schwarz in each sine mode: ||y_x(T)||<=M*sqrt(T/2).
% Integrating ||(-Delta)^(1/2)S(t)|| also gives 2*M/(pi*sqrt(e)).
h1_bound = forcing_bound*min(sqrt(T/2),2/(P*exp(intval(1)/2)));
checks.l2_norm = summarize_gap(l2_bound-target_norm);
checks.h1_seminorm = summarize_gap(h1_bound-target_gradient_norm);
checks.target_l2_norm = [inf(target_norm),sup(target_norm)];
checks.target_h1_seminorm = [inf(target_gradient_norm),sup(target_gradient_norm)];
checks.l2_bound = [inf(l2_bound),sup(l2_bound)];
checks.h1_seminorm_bound = [inf(h1_bound),sup(h1_bound)];

statuses = {checks.lower_comparison.status,checks.upper_comparison.status, ...
    checks.l2_norm.status,checks.h1_seminorm.status};
checks.all_conditions_satisfied = all(strcmp(statuses,'satisfied'));
checks.certified_nonreachability = any(strcmp(statuses,'nonreachable'));
checks.T = [inf(T),sup(T)];
checks.target_coefficients = [inf(coefficients),sup(coefficients)];
end

function value = distance_over_sine(X,P)
% x/sin(pi*x) is increasing on [0,1/5]; its limit at zero is 1/pi.
a = inf(X); b = sup(X);
lower = intval(ones(size(a)))/P;
upper = lower;
positive = a>0;
lower(positive) = intval(a(positive))./sin(P*intval(a(positive)));
positive = b>0;
upper(positive) = intval(b(positive))./sin(P*intval(b(positive)));
value = infsup(inf(lower),sup(upper));
end

function value = sine_ratio(coefficients,Z)
previous = intval(zeros(size(Z)));
current = intval(ones(size(Z)));
value = intval(zeros(size(Z)));
for k = 1:size(coefficients,1)
    value = value+coefficients(k)*current;
    next = 2*Z.*current-previous;
    previous = current;
    current = next;
end
end

function result = summarize_gap(gap)
% Nonnegative gap: necessary condition satisfied. Negative upper bound:
% an actual violation is certified. Overlap with zero is inconclusive.
result.minimum_gap = [min(inf(gap)),min(sup(gap))];
if all(inf(gap)>=0)
    result.status = 'satisfied';
elseif any(sup(gap)<0)
    result.status = 'nonreachable';
else
    result.status = 'inconclusive';
end
end

function point = violation_point(gap,X)
index = find(sup(gap)<0,1);
point = [];
if ~isempty(index)
    point = mid(X(index));
end
end
