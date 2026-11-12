function tests = test_cap_projected_targets
% Independent checks of the exact-target helper, not new non-reachability proofs.
tests = functiontests(localfunctions);
end

function testEndpointAgainstTaylorIntegral(~)
% Independently integrate exp(t*F)*f using its interval Taylor series.
% The remainder is bounded in the induced infinity norm.
T = intval(3)/50;
region = [intval(0),intval(1)/2];
level = intval(12)/25;
target = iv_coupled_constant_control_endpoint(2,T,region,level);
P = intval('pi');
for k = 1:2
    lambda = (intval(k)*P)^2;
    F = [-lambda-1,2;-2,-intval(3)/2*lambda-1];
    forcing = [intval(0);2*level*(1-cos(P*k/2))/(P*k)];
    term = T*forcing;
    integral = term;
    for j = 1:49
        term = (T*F*term)/(j+1);
        integral = integral+term;
    end
    matrix_norm = max(sum(abs(F),2));
    forcing_norm = max(abs(forcing));
    bound = T*forcing_norm*exp(T*matrix_norm);
    for j = 1:50
        bound = bound*(T*matrix_norm)/j;
    end
    bound = bound/51;
    reference = integral+midrad(zeros(2,1),sup(bound)*ones(2,1));
    difference = target(k,:).'-reference;
    assert(all(inf(difference)<=0) && all(sup(difference)>=0));
    assert(max(sup(abs(difference)))<1e-12);
end
end

function testZeroControlAndRegionAdditivity(~)
T = intval(1)/5;
region = [intval(0),intval(1)/2];
for target = {iv_coupled_constant_control_endpoint(2,intval(0),region,intval(1)), ...
        iv_coupled_constant_control_endpoint(2,T,region,intval(0))}
    value = target{1};
    assert(all(inf(value(:))<=0) && all(sup(value(:))>=0));
    assert(max(sup(abs(value(:))))<1e-14);
end
whole = iv_coupled_constant_control_endpoint(4,T,region,intval(1));
split = [intval(0),intval(1)/4;intval(1)/4,intval(1)/2];
parts = iv_coupled_constant_control_endpoint(4,T,split,intval(1));
difference = whole-parts;
assert(all(inf(difference(:))<=0) && all(sup(difference(:))>=0));
assert(max(sup(abs(difference(:))))<1e-12);
end

function testTargetDataAndWitnessConsistency(~)
% Current witness inputs are required even before the first full campaign.
root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
coefficients = readmatrix(fullfile(root,'data','witnesses','coupled_split_sine.csv'));
assert(isequal(size(coefficients),[16,2]) && all(isfinite(coefficients(:))));
assert(norm(coefficients(:))>0);
for name = {'heat_positive_coarse.mat','heat_min_time_coarse.mat','heat_intro_seed.mat'}
    saved = load(fullfile(root,'data','witnesses',name{1}),'p','Nx');
    assert(iscolumn(saved.p) && length(saved.p)==saved.Nx-1);
    assert(all(isfinite(saved.p)) && norm(saved.p)>0);
end
saved = load(fullfile(root,'data','witnesses','heat_projected_sine.mat'),'sine_coefficients');
assert(isvector(saved.sine_coefficients) && length(saved.sine_coefficients)==24);
assert(all(isfinite(saved.sine_coefficients)) && norm(saved.sine_coefficients)>0);
end
