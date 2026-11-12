function test_heat_sine_filters
% Standalone regression checks; does not modify or extend run_cap_tests.
% Initialize INTLAB before calling this function.
root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root,'examples','heat_checks'));

% Nonnegative forcing gives y>=0: a strictly negative target is excluded.
negative = iv_heat_sine_filter_checks(intval(1),-intval(1)/100,100);
assert(strcmp(negative.lower_comparison.status,'nonreachable'));
assert(negative.certified_nonreachability);
assert(negative.lower_comparison.violation_point>0 ...
    && negative.lower_comparison.violation_point<1);

% A large positive target violates comparison and both norm bounds.
large = iv_heat_sine_filter_checks(intval(1),intval(1),100);
assert(strcmp(large.upper_comparison.status,'nonreachable'));
assert(strcmp(large.l2_norm.status,'nonreachable'));
assert(strcmp(large.h1_seminorm.status,'nonreachable'));
assert(large.certified_nonreachability);

% The zero target is genuinely reachable, with equality in lower comparison.
zero = iv_heat_sine_filter_checks(intval(1),intval(0),100);
assert(zero.all_conditions_satisfied && ~zero.certified_nonreachability);
assert(all(zero.lower_comparison.minimum_gap==0));

% Parameter uncertainty spanning positive and negative targets is inconclusive.
uncertain = iv_heat_sine_filter_checks(intval(1),midrad(0,.001),100);
assert(strcmp(uncertain.lower_comparison.status,'inconclusive'));
assert(~uncertain.certified_nonreachability);

% For a pure sine, the exact norm and derivative norm are known independently.
P = intval('pi');
small = iv_heat_sine_filter_checks(intval(1),intval(1)/50,100);
exact_l2 = intval(1)/(50*sqrt(intval(2)));
exact_h1 = P/(50*sqrt(intval(2)));
assert(small.target_l2_norm(1)<=sup(exact_l2) ...
    && small.target_l2_norm(2)>=inf(exact_l2));
assert(small.target_h1_seminorm(1)<=sup(exact_h1) ...
    && small.target_h1_seminorm(2)>=inf(exact_h1));

for name = {'heat_intro','heat_positive','heat_min_time','heat_projected'}
    result = verify_heat_sine_filters(name{1});
    assert(result.all_conditions_satisfied && ~result.certified_nonreachability);
end
disp('All scalar comparison/norm checks passed, including rejection and inconclusive cases.');
end
