% Check comparison and norm conditions for the scalar sine targets.
% Initialize INTLAB before running. Outputs: results/heat_checks/<case>.*.
addpath(fileparts(mfilename('fullpath')));
case_names = {'heat_intro','heat_positive','heat_min_time','heat_projected'};
checks = cell(size(case_names));
for j = 1:length(case_names)
    checks{j} = verify_heat_sine_filters(case_names{j});
end
