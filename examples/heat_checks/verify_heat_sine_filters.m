function checks = verify_heat_sine_filters(case_name,subdivisions)
%VERIFY_HEAT_SINE_FILTERS Check the scalar sine targets used in the article.
% INTLAB must be initialized. This function does not run the FEM certifier.
% Inputs are stated explicitly below, independently of the campaign scripts.
if nargin<1
    case_name = 'heat_positive';
end
if nargin<2
    subdivisions = 400;
end
switch case_name
    case 'heat_intro'
        T = intval('1.15');
        coefficients = intval('0.037');
    case 'heat_positive'
        T = intval(1);
        coefficients = intval('0.035');
    case 'heat_min_time'
        T = intval(13)/100;
        coefficients = intval(1)/50;
    case 'heat_projected'
        T = intval(13)/100;
        k = intval((1:6)'); P = intval('pi');
        coefficients = 2*(intval(3)/4)*(cos(P*k/5)-cos(2*P*k/5) ...
            +cos(4*P*k/5)-cos(P*k)).*(1-exp(-(P*k).^2*T))./(P*k).^3;
        coefficients(5) = intval(0);
    otherwise
        error('Choose heat_intro, heat_positive, heat_min_time or heat_projected.');
end

checks = iv_heat_sine_filter_checks(T,coefficients,subdivisions);
checks.case_name = case_name;
project_root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
output_dir = fullfile(project_root,'results','heat_checks');
if ~isfolder(output_dir)
    mkdir(output_dir);
end
save(fullfile(output_dir,[case_name '.mat']),'checks');
fid = fopen(fullfile(output_dir,[case_name '.json']),'w');
assert(fid>=0);
fprintf(fid,'%s\n',jsonencode(checks));
fclose(fid);

fprintf('\n%s: comparison lower / upper = %s / %s.\n',case_name, ...
    checks.lower_comparison.status,checks.upper_comparison.status);
fprintf('L2 norm / H1 seminorm = %s / %s.\n', ...
    checks.l2_norm.status,checks.h1_seminorm.status);
fprintf('Certified non-reachability by these obstructions: %d.\n', ...
    checks.certified_nonreachability);
end
