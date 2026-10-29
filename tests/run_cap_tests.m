function results = run_cap_tests()
%RUN_CAP_TESTS Run the small CAP regression suite after starting INTLAB.
% From any working directory: addpath('/path/to/paper-cap/tests'); run_cap_tests

startintlab;

% Use MATLAB's native cleanup object.
cleanup_file = which('onCleanup');
if contains(cleanup_file, [filesep 'newfun' filesep 'onCleanup' filesep])
    rmpath(fileparts(cleanup_file));
end


test_dir = fileparts(mfilename('fullpath'));
project_dir = fileparts(test_dir);
previous_path = path;
restore_path = onCleanup(@()path(previous_path)); %#ok<NASGU>
addpath(genpath(fullfile(project_dir,'src')),'-begin');
assert(exist('intval','file') ~= 0, ...
    'Start INTLAB before running the CAP tests.');
results = runtests(fullfile(test_dir,'cap'));
disp(results);
assert(~isempty(results) && all([results.Passed]), ...
    'Some CAP regression tests failed or did not complete.');
end
