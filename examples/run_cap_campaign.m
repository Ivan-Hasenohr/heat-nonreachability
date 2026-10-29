function campaign_dir = run_cap_campaign(include_table)
% Run the paper certifications from a dated copy of the project.
% Start INTLAB before calling this function.
% Use run_cap_campaign(false) to skip the L2 table.

if nargin < 1
    include_table = true;
end
validateattributes(include_table,{'logical'},{'scalar'});
assert(exist('intval','file') ~= 0,'Initialize INTLAB before running the campaign.');

%% 1. Choose the examples

project_root = fileparts(fileparts(mfilename('fullpath')));

% Each row contains: case name, script in examples/, output in results/.
cases = {
    'theorem6',      'heat_positive/run_heat_intro_certificate.m',    'theorem6.mat';
    'corollary18',   'heat_positive/run_heat_positive_certificate.m', 'corollary18.mat';
    'corollary19',   'heat_min_time/run_heat_min_time_certificate.m', 'corollary19.mat';
    'proposition20', 'heat_tent/run_heat_tent_certificate.m',         'proposition20.mat';
    'proposition22', 'coupled_ball/run_coupled_ball_certificate.m',   'proposition22.mat';
    'proposition24', 'coupled_partial/run_coupled_partial.m',         'proposition24.mat';
    'table_l2',      'table_l2/compute_heat_l2_table.m',               'table_l2/table_L2_all.mat'
};

if ~include_table
    cases(end,:) = [];
end
n_cases = size(cases,1);

% Check the input files before creating the campaign directory.
assert(isfolder(fullfile(project_root,'src')),'The src directory is missing.');
assert(isfolder(fullfile(project_root,'data','witnesses')), ...
    'The witness directory is missing.');
for i = 1:n_cases
    assert(isfile(fullfile(project_root,'examples',cases{i,2})), ...
        'Missing campaign script: %s',cases{i,2});
end

%% 2. Copy the sources and witnesses

campaign_dir = fullfile(project_root,'results','campaigns', ...
    datestr(now,'yyyymmdd_HHMMSS_FFF'));
assert(~isfolder(campaign_dir),'The campaign directory already exists.');
code_root = fullfile(campaign_dir,'code');

mkdir(code_root);
mkdir(fullfile(code_root,'data'));
mkdir(fullfile(code_root,'results','table_l2'));
mkdir(fullfile(campaign_dir,'logs'));

copyfile(fullfile(project_root,'src'),fullfile(code_root,'src'));
copyfile(fullfile(project_root,'examples'),fullfile(code_root,'examples'));
copyfile(fullfile(project_root,'data','witnesses'), ...
    fullfile(campaign_dir,'initial_witnesses'));
copyfile(fullfile(campaign_dir,'initial_witnesses'), ...
    fullfile(code_root,'data','witnesses'));

% Restore the caller's path and directory when this function exits.
old_path = path;
old_directory = pwd;
restore_path = onCleanup(@()path(old_path)); %#ok<NASGU>
restore_directory = onCleanup(@()cd(old_directory)); %#ok<NASGU>
addpath(genpath(fullfile(code_root,'src')),'-begin');

fprintf('Campaign directory:\n%s\n',campaign_dir);
fprintf('New result files:\n%s\n',fullfile(code_root,'results'));

%% 3. Prepare the summary

summary.started_at = datestr(now,31);
summary.matlab_version = version;
summary.intlab_constructor = which('intval');
summary.source_project = project_root;
summary.include_table = include_table;
summary.completed = false;
summary.all_certified = false;

for i = 1:n_cases
    summary.cases(i).name = cases{i,1};
    summary.cases(i).status = 'pending';
    summary.cases(i).wall_seconds = NaN;
    summary.cases(i).result_file = '';
    summary.cases(i).message = '';
end

summary_file = fullfile(campaign_dir,'summary.mat');
save(summary_file,'summary');

%% 4. Run each case and save its status

for i = 1:n_cases
    case_name = cases{i,1};
    diary(fullfile(campaign_dir,'logs',[case_name,'.txt']));
    diary on;
    fprintf('\n[%d/%d] %s -- started %s\n',i,n_cases,case_name,datestr(now,31));
    started = tic;

    try
        % A separate workspace prevents variables leaking between scripts.
        execute_script(fullfile(code_root,'examples',cases{i,2}));

        result_file = fullfile(code_root,'results',cases{i,3});
        certified = check_campaign_result(result_file,strcmp(case_name,'table_l2'));

        if certified
            summary.cases(i).status = 'certified';
        else
            summary.cases(i).status = 'not_certified';
        end
        summary.cases(i).result_file = fullfile('code','results',cases{i,3});
    catch failure
        % Keep the error in the summary and continue with the next case.
        summary.cases(i).status = 'error';
        summary.cases(i).message = failure.message;
        fprintf(2,'%s\n',getReport(failure,'extended','hyperlinks','off'));
    end

    summary.cases(i).wall_seconds = toc(started);
    fprintf('[%d/%d] %s -- %s (%.1f s)\n',i,n_cases,case_name, ...
        summary.cases(i).status,summary.cases(i).wall_seconds);
    diary off;
    save(summary_file,'summary');
end

%% 5. Save and display the final summary

summary.completed = true;
summary.finished_at = datestr(now,31);
summary.all_certified = all(strcmp({summary.cases.status},'certified'));
save(summary_file,'summary');

disp(struct2table(summary.cases));
fprintf('All selected results certified: %d\n',summary.all_certified);
fprintf('Summary saved to:\n%s\n',summary_file);
end

function execute_script(script_file)
% Run one script without sharing the launcher's workspace.
run(script_file);
end
