% Compare the raw and regularized witnesses for Proposition 20.
% Both evaluations use the original, unregularized separation functional.
% The raw fine-mesh evaluation can take several minutes.

project_root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(genpath(fullfile(project_root,'src')));

% To use a campaign directly, set result_dir to its code/results directory.
result_dir = fullfile(project_root,'results');
source_result_file = fullfile(result_dir,'proposition20.mat');
source_raw_witness = fullfile(project_root,'data','witnesses', ...
    'proposition20_unregularized.mat');
assert(isfile(source_result_file), ...
    'Run run_heat_tent_certificate or extract the reference results archive first.');

reg = load(source_result_file);
raw = load(source_raw_witness,'p');
Nx = reg.args{11};

% Transfer and normalize exactly as in the certification scripts.
p_f = resample_nodal_witness(raw.p,Nx);
p_f = p_f/(norm(p_f)*sqrt(1/Nx));
p_f_reg = reg.p;
norm_d2_pf = iv_cubic_hermite_d2_norm(iv_natural_cubic_from_nodes(p_f));
norm_d2_pf_reg = iv_cubic_hermite_d2_norm(iv_natural_cubic_from_nodes(p_f_reg));

[certified_pf,J_pf,~,~,~,bounds_pf] = certify_heat_nonreachability( ...
    p_f,reg.y_0,reg.y_f_base,reg.args);
J_pf_reg = reg.J;
bounds_pf_reg = reg.bounds;
J_discrete_pf = mid(bounds_pf.discrete_enclosure);
J_discrete_pf_reg = mid(bounds_pf_reg.discrete_enclosure);

output_file = fullfile(result_dir,'tent_witness_comparison.mat');
save(output_file,'p_f','p_f_reg','Nx','certified_pf','J_pf','J_pf_reg', ...
    'bounds_pf','bounds_pf_reg','norm_d2_pf','norm_d2_pf_reg', ...
    'J_discrete_pf','J_discrete_pf_reg','source_result_file','source_raw_witness');

format infsup
norm_d2_pf
norm_d2_pf_reg
J_pf
J_pf_reg
J_discrete_pf
J_discrete_pf_reg
fprintf('Comparison saved to:\n%s\n',output_file);
