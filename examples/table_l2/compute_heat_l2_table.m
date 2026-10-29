% Compute all numerical rows of the L2 table.
project_root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(genpath(fullfile(project_root,'src')));

Nx_list = [100, 200, 500, 1000, 2000];
Nt_list = [2000, 8000, 50000, 200000, 800000];
intval_T = intval('1');

n_rows = length(Nx_list);
n_modes = 8;

lambda_num_all = NaN(n_rows,n_modes);
lambda_cert_all = intval(NaN(n_rows,n_modes));
sups_all = NaN(n_rows,n_modes);
J_all = intval(NaN(n_rows,n_modes));
p_all = cell(n_rows,n_modes);
args_all = cell(n_rows,1);

output_dir = fullfile(project_root,'results','table_l2');
if ~isfolder(output_dir)
    mkdir(output_dir);
end

output_file = fullfile(output_dir,'table_L2_all.mat');

started = tic;

for row = 1:n_rows
    Nx = Nx_list(row);
    Nt = Nt_list(row);

    fprintf('Row %d/%d: Nx = %d, Nt = %d\n', ...
        row,n_rows,Nx,Nt);



[args,args_float] = init_heat_fem(Nx,Nt,intval_T,false);


% The L2 table uses controls over the entire spatial domain.
args{5} = intval(eye(Nx-1));
args{13} = infsup(0,1);

intval_T = args{1};
intval_Dt = args{3};
intval_B = args{5};
intval_alpha = args{6};
intval_M = args{9};
intval_Dx = args{10};
Nx = args_float{11};
intval_kappa = args{12};

x = transpose(linspace(1/Nx,1-1/Nx,Nx-1));

lambda_num = NaN(1,n_modes);
sups = NaN(1,n_modes);

for k = 1:n_modes
p = sqrt(2)*sin(k*pi*x);
y_k = (transpose(1:k) == k);
norm_y_k = 1/sqrt(intval(2));
p_base = iv_natural_cubic_from_nodes(p);

control_support = iv_heat_l2_control_support(p,args);

[base_error_bound, interpolation_error_bound] = compute_amplitude_search_bounds(p_base,intval_T,intval_Dt,intval_B,intval_alpha,intval_M,intval_Dx,intval_kappa);

% Include the two zero Dirichlet boundary values.
pairing_y_f = iv_p1_sine_series_pairing([intval(0); intval(p(:)); intval(0)], y_k, intval_Dx);

lambda = inf((control_support + base_error_bound)/(pairing_y_f-interpolation_error_bound));

err = iv_heat_discretization_bound(p_base, 0,lambda*norm_y_k,args);

i = 0;

if sup(control_support - lambda*pairing_y_f + err) >= 0
while sup(control_support - lambda*pairing_y_f + err) >= 0  && i < 1000
    i = i+1;
    lambda = lambda + 1e-5;
    err = iv_heat_discretization_bound(p_base, 0,lambda*norm_y_k,args);
end
else 
while sup(control_support - lambda*pairing_y_f + err) < 0  && i < 1000
    i = i+1;
    lambda = lambda - 1e-4;
    err = iv_heat_discretization_bound(p_base, 0,lambda*norm_y_k,args);
end
lambda = lambda + 1e-4;
end
lambda_num(k) = 1e-4*ceil(1e4*lambda);

% Certify the exact decimal amplitude reported in the table.
lambda_cert = intval(ceil(1e4*lambda))/10000;
err = iv_heat_discretization_bound(p_base, 0, lambda_cert*norm_y_k, args);

sups(k) = sup(control_support - lambda_cert*pairing_y_f + err);

if ~(isfinite(sups(k)) && sups(k) < 0)
    error("The rounded amplitude was not certified.")
end

disp([row, k, lambda_num(k)])

J_opt = control_support - lambda_cert*pairing_y_f + infsup(-err,err);

lambda_cert_all(row,k) = lambda_cert;
J_all(row,k) = J_opt;
p_all{row,k} = p;


end

    % Store the completed row and its verified model parameters.
    lambda_num_all(row,:) = lambda_num;
    sups_all(row,:) = sups;
    args_all{row} = args;

    % Save all completed rows after each spatial mesh.
    completed_rows = row;
    elapsed_seconds = toc(started);

    save(output_file, ...
        'Nx_list','Nt_list','intval_T','n_modes', ...
        'lambda_num_all','lambda_cert_all','sups_all','J_all', ...
        'p_all','args_all','elapsed_seconds','completed_rows');

    fprintf('Saved row %d/%d.\n',completed_rows,n_rows);
end

disp('Certified amplitudes:');
disp(lambda_num_all);
disp('Results saved to:');
disp(output_file);

function [base_error_bound, interpolation_error_bound] = compute_amplitude_search_bounds(p_opt_base,T,Dt,B,alpha,M,h,kappa)


norm_d2_p_f = iv_cubic_hermite_d2_norm(p_opt_base);

norm_A_star_p_f = kappa*norm_d2_p_f;
interpolation_error_bound = h^2/intval('pi')^2 * norm_d2_p_f;

norm_y_0 = 0;

norm_B = norm(B);

C_0 = 2/(intval('pi')*kappa);
C_alpha = (2+2/intval('sqrt3'));
a_0 = kappa;
a_1 = kappa;
C_1 = a_1^(3/2)*C_0^2/(4*sqrt(a_0));
C_2 = C_1*(7+4*intval('log2')/cos(alpha)+C_alpha);

err_1 = (0.5*M*T*norm_B*Dt + (norm_y_0 + M*T*norm_B)*(C_2*h^2+C_alpha/(2*exp(intval(1))*cos(alpha))*Dt))*norm_A_star_p_f;

base_error_bound = err_1 + ((norm_y_0 + M*T*norm_B)*C_alpha)*interpolation_error_bound;
end