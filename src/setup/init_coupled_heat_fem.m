function [args, args_float] = init_coupled_heat_fem(Nx,Nt,intval_T,build_float_table)
%INIT_COUPLED_HEAT_FEM Build coupled heat FEM data for certification and search.
% Components are concatenated: first component, then second component.
% Diffusion, reaction and control parameters are set below.
% Controls act on the second component only.
% intval_T must enclose the intended positive final time.
% build_float_table defaults to false; enable it for witness searches.
% Certification does not require the floating propagation table.

if nargin < 4
    build_float_table = false;
end

validateattributes(Nx, {'double'}, ...
    {'real','finite','scalar','integer','>=',3});
validateattributes(Nt, {'double'}, ...
    {'real','finite','scalar','integer','>=',2});
validateattributes(build_float_table, {'logical'}, {'scalar'});

assert(isa(intval_T,'intval') && isequal(size(intval_T),[1,1]) ...
    && isreal(intval_T) && isfinite(inf(intval_T)) ...
    && isfinite(sup(intval_T)) && inf(intval_T) > 0, ...
    'The final time must be a finite, positive scalar intval.');

%% Parameters 

%Space discretisation parameters
Dx = 1/Nx;
intval_Dx = 1/intval(Nx);

%Time discretisation parameters
T = mid(intval_T); % Floating final time.
Dt = T/Nt; % Floating time step.


% Operator parameters


%Diffusion parameters
intval_kappa_1 = midrad(1,0);
intval_kappa_2 = midrad(1.5,0);
intval_kappa = [intval_kappa_1; intval_kappa_2];
kappa = mid(intval_kappa);

% Reaction matrix entries.
a = intval('-1');
b = intval('2');
c = intval('-2');
d = intval('-1');

lambda_max_S_1 = (a+d+sqrt((a-d)^2 + (b+c)^2))/intval('2');
lambda_max_S_2 = a*intval_kappa_1 + d*intval_kappa_2 + sqrt((a*intval_kappa_1-d*intval_kappa_2)^2 + (b*intval_kappa_1 + c*intval_kappa_2)^2);
lambda_max_S_3 = a*intval_kappa_1 + d*intval_kappa_2 + sqrt((a*intval_kappa_1-d*intval_kappa_2)^2 + (c*intval_kappa_1 + b*intval_kappa_2)^2);

% Proceed only when the parameter conditions are certified with strict margins.
if ~(lambda_max_S_1 < intval('pi')^2 * min(intval_kappa_1, intval_kappa_2) && max(lambda_max_S_2, lambda_max_S_3) < intval('pi')^2  * min(intval_kappa_1^2, intval_kappa_2^2))
    error("The parameter conditions are not certified with strict margins.")
end

intval_reaction_matrix = [a b; c d];
reaction_matrix = mid(intval_reaction_matrix);

% Control domain.
control_left = 0;
control_right = 0.5;
control_intervals = infsup(control_left,control_right);

% Pointwise lower and upper control bounds.
intval_control_bounds = [intval('-1') intval('2')];
control_bounds = mid(intval_control_bounds);

% Bound on the absolute control amplitude.
intval_M = max(abs(intval_control_bounds));
M = mid(intval_M);

%% Computations

%Computation of operator A
intval_int_range = intval(linspace(1,Nx-1,Nx-1));
intval_inv_D_diag =(3/2)./( 1 + 0.5*cos(intval('pi')/intval(Nx)*intval_int_range));

intval_V = intval('sqrt2')*sin(intval('pi')/intval(Nx)*intval_int_range.*transpose(intval_int_range));

intval_M_inv = intval_V*diag(intval_inv_D_diag)*intval_V;

toeplitz_K =  intval(zeros(Nx-1,1));
toeplitz_K(1) = -2*Nx;
toeplitz_K(2) = 1*Nx;

intval_K = toeplitz(toeplitz_K);


intval_A = intval(zeros(2*Nx-2,2*Nx-2));

intval_A(1:Nx-1,1:Nx-1) = intval_kappa_1*intval_M_inv*intval_K;
intval_A(Nx:2*Nx-2,Nx:2*Nx-2) =  intval_kappa_2*intval_M_inv*intval_K;


intval_A(1:Nx-1,1:Nx-1) = intval_A(1:Nx-1,1:Nx-1) + a*eye(Nx-1);
intval_A(1:Nx-1,Nx:2*Nx-2) = b*eye(Nx-1);
intval_A(Nx:2*Nx-2,1:Nx-1) = c*eye(Nx-1);
intval_A(Nx:2*Nx-2,Nx:2*Nx-2) = intval_A(Nx:2*Nx-2,Nx:2*Nx-2) + d*eye(Nx-1);

fl_A = mid(intval_A);


% Coercivity and continuity constants for the form associated with -A.

a_0 = min(intval_kappa_1,intval_kappa_2) ...
    - max(lambda_max_S_1/intval('pi')^2,0);

sigma_max = sqrt((a^2+b^2+c^2+d^2+sqrt((a^2+b^2-c^2-d^2)^2 + 4*(a*c+b*d)^2))/2);
a_1 = max(intval_kappa_1, intval_kappa_2) + sigma_max/intval('pi')^2;

% Coercivity lower bound, followed by continuity upper bound.
form_bounds = [a_0 a_1];


% Nodal control mask for the floating search.
% Verified control-support evaluations use the intervals in args{15}.

control_mask = zeros(1,Nx+1);

for i = 1:length(control_intervals)
    interval_left = inf(control_intervals(i));
    interval_right = sup(control_intervals(i));
    control_mask = control_mask | ...
        (abs(linspace(0,1,Nx+1) - (interval_left+interval_right)/2) < (interval_right-interval_left)/2);
end

B = zeros(2*Nx-2,Nx-1);

% Inject controls into the second state component only.
B(Nx:2*Nx-2,1:Nx-1) = diag(double(control_mask(2:Nx)));

% Legacy propagation-table fields.
% The floating table stores implicit-Euler propagation data, not exponentials.
% No verified propagation table is constructed here.
intval_euler_powers = NaN;
euler_powers = NaN;

if build_float_table
    euler_powers = fl_euler_power_table(Nt, Dt, fl_A);
end

norm_L_T = -1; % Sentinel: the control-map norm is not computed here.


    
%% Initialisation of args and args_float

% Data for verified computations.
args = make_control_args(intval_T,Nt,intval_T/intval(Nt), ...
    intval_A,intval(B),form_bounds, norm_L_T, intval_euler_powers,intval_M, ...
    intval_Dx,Nx,intval_kappa,intval_reaction_matrix,intval_control_bounds,control_intervals); 

% Data for floating witness searches.
args_float = make_control_args(T,Nt,Dt,fl_A,B,mid(form_bounds), ...
    norm_L_T, euler_powers,M,Dx,Nx,kappa,mid(reaction_matrix),control_bounds,control_intervals); 

end