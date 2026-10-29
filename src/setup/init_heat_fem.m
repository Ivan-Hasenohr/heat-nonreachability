function [args, args_float] = init_heat_fem(Nx,Nt,intval_T,build_float_table)
%INIT_HEAT_FEM Build scalar heat FEM data for certification and search.
% The diffusion coefficient, control bound and control intervals are set below.
% intval_T must enclose the intended positive final time.
% Floating data use midpoint approximations of the verified parameters.
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


% Space discretization parameters.
Dx = 1/Nx; 
intval_Dx = 1/intval(Nx);

%Time discretisation parameters

T = mid(intval_T); % Floating final time.
Dt = T/Nt; % Floating time step.
intval_Dt = intval_T/intval(Nt); % Verified time step.

%Other parameters:

intval_kappa = intval('1'); %Diffusion coefficient
kappa = mid(intval_kappa);

% Outward enclosures of the exact control intervals.
control_intervals = [infsup(inf(intval(1)/5), sup(intval(2)/5)), infsup(inf(intval(4)/5), 1)];


%Verified computation of operator A
intval_int_range = intval(linspace(1,Nx-1,Nx-1));
intval_inv_D_diag =(3/2)*intval(Nx)./( 1 + 0.5*cos(intval('pi')/intval(Nx)*intval_int_range));


intval_V = sqrt(2/intval(Nx))*sin(intval('pi')/intval(Nx)*intval_int_range.*transpose(intval_int_range));


intval_M_inv = intval_V*diag(intval_inv_D_diag)*intval_V;
fl_M_inv = mid(intval_M_inv); % Floating approximation of the inverse mass matrix.


toeplitz_K =  intval(zeros(Nx-1,1));
toeplitz_K(1) = 2*Nx;
toeplitz_K(2) = -1*Nx;

intval_K = toeplitz(toeplitz_K);
fl_K = mid(intval_K);


intval_A = -intval_kappa*intval_M_inv*intval_K;
fl_A = -kappa*(fl_M_inv*fl_K);

% Nodal control mask for the floating search.
% This matrix does not describe partially controlled mesh cells.
% Verified control-support evaluations use the intervals in args{13}.

% Nodal control mask used only by the floating search.
x_nodes = (1:Nx-1)'/Nx;
control_mask = false(Nx-1,1);

for j = 1:length(control_intervals)
    control_mask = control_mask | ...
        (x_nodes >= inf(control_intervals(j)) & ...
         x_nodes <= sup(control_intervals(j)));
end

B = diag(double(control_mask));


%Other parameters:

alpha = 0; % Scalar Dirichlet diffusion is self-adjoint; a_0 = a_1.
intval_alpha = intval(0);

% Legacy propagation-table fields.
% The floating table stores implicit-Euler propagation data, not exponentials.
% No verified propagation table is constructed here.
intval_euler_powers = NaN;
euler_powers = NaN;

if build_float_table
    euler_powers = fl_euler_power_table(Nt, Dt, fl_A);
end

norm_L_T = -1; % Sentinel: the control-map norm is not computed here.



%% Additional parameters specific to the control problem

M=1; %Control bound. Must verify M > 0.

%% Initialisation of args and args_float

% Data for verified computations.
args = make_control_args(intval_T,Nt,intval_Dt,intval_A,intval(B), ...
    intval_alpha, norm_L_T, intval_euler_powers,M,intval_Dx,Nx,intval_kappa, control_intervals); 

% Data for floating witness searches.
args_float = make_control_args(T,Nt,Dt,fl_A,B,alpha, norm_L_T, ...
    euler_powers,M,Dx,Nx,kappa, control_intervals); 
end