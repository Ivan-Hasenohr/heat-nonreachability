function control_support = iv_heat_l2_control_support(p_f,args)
% Computes the discrete control support for a spatial L2 bound M at each time.
% Processes the implicit Euler adjoint iterates in blocks to limit memory use.



Nt = args{2};
Dt = args{3};
Dx = args{10};
Nx = args{11};
control_intervals =  args{13};
if any(inf(control_intervals(2:end)) - sup(control_intervals(1:end-1)) <= 2*Dx) || any(rad(control_intervals) <= Dx/2)
    warning("The control intervals are too narrow or too close for the current boundary-cell treatment.")
end

%% Computations of S_T^star (p_f), first block

A_t = args{4}';
C = verifyinv(eye(Nx-1) - Dt*A_t);

%We use the symmetry in time, and compute the first time block, for
%i=1:2^r. Cr = C^(2^r) will be used to pass from one interval to the next

[adjoint_block, r, Cr] = create_first_block(C, p_f, Nt);
block_size = 2^r;


%% Loop 

num_blocks = ceil(Nt / block_size);

control_support = 0;

for i = 0:num_blocks-1

    start_idx = i * block_size;
    end_idx = min((i+1) * block_size - 1, Nt - 1);
    num_columns = end_idx - start_idx + 1;

    adjoint_values = adjoint_block(:, 1:num_columns);

    %Computation of the support function on this block:
    control_support = control_support + compute_block_support(adjoint_values, args);

    if i < num_blocks - 1
        adjoint_block = update_block(adjoint_block, Cr, r, i, Nt);
    end
end

end


function block_support = compute_block_support(adjoint_values, args)


Dt = args{3};
M = args{9};
Dx = args{10};
Nx = args{11};
control_intervals = args{13};

%% Computations of \sigma_{BU}(S_t_star(p_f))

%General idea:

% B is of the form union of disjoints intervals. First, we treat the case
% of mesh-intervals strictly inside of one of the B-intervals. Then, for each
% interval, we compute the value of the integral on the mesh-intervals on
% the boundary of the B-intervals.

% For each mesh-interval, the left and right values, taking into account
% the 0 endpoint values

[~,n] = size(adjoint_values);

right_value = intval(zeros(Nx,n));
right_value(1:Nx-1,:) = adjoint_values; 

left_value = intval(zeros(Nx,n));
left_value(2:Nx,:) = adjoint_values; 

coeff_1 = sqrt(2-intval('sqrt3'))/2;
coeff_2 = coeff_1*(2+intval('sqrt3'));

mix_values = (coeff_1*left_value + coeff_2*right_value).^2 + (coeff_2*left_value + coeff_1*right_value).^2; % = left^2 + left*right + right^2 but is guaranteed positive


% First case: mesh-intervals strictly in B-intervals

cell_indices = repmat((1:Nx)',1,n);

% Creates and fills the mask for mesh-intervals strictly in B-intervals
interior_cell_mask = zeros(size(left_value));

for i=1:length(control_intervals)
    control_intervals(i) = max(0, min(1,control_intervals(i)));
    j = max(1, floor(inf(intval(inf(control_intervals(i)))*intval(Nx))) + 1);
    k = min(Nx, ceil(sup(intval(sup(control_intervals(i)))*intval(Nx))));
    interior_cell_mask = interior_cell_mask | ((j < cell_indices) & (cell_indices < k));
end

% Computation of the integral for the first case


interior_contribution = Dx/3*sum(interior_cell_mask.*mix_values);

% Second case: B-interval boundaries: 
%   - for both inf and sup boundaries, we calculate in which 
%       mesh-invertal it is, and calculate accordingly


boundary_contribution = 0;




for i=1:length(control_intervals) %for every B-interval
    % left boundary of B-interval
    interval_left = inf(control_intervals(i)) ;
    interval_right = sup(control_intervals(i));

    j = max(1, floor(inf(intval(interval_left)*intval(Nx))) + 1);%index of mesh-interval containing inf(B-interval)

    % Left and right values on the inf interval:
    left_cell_left_value = left_value(j,:);
    left_cell_right_value = right_value(j,:);

    k = min(Nx, ceil(sup(intval(interval_right)*intval(Nx)))); %index of mesh-interval containing sup(B-interval)

    % Left and right values on the sup interval:
    right_cell_left_value = left_value(k,:);
    right_cell_right_value = right_value(k,:);
    
    if any(j > k) 
        error("The control interval has inconsistent boundary-cell indices.")
    elseif any(j == k)
        error("Both control boundaries lie in the same mesh cell; this case is not supported.")
    end
    
    % values on the inf and sup boundaries of the B-interval
    value_at_interval_left = (left_cell_left_value*(j*Dx-interval_left) + left_cell_right_value*(interval_left-(j-1)*Dx))/Dx;
    value_at_interval_right = (right_cell_left_value*(k*Dx-interval_right) + right_cell_right_value*(interval_right-(k-1)*Dx))/Dx;
    
    %First the inf case:

    if ~(inf(j*Dx) <= interval_left && interval_left <= sup(j*Dx)) % If the B-inf is not already on the rightside

        left_boundary_contribution = (j*Dx -interval_left)/3*((coeff_1*value_at_interval_left + coeff_2*left_cell_right_value).^2 + (coeff_2*value_at_interval_left + coeff_1*left_cell_right_value).^2);    %Compute the inegral
    else
        % Enclose the boundary contribution using a rectangle bound.
        left_boundary_contribution = infsup(zeros(1,n), sup(max(intval(0),j*Dx-interval_left).*max(left_cell_left_value.^2,left_cell_right_value.^2)));
    end

    %Second the sup case:

    if ~(inf((k-1)*Dx) <= interval_right && interval_right <= sup((k-1)*Dx))% If the B-sup is not already on the leftside

        right_boundary_contribution = (interval_right-(k-1)*Dx)/3*((coeff_1*right_cell_left_value + coeff_2*value_at_interval_right).^2 + (coeff_2*right_cell_left_value + coeff_1*value_at_interval_right).^2); %Compute the inegral
    else
        % Enclose the boundary contribution using a rectangle bound.
        right_boundary_contribution = infsup(zeros(1,n), sup(max(intval(0),interval_right-(k-1)*Dx).*max(right_cell_left_value.^2,right_cell_right_value.^2)));
    end
    
    boundary_contribution = boundary_contribution + right_boundary_contribution + left_boundary_contribution; %Sum everything


 end
% Take the spatial L2 norm at each time step, then sum over time.
block_support = M*Dt*sum(sqrt(interior_contribution + boundary_contribution),"all");


end






function AX_updated = update_block(AX_block, block_propagator, block_exponent, i, n)
    block_size = 2^block_exponent;

    next_start = (i+1) * block_size;
    next_end = min(next_start + block_size - 1, n - 1);
    num_valid_columns = next_end - next_start + 1;

    % Update the columns needed for the next block.
    AX_updated = block_propagator * AX_block(:, 1:num_valid_columns);
end

function [AX_block, block_exponent, block_propagator] = create_first_block(A, x, num_steps)
    d = size(A, 1);
    bytes_per_entry = 16;

    info = memory;
    free_bytes = floor(0.5 * info.MemAvailableAllArrays);

    max_cols = floor(free_bytes / (d * bytes_per_entry));

    % Choose the block size subject to the memory limit and the size cap.
    block_exponent =  min(13,min(floor(log2(max_cols)), floor(log2(num_steps)) + 1));
    block_size = 2^block_exponent;

    % Initialise all columns with A*x.
    AX_block = intval(repmat(A*x, 1, block_size));

    % Generate consecutive powers using binary column offsets.
    Ak = intval(A);
    for k = 1:block_exponent
        mask = bitget(0:block_size-1, k); % Add the binary powers corresponding to column offsets.
        idx = find(mask);
        AX_block(:, idx) = Ak * AX_block(:, idx);
        Ak = Ak^2;
    end

    block_propagator = Ak;  % A^(2^r), used to advance to the next block.
end