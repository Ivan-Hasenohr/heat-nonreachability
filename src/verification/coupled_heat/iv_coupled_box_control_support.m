function control_support = iv_coupled_box_control_support(p_f,args)
% Enclose the discrete control support sum for a <= u <= b.
% The control acts on the second component of the coupled heat equation.
% Use piecewise affine spatial functions with zero boundary values.
% Time indices are q = 1,...,Nt.
% Continuous-model error bounds are handled elsewhere.

Nt = args{2};
Dt = args{3};
Dx = args{10};
Nx = args{11};
control_intervals =  args{15};

if any(inf(control_intervals(2:end)) - sup(control_intervals(1:end-1)) <= 2*Dx) || any(rad(control_intervals) <= Dx/2)
    warning("The control intervals may be too narrow or too close for this mesh.")
end

%% Computations of S_T^star (p_f), first block

A_star = args{4}';
C = verifyinv(eye(2*Nx-2) - Dt*A_star);

% Build the first block [C*p_f, ..., C^block_size*p_f].
% block_propagator encloses C^block_size and advances each block to the next one.

[adjoint_block, r, block_propagator] = create_first_block(C, p_f, Nt);
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
        adjoint_block = update_block(adjoint_block, block_propagator, r, i, Nt);
    end
end

end


function block_support = compute_block_support(adjoint_values, args)
% Computes the support function for constraints of the type a <= u <= b
% Only control bounds satisfying a <= 0 <= b are supported.


Dt = args{3};
h = args{10};
Nx = args{11};
control_intervals = args{15};
control_bounds = args{14};
a = min(control_bounds); b = max(control_bounds);


if ~(sup(a) <= 0 && inf(b) >= 0)
    error("This implementation requires a <= 0 <= b.")
end
%% Computations of \sigma_{BU}(S_t_star(p_f))


% For each mesh-interval, the left and right values, taking into account
% the 0 endpoint values

[~,num_columns] = size(adjoint_values);

right_value = intval(zeros(Nx,num_columns));
right_value(1:Nx-1,:) = adjoint_values(Nx:2*Nx-2,:); 

left_value = intval(zeros(Nx,num_columns));
left_value(2:Nx,:) = adjoint_values(Nx:2*Nx-2,:); 

% First case: mesh-intervals strictly in B-intervals

% Interior cells of the control region.
interior_cell_mask = zeros(size(left_value));

for i=1:length(control_intervals)
    control_intervals(i) = max(0, min(1,control_intervals(i)));
    
    j = max(1, floor(inf(intval(inf(control_intervals(i)))*intval(Nx))) + 1);
    k = min(Nx, ceil(sup(intval(sup(control_intervals(i)))*intval(Nx))));
    cell_indices = repmat((1:Nx)', 1, num_columns);
    interior_cell_mask = interior_cell_mask | ((j < cell_indices) & (cell_indices < k));
end


interior_contribution = Dt*sum(sum(interior_cell_mask.*integrate_affine_box_support(left_value, right_value, h, a, b)), "omitnan");

%% All cases inside each interval:

boundary_contribution = 0;

for i=1:length(control_intervals) 
    
    
    interval_left = inf(control_intervals(i)) ;
    interval_right = sup(control_intervals(i));

    j = max(1, floor(inf(intval(interval_left)*intval(Nx))) + 1); % Left boundary cell

    % Left boundary cell of the control interval.
    left_cell_left_value = left_value(j,:);
    left_cell_right_value = right_value(j,:);

    k = min(Nx, ceil(sup(intval(interval_right)*intval(Nx)))); % Right boundary cell

    % Right boundary cell of the control interval.
    right_cell_left_value = left_value(k,:);
    right_cell_right_value = right_value(k,:);
    
    if any(j > k) 
        error("The left boundary cell index exceeds the right boundary cell index.")
    elseif any(j == k)
        error("Control intervals contained in a single mesh cell are not supported.")
    end
    
    
    %% First the B-inf case:

    % Interpolate the field at the left control boundary.
    value_at_interval_left = (left_cell_left_value*(j*h-interval_left) + left_cell_right_value*(interval_left-(j-1)*h))/h;

    if ~(inf(j*h) <= interval_left && interval_left <= sup(j*h)) % If the B-inf is not already on the rightside        
        left_boundary_contribution = Dt*sum(integrate_affine_box_support(value_at_interval_left, left_cell_right_value, j*h-interval_left, a, b), "all");

    else
        % Bound the left boundary contribution by a rectangle.
        left_boundary_contribution = infsup(0, sup(Dt*max(abs(a),abs(b))*sum(max(intval(0), j*h-interval_left) .* max(abs(left_cell_left_value), abs(left_cell_right_value)), "all")));
    end



    %% Second the B-sup case:

   % Interpolate the field at the right control boundary.
    value_at_interval_right = (right_cell_left_value*(k*h-interval_right) + right_cell_right_value*(interval_right-(k-1)*h))/h;

    if ~(inf((k-1)*h) <= interval_right && interval_right <= sup((k-1)*h))% If the B-sup is not already on the leftside
         right_boundary_contribution = Dt*sum(integrate_affine_box_support(right_cell_left_value, value_at_interval_right, interval_right-(k-1)*h, a, b), "all");
    else
        % Bound the right boundary contribution by a rectangle.
        right_boundary_contribution = infsup(0, sup(Dt*max(abs(a),abs(b))*sum(max(intval(0), interval_right-(k-1)*h) .* max(abs(right_cell_left_value), abs(right_cell_right_value)), "all")));
    end

    
    %% Final sum for the B-interval

    boundary_contribution = boundary_contribution + right_boundary_contribution + left_boundary_contribution; %Sum everything
    
    
 end

block_support = interior_contribution + boundary_contribution;
end


function integrals = integrate_affine_box_support(left_value, right_value, h, a, b)
% Enclose the integral of b*max(f,0) + a*min(f,0) over [0,h].
% The affine function f has endpoint values left_value and right_value.
% Requires a <= 0 <= b.
% Intermediate contributions are scaled by h/2 at the end.

if ~(sup(a) <= 0 && inf(b) >= 0)
    error("This implementation requires a <= 0 <= b.")
end

% Requires many masks to treat the cases where the values are positive,
% negative or mixed
integrals= intval(zeros(size(left_value)));

max_value = max(left_value, right_value); % max of left and right value of mesh-intervals 
min_value = min(left_value, right_value); % min of left and right value of mesh-intervals

endpoint_sum = left_value + right_value;


% Mask creation for mesh-intervals
mask_zero = (sup(min_value) >= 0) & (inf(max_value) <= 0); % if 0 is in both min and max values 
mask_pos = (1-mask_zero) & (min_value >= 0); % if not and if both are positive or zero
mask_neg = (1-mask_zero) & (max_value <= 0); % if not and if both are negative or zero
mask_mix = (1-mask_neg) & (1-mask_pos) & (1-mask_zero); % the rest (that is sup(max_value) > 0 > inf(min_value)


%% Zero, positive or negative case: [a*(signed area of inf trapeze), b*(signed area of sup trapeze)] inside B

integrals(mask_pos) = b*(endpoint_sum(mask_pos));
integrals(mask_neg) = a*(endpoint_sum(mask_neg));
integrals(mask_zero) = hull(intval(0), hull(a*intval(inf(endpoint_sum(mask_zero))), b*intval(sup(endpoint_sum(mask_zero)))));

%% Mixed case

% Compute an auxiliary value from endpoint values closest to zero.
% Zero is included below to ensure a conservative lower bound.

% 0 if 0 in left_value, inf(left_value) if left_value >= 0, sup(left_value)
% if left_value <= 0:
left_nearest_zero = sign(left_value).*intval(inf(abs(left_value))); 

%Same for right_value:
right_nearest_zero = sign(right_value).*intval(inf(abs(right_value)));


min_nearest_zero = min(left_nearest_zero,right_nearest_zero);
max_nearest_zero = max(left_nearest_zero,right_nearest_zero);

% Auxiliary expression using endpoint values nearest zero.
% The final enclosure includes zero to retain a conservative lower bound.

mixed_auxiliary_value = b*max_nearest_zero.^2./(max_nearest_zero-min_nearest_zero) - a*min_nearest_zero.^2./abs(min_nearest_zero-max_nearest_zero);


% Evaluate the four combinations of endpoint lower and upper bounds.
% The first suffix refers to the left endpoint, the second to the right.

left_lower = intval(inf(left_value));
left_upper = intval(sup(left_value));

right_lower = intval(inf(right_value));
right_upper = intval(sup(right_value));

mixed_corner_ll = b*max(left_lower,right_lower).^2./abs(left_lower-right_lower) - a*min(left_lower,right_lower).^2./abs(left_lower-right_lower);
mixed_corner_lu = b*max(left_lower,right_upper).^2./abs(left_lower-right_upper) - a*min(left_lower,right_upper).^2./abs(left_lower-right_upper);
mixed_corner_ul = b*max(left_upper,right_lower).^2./abs(left_upper-right_lower) - a*min(left_upper,right_lower).^2./abs(left_upper-right_lower);
mixed_corner_uu = b*max(left_upper,right_upper).^2./abs(left_upper-right_upper) - a*min(left_upper,right_upper).^2./abs(left_upper-right_upper);

mixed_upper_bound = max(max(mixed_corner_ll,mixed_corner_lu),max(mixed_corner_ul,mixed_corner_uu));

% Include constant functions at corners with equal endpoint values.
equal_left_lower_corner = (left_lower == right_lower) | (left_lower == right_upper);

equal_left_upper_corner = (left_upper == right_lower) | (left_upper == right_upper);

equal_corner_bound = 2*max(equal_left_lower_corner.*max(a*left_lower, b*left_lower), equal_left_upper_corner.*max(a*left_upper, b*left_upper));

mixed_upper_bound = max(mixed_upper_bound, equal_corner_bound);

% Limit the upper bound using convexity and endpoint support bounds.
mixed_upper_bound = min(mixed_upper_bound, max(a*left_lower, b*left_upper) + max(a*right_lower, b*right_upper));

% Final addition to the integrals

integrals(mask_mix) = hull(intval(0), hull(mixed_auxiliary_value(mask_mix), mixed_upper_bound(mask_mix)));

integrals = h/2*integrals;
end




function updated_block = update_block(power_block, block_propagator, block_exponent, block_index, n)
    block_size = 2^block_exponent;

    next_start = (block_index+1) * block_size;
    next_end = min(next_start + block_size - 1, n - 1);
    num_valid_columns = next_end - next_start + 1;

    % Update the columns needed for the next block.
    updated_block = block_propagator * power_block(:, 1:num_valid_columns);
end

function [power_block, block_exponent, block_propagator] = create_first_block(propagator, x, n)
    state_dimension = size(propagator, 1);
    bytes_per_entry = 16;

    info = memory;
    memory_budget_bytes = floor(0.5 * info.MemAvailableAllArrays);

    max_columns = floor(memory_budget_bytes / (state_dimension * bytes_per_entry));

    % Limit the block size using the memory estimate and a cap of 2048 columns.
    block_exponent =  min(11,min(floor(log2(max_columns)), floor(log2(n)) + 1));
    block_size = 2^block_exponent;

    % Initialize each column with propagator*x.
    power_block = intval(repmat(propagator*x, 1, block_size));

    % Apply the remaining powers using binary decomposition.
    binary_power = intval(propagator);
    for k = 1:block_exponent
        mask = bitget(0:block_size-1, k); % Binary digits of the remaining exponent
        idx = find(mask);
        power_block(:, idx) = binary_power * power_block(:, idx);
        binary_power = binary_power^2;
    end

    block_propagator = binary_power; % Enclosure of propagator^block_size.
end