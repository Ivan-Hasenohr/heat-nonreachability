function control_support = iv_heat_positive_control_support(p_f,args)
% Enclose the discrete control support sum for 0 <= u <= M.
% Use continuous piecewise affine spatial functions with zero boundary values.
% Sum the time contributions for powers C^q, q = 1,...,Nt.
%
% p_f: interior nodal coefficients, with size (Nx-1)-by-1.
% args{2}: Nt; args{3}: Dt; args{4}: discrete evolution matrix.
% args{9}: M; args{10}: Dx; args{11}: Nx.
% args{13}: sorted, disjoint control intervals encoded by interval endpoints.
%
% Continuous-time and spatial discretization error bounds are handled elsewhere.

Nt = args{2};
Dt = args{3};
Dx = args{10};
Nx = args{11};
control_intervals =  args{13};
if any(inf(control_intervals(2:end)) - sup(control_intervals(1:end-1)) <= 2*Dx) || any(rad(control_intervals) <= Dx/2)
    warning("The control intervals may be too narrow or too close for this mesh.")
end

%% Computations of S_T^star (p_f), first block

A_t = args{4}';
C = verifyinv(eye(Nx-1) - Dt*A_t);

% Build the first block [C*p_f, ..., C^block_size*p_f].
% Cr encloses C^block_size and advances each block to the next one.

[adjoint_block, r, Cr] = create_first_block(C, p_f, Nt);
block_size = 2^r;


%% Loop 

num_blocks = ceil(Nt / block_size);

control_support = 0;

for i = 0:num_blocks-1
    % Column positions are zero-based; position q represents C^(q+1)*p_f.
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

% Integrate the positive part of the piecewise affine function.
% Each control interval is split into interior cells and two boundary cells.
% Boundary cells are excluded from the interior mask.
% Homogeneous Dirichlet values are added at x = 0 and x = 1.

[~,n] = size(adjoint_values);

right_value = intval(zeros(Nx,n));
right_value(1:Nx-1,:) = adjoint_values; 

left_value = intval(zeros(Nx,n));
left_value(2:Nx,:) = adjoint_values; 

% First case: mesh-intervals strictly in B-intervals

% Creates and fills the mask for mesh-intervals strictly in B-intervals
interior_cell_mask = zeros(size(left_value));

for i=1:length(control_intervals)
    control_intervals(i) = max(0, min(1,control_intervals(i)));

    j = max(1, floor(inf(intval(inf(control_intervals(i)))*intval(Nx))) + 1);
    k = min(Nx, ceil(sup(intval(sup(control_intervals(i)))*intval(Nx))));

    cell_indices = repmat((1:Nx)', 1, n);

    interior_cell_mask = interior_cell_mask | ((j < cell_indices) & (cell_indices < k));
end

% Computation of the integral for the first case
% Requires many masks to treat the cases where the values are positive,
% negative or mixed


max_value = max(left_value, right_value); % max of left and right value of mesh-intervals 
min_value = min(left_value, right_value); % min of left and right value of mesh-intervals


% Mask creation for mesh-intervals strictly in B-intervals
mask_zero = (sup(min_value) >= 0) & (inf(max_value) <= 0) & interior_cell_mask; % if 0 is in both min and max values 
mask_pos = (1-mask_zero) & (sup(min_value) >= 0) & interior_cell_mask; % % Both endpoint upper bounds are nonnegative; exclude mask_zero.
mask_neg = (1-mask_zero) & (max_value <= 0) & interior_cell_mask; % if not and if both are completely negative or zero
mask_mix = (1-mask_neg) & (1-mask_pos) & (1-mask_zero) & interior_cell_mask; % % Remaining interior cells after excluding the other three masks.

if ~all(mask_zero+mask_pos+mask_neg+mask_mix == interior_cell_mask, "all") %If there's a mistake on the partition of possibilities (shouldn't happen)
    error("The sign masks do not form a partition of the interior-cell mask.")
end

interior_contribution = M*Dt*Dx/2*sum((mask_zero+mask_pos).*max(0,(left_value+right_value)), "all"); % Positive or zero case: max(0, area of trapeze)

% Integrate mixed cells using endpoint bounds.
interior_contribution = interior_contribution + M*Dt*integral_pos_interval(left_value(mask_mix), right_value(mask_mix), Dx);


% Second case: B-interval boundaries: 
%   - for both inf and sup boundaries, we calculate in which mesh-invertal
%       it is, and calculate accordingly

boundary_contribution = 0;

for i=1:length(control_intervals) %for every B-interval
    
    % left boundary of B-interval
    interval_left = inf(control_intervals(i)) ;
    interval_right = sup(control_intervals(i));

    j = max(1, floor(inf(intval(interval_left)*intval(Nx))) + 1); % Left boundary cell

    % Left and right values on the inf interval:
    left_cell_left_value = left_value(j,:);
    left_cell_right_value = right_value(j,:);

    k = min(Nx, ceil(sup(intval(interval_right)*intval(Nx)))); % Right boundary cell

    % Left and right values on the sup interval:
    right_cell_left_value = left_value(k,:);
    right_cell_right_value = right_value(k,:);
    
    if any(j > k) 
        error("The left boundary cell index exceeds the right boundary cell index.")
    elseif any(j == k)
        error("Control intervals contained in a single mesh cell are not supported.")
    end
    
    % values on the inf and sup boundaries of the B-interval
    value_at_interval_left = (left_cell_left_value*(j*Dx-interval_left) + left_cell_right_value*(interval_left-(j-1)*Dx))/Dx;
    value_at_interval_right = (right_cell_left_value*(k*Dx-interval_right) + right_cell_right_value*(interval_right-(k-1)*Dx))/Dx;
    
    %First the inf case:

    if ~(inf(j*Dx) <= interval_left && interval_left <= sup(j*Dx)) % If the B-inf is not already on the rightside
        left_boundary_contribution = integral_pos_interval(value_at_interval_left, left_cell_right_value,j*Dx -interval_left); %Compute the inegral
    else
        % Bound the boundary contribution by a rectangle over the remaining width.
        left_boundary_contribution = infsup(0, sup(sum(max(intval(0), j*Dx-interval_left) .* max(intval(0), max(left_cell_left_value, left_cell_right_value)), "all")));
    end

    %Second the sup case:

    if ~(inf((k-1)*Dx) <= interval_right && interval_right <= sup((k-1)*Dx))% If the B-sup is not already on the leftside
        right_boundary_contribution = integral_pos_interval(right_cell_left_value, value_at_interval_right, interval_right-(k-1)*Dx);%Compute the inegral
    else
        % Bound the boundary contribution by a rectangle over the remaining width.
        right_boundary_contribution = infsup(0, sup(sum(max(intval(0), interval_right-(k-1)*Dx) .* max(intval(0), max(right_cell_left_value, right_cell_right_value)), "all")));
    end
    
    boundary_contribution = boundary_contribution + right_boundary_contribution + left_boundary_contribution; %Sum everything
    
 end
% Combine the interior and boundary contributions.
block_support = interior_contribution + M*Dt*boundary_contribution;
end


function val = integral_pos_interval(a,b, dx)
%Computes, the value of  int_0^Dx (a + (b-a)*x/Dx)_+ dx where x_+ denote the
%positive part of x. a and b can be (should be) vectors. In which case the
%returned value will be the sum of all integrals

a = intval(a);
b = intval(b);

%Taking the inf and sup of both a and b
maxi_sup = intval(max(sup(a),sup(b)));
mini_sup = intval(min(sup(a),sup(b)));

maxi_inf = intval(max(inf(a),inf(b)));
mini_inf = intval(min(inf(a),inf(b)));

% Lower case first ( inf(a), inf(b) )
mask_pospos = mini_inf >= 0;
mask_posneg = ~(maxi_inf <= 0) & ~mask_pospos;

% Inactive quotient terms may produce NaN when the endpoint bounds coincide.
% The sums below use omitnan.
val_inf = dx/2*(sum((maxi_inf+mini_inf).*mask_pospos, "omitnan") + sum(max(0,maxi_inf.^2./(maxi_inf-mini_inf)).*mask_posneg, "omitnan"));

%Higher case next (sup(a), sup(b))
mask_pospos = mini_sup >= 0;
mask_posneg = ~(maxi_sup <= 0) & ~mask_pospos;

val_sup = dx/2*(sum(sum((maxi_sup+mini_sup).*mask_pospos, "omitnan") + sum(max(0,maxi_sup.^2./(maxi_sup-mini_sup)).*mask_posneg, "omitnan")));

val = infsup(inf(val_inf),sup(val_sup));

end


function AX_updated = update_block(AX_block, block_propagator, block_exponent, i, num_steps)
    block_size = 2^block_exponent;

    next_start = (i+1) * block_size;
    next_end = min(next_start + block_size - 1, num_steps - 1);
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

    % Limit the block size using the memory estimate and a cap of 2048 columns.
    block_exponent =  min(11,min(floor(log2(max_cols)), floor(log2(num_steps)) + 1));
    block_size = 2^block_exponent;

    % Initialize each column with A*x.
    AX_block = intval(repmat(A*x, 1, block_size));

    % Apply the remaining powers using binary decomposition.
    Ak = intval(A);
    for k = 1:block_exponent
        mask = bitget(0:block_size-1, k); % Binary digits of the remaining exponent
        idx = find(mask);
        AX_block(:, idx) = Ak * AX_block(:, idx);
        Ak = Ak^2;
    end

    block_propagator = Ak; % Enclosure of A^block_size.
end