function plot_cubic_polynomials(splines, x_min, x_max, num_points)
    hold on
    [n,n_spline] = size(splines);
    Nx = fix(n/4);
    x = transpose(linspace(x_min, x_max, num_points));
    
    sampled_values = intval(zeros(num_points, n_spline));
    first_sample_index = 2;
    if x_min == 0
    sampled_values(1,:) = splines(1,:)*x_min^3 + splines(2,:)*x_min^2 + splines(3,:)*x_min + splines(4,:);%mid()
    else 
        first_sample_index = 1;
    end
    for i = first_sample_index:num_points
        if floor(x(i)*Nx) ~= x(i)*Nx
            k = 1 + floor(x(i)*Nx);
        else 
            k = floor(x(i)*Nx);
        end
        sampled_values(i,:) = splines(4*k-3,:)*x(i)^3 + splines(4*k-2,:)*x(i)^2 + splines(4*k-1,:)*x(i) + splines(4*k,:);
    end

    plot(x,sampled_values)
end