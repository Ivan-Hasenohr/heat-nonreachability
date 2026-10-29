function [p,p_spline, p_base] = resample_polynomial_witness(spline, Nx)
    x = linspace(1/Nx,(Nx-1)/Nx,Nx-1);
    p = zeros(Nx-1, 1);

    Nx_source = fix(length(spline)/4);

    for i = 1:Nx-1
        if floor(x(i)*Nx_source) ~= x(i)*Nx_source
            k = 1 + floor(x(i)*Nx_source);
        else 
            k = floor(x(i)*Nx_source);
        end
        p(i) = mid(spline(4*k-3))*x(i)^3 + mid(spline(4*k-2))*x(i)^2 + mid(spline(4*k-1))*x(i) + mid(spline(4*k));
    end
    p_base = iv_natural_cubic_from_nodes(p);
    p_spline = iv_cubic_hermite_to_polynomials(p_base);
end