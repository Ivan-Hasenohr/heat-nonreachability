function [p,p_spline, p_base] = resample_nodal_witness(p_f, Nx)
    p_f_base = iv_natural_cubic_from_nodes(p_f);
    p_f_spline = iv_cubic_hermite_to_polynomials(p_f_base);
    
    x = linspace(1/Nx,(Nx-1)/Nx,Nx-1);
    p = zeros(Nx-1, 1);

    Nx_source = fix(length(p_f_spline)/4);

    for i = 1:Nx-1
        if floor(x(i)*Nx_source) ~= x(i)*Nx_source
            k = 1 + floor(x(i)*Nx_source);
        else 
            k = floor(x(i)*Nx_source);
        end
        p(i) = mid(p_f_spline(4*k-3))*x(i)^3 + mid(p_f_spline(4*k-2))*x(i)^2 + mid(p_f_spline(4*k-1))*x(i) + mid(p_f_spline(4*k));
    end
    p_base = iv_natural_cubic_from_nodes(p);
    p_spline = iv_cubic_hermite_to_polynomials(p_base);
end