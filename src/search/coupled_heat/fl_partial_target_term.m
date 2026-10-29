function target_term = fl_partial_target_term(p_f,y_f)
Nx = fix(length(p_f)/2)+1;

if norm(p_f(Nx:2*Nx-2)) ~= 0 
    target_term = Inf;
else
    target_term = transpose(p_f)*y_f/Nx;
end

end