function target_term = fl_ball_target_term(p_f,y_tilde_T_float, target_radius_fl)

Nx = fix(length(p_f)/2)+1;

target_term = transpose(p_f)*y_tilde_T_float/Nx + target_radius_fl*norm(p_f)/sqrt(Nx);
end