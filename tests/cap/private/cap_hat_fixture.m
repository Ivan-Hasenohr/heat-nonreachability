function [p,base,norm_squared,pairing] = cap_hat_fixture()
% Independent rational reference for the natural spline of a hat on Nx=8.
% See tests/README_CAP.md for its second-derivative construction.
p = [0;0;0;1;0;0;0];
scaled_derivatives = intval([-3;6;-21;78;0;-78;21;-6;3])/97;
base = [intval(p);scaled_derivatives];
norm_squared = intval(71989)/658630;
pairing = intval(731)/7760; % Integral of the P1 hat times its natural spline.
end
