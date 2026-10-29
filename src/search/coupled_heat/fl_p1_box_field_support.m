function control_support = fl_p1_box_field_support(u,args)
% Floating support integral of a piecewise affine control-space field.
% u contains nodal adjoint-field values, not primal control values.
% Zero endpoint values are appended before spatial integration.
% The input field already includes the nodal control mask.
% This routine does not integrate directly over the exact control intervals.
% Its value can differ from the nodal-quadrature support used by the search.
% No verified error enclosure is provided.

Dt = args{3};
bounds = args{14};
a = min(bounds); b = max(bounds);
Dx = args{10};

[n,Nt] = size(u);
Nx = n+1;

right_value = zeros(Nx,Nt);
right_value(1:Nx-1,:) = u;

left_value = zeros(Nx,Nt);
left_value(2:Nx,:) = u;

max_value = max(left_value, right_value);
min_value = min(left_value, right_value);

mask_zero = (min_value >= 0).*(max_value <= 0);
mask_pos = (1-mask_zero).*(min_value >= 0);
mask_neg = (1-mask_zero).*(max_value <= 0);
mask_mix = (1-mask_neg).*(1-mask_pos).*(1-mask_zero);

control_support = b*Dt*Dx/2*sum((mask_zero+mask_pos).*(left_value+right_value), "all");
control_support = control_support + a*Dt*Dx/2*sum((mask_neg).*(left_value+right_value), "all");

control_support = control_support + b*Dt*Dx/2*sum(mask_mix.*max_value.^2./(max(max_value-min_value,1e-100)), "all");
control_support = control_support + a*Dt*Dx/2*sum(mask_mix.*min_value.^2./(min(min_value-max_value,-1e-100)), "all");

end