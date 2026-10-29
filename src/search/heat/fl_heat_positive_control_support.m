function control_support = fl_heat_positive_control_support(p,args)


Dt = args{3};
M = args{9};
Dx = args{10};

control_support = M*Dt*Dx*sum(max(0,fl_apply_control_adjoint(p,args)),"all");

end