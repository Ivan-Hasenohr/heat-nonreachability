function v = fl_project_positive_controls(u,~,args)
%FL_PROX_G_CLOHESSY Fonction pour Chambolle-Pock 2011 pour un ensemble de
%contrôle sous la forme 0 <= u(x,t) <= M
%   Detailed explanation goes here
M =  args{9};
v = min(M,max(0,u));
end

