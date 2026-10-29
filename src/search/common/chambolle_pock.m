function [x,y] = chambolle_pock(prox_F_star, prox_G, K, K_t, sigma, tau, theta, x_0, y_0,continue_search)
%CHAMBOLLE_POCK Implements the 2011 Chambolle-Pock algorithm
%
%Notations are the same that the ones used in the 2011 paper
%Beware the sign of y : it is not the usual convention used in optimal
%control
%
% Approximates a saddle-point of the primal-dual problem : 
%       inf_{x in E} F(Kx) + G(x)       (P)
%       sup_{y in X} -F^*(y) - G^*(-K^*y)     (D)
%
%
% Given : 
%   - prox_F_star : proximal operator for F^*, as defined in Chambolle_Pock_2011
%   - prox_G : proximal operator for G, as defined in Chambolle_Pock_2011
%   - K : operator, as defined in Chambolle_Pock_2011
%   - K_t = K^*, as defined in Chambolle_Pock_2011
%   - sigma, tau, theta : as defined in Chambolle_Pock_2011. Must verify
%   sigma*tau*norm(K)^2 < 1
%   - x_0, y_0 : initial conditions of x and y
%   - continue_search : function returning a boolean determining when the while
%   loop should stop : continue_search(x,y,x_old,y_old,i), where :
%       - i is the number of iterations having been computed
%       - x,y are the numbers computed at the i-th iteration
%       - x,y are the numbers computed at the (i-1)-th iteration
%
% Returns an approximation of a saddle-point of the 
% primal-dual problem (P)&(D). For further information, see
% Chambolle_Pock_2011.

x=x_0;
x_bar = x_0; 
y = y_0;
i = 0;
x_old = x;
y_old = y;
while continue_search(x,y,x_old, y_old,i)    
    x_old = x;
    y_old = y;
    y = prox_F_star(y+sigma*K(x_bar),sigma);
    x = prox_G(x - tau*K_t(y),tau);
    x_bar = x + theta*(x - x_old);
    
    i=i+1;
end
end

