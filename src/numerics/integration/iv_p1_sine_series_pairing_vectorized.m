function pairing = iv_p1_sine_series_pairing_vectorized(p_f, sine_coeffs, Dx)
% Computes the scalar product of a P1 with a sum of sinus

% Use interval column vectors before performing arithmetic.
p_f = intval(p_f(:));
sine_coeffs = intval(sine_coeffs(:));
Dx = intval(Dx);

Nx = length(p_f)-1;
iv_pi = intval('pi');
mode_indices = (1:length(sine_coeffs));
node_indices = (0:Nx)';
pairing = 1/(iv_pi^2*Dx)*sum((p_f(2:end) - p_f(1:end-1)).*sum(sine_coeffs'.*(sin(mode_indices.*iv_pi.*node_indices(2:end)*Dx) - sin(mode_indices.*iv_pi.*node_indices(1:end-1)*Dx))./(mode_indices.^2),2));
% Add the boundary term for each sine mode.
pairing = pairing + 1/iv_pi*sum(sine_coeffs./mode_indices'.*(p_f(1) - (-1).^mode_indices'.*p_f(end)));
end