function pairing = iv_p1_sine_series_pairing(p_f, sine_coeffs, Dx)
    % Computes the scalar product of a P1 with a sum of sinus
    % Would require vectorising 
    % Promote inputs before differences and products.
    p_f = intval(p_f(:));
    sine_coeffs = intval(sine_coeffs(:));
    Dx = intval(Dx);
    pairing = intval('0');
    for i = 1:length(sine_coeffs) %for each basis value, sums the value of the scalar product between P1 and sin(i*pi*\cdot)
        pairing = pairing + sine_coeffs(i)*pair_p1_with_sine_mode(p_f, i,Dx); 
    end
end

function s = pair_p1_with_sine_mode(p_f,k,Dx)
    s = intval(0);
    for i = 1:length(p_f)-1
        s = s + integrate_affine_sine_cell(p_f(i),p_f(i+1),i-1,k,Dx);
        %intval_int_sin_affine_Dx(k,(p_f(i+1) - p_f(i))/Dx,(i+1)*p_f(i) - i*p_f(i+1),i,Dx);
    end
end

function cell_pairing = integrate_affine_sine_cell(left_value,right_value,cell_index,mode_index,Dx) 

%cell_index starts at 0
iv_pi = intval('pi');
cell_pairing = ((right_value-left_value)*(sin(iv_pi*(cell_index+1)*Dx*mode_index) - sin(iv_pi*cell_index*Dx*mode_index)) + iv_pi*Dx*mode_index*left_value*cos(iv_pi*Dx*cell_index*mode_index) -iv_pi*Dx*mode_index*right_value*cos(iv_pi*(cell_index+1)*Dx*mode_index))/(iv_pi^2*Dx*mode_index^2);
end 