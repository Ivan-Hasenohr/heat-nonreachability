function norm_bound = iv_sine_series_d2_norm_bound(sine_coeffs)
    weighted_sum = intval(0);
    for i = 1:length(sine_coeffs)
        weighted_sum = weighted_sum + intval(sine_coeffs(i))^2*intval(i)^4;
    end
    norm_bound = sup(intval('pi')^2*sqrt(weighted_sum/intval(2)));
end