function plot_sine_series(n_points, sine_coeffs)
x = transpose(linspace(0,1,n_points+1));
sampled_values = evaluate_sine_series(x,sine_coeffs);
hold on
plot(x,sampled_values);
end

function to_plot = evaluate_sine_series(x, z_base)

    to_plot = zeros(length(x), 1);
    for i = 1:length(z_base)
        to_plot = to_plot + z_base(i)*sin(i*pi*x);
    end
end