function tests = test_cap_numerics
tests = functiontests(localfunctions);
end

function testScalarEulerPowers(~)
% Exact backward-Euler factor: (1-(1/8)*(-2))^(-1)=4/5.
for steps = [2,3,7,8,9,15,16,17]
    args = cell(1,4);
    args{2} = steps; args{3} = intval(1)/8; args{4} = intval(-2);
    expected = 3*(intval(4)/5)^steps;
    cap_assert_close(iv_euler_state_at_T(intval(3),args),expected);
    cap_assert_close(iv_euler_adjoint_at_T(intval(3),args),expected);
    powers = iv_euler_power_table(steps,args{3},args{4});
    assert(isequal(size(powers),[1,1,steps+1]));
    for exponent = 0:steps
        cap_assert_close(powers(:,:,exponent+1),(intval(4)/5)^exponent);
    end
end
end

function testStateAndAdjointAreDistinct(~)
for steps = [3,8,9,17]
    args = cell(1,4);
    args{2} = steps; args{3} = intval(1)/2;
    args{4} = intval([-2,1;0,-2]); p = intval([1;-2]);
    % Closed form for this triangular resolvent, not a production power call.
    power = (intval(1)/2)^steps*[intval(1),intval(steps)/4;intval(0),intval(1)];
    cap_assert_close(iv_euler_state_at_T(p,args),power*p);
    cap_assert_close(iv_euler_adjoint_at_T(p,args),power'*p);
    powers = iv_euler_power_table(steps,args{3},args{4});
    cap_assert_close(powers(:,:,steps+1),power);
end
end

function testPolynomialIntegrals(~)
q = cap_quadratic(4); hat = intval([0;1;0]);
cap_assert_close(iv_p1_inner_product(hat,hat),intval(1)/6);
cap_assert_close(iv_p1_cubic_hermite_inner_product(hat,q),intval(23)/384);
cap_assert_close(iv_cubic_hermite_inner_product(q,q),intval(1)/30);
cap_assert_close(iv_cubic_hermite_d2_norm(q),intval(2));
end

function testNaturalSplineAgainstRationalReference(~)
[p,base,norm_squared,pairing] = cap_hat_fixture();
cap_assert_close(iv_natural_cubic_from_nodes(p),base);
cap_assert_close(iv_cubic_hermite_inner_product(base,base),norm_squared);
cap_assert_close(iv_cubic_hermite_d2_norm(base),sqrt(intval(712704)/97));
cap_assert_close(iv_p1_cubic_hermite_inner_product(intval(p),base),pairing);
end

function testSinePairings(~)
hat = intval([0;0;1;0;0]); % Includes the Dirichlet endpoint values.
for mode = [1,2,3]
    frequency = intval(mode)*intval('pi');
    if mod(mode,2)==0
        expected = intval(0);
    else
        expected = 8*(1-cos(frequency/4))*sin(frequency/2)/frequency^2;
    end
    coefficients = intval(zeros(mode,1)); coefficients(mode) = 1;
    cap_assert_close(iv_p1_sine_series_pairing(hat,coefficients,intval(1)/4),expected);
    cap_assert_close(iv_p1_sine_series_pairing_vectorized(hat,coefficients,intval(1)/4),expected);
end
end
