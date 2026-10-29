function tests = test_cap_error_bounds
tests = functiontests(localfunctions);
end

function testScalarPolynomialBudget(~)
% q=x(1-x): ||q''||=2, independently of the Hermite integration code.
args=init_heat_fem(8,16,intval(1)/8,false);
args{12}=intval(2); q=cap_quadratic(8);
ca=2+2/intval('sqrt3');
c2=(7+4*intval('log2')+ca)/(2*intval('pi')^2);
c3=ca/(2*exp(intval(1)));
for kind=[1,2]
    for amplitude=[0,1,2]
        for measure=[1/4,1]
            args{9}=intval(amplitude); args{13}=infsup(0,measure);
            Q=intval(amplitude);
            if kind==1, Q=Q*sqrt(intval(measure)); end
            time=Q*args{1}*args{3}*(intval(1)/2+c3)*4;
            space=Q*args{1}*c2*args{10}^2*4 ...
                +(Q*args{1}*ca+intval(7)/8)*args{10}^2/intval('pi')^2*2;
            [err,b]=iv_heat_discretization_bound(q,0,intval(7)/8,args,kind);
            compare_upper(err,time+space);
            compare_upper(b.time,time); compare_upper(b.space,space);
            assert(isequal(err,b.discretization));
        end
    end
end
end

function testCoupledPolynomialBudget(~)
% Integrate ||diag(kappa)*p''+R'*p||^2 explicitly for p=u*x(1-x).
args=init_coupled_heat_fem(8,16,intval(1)/8,false);
args{14}=intval([-2,2]); args{15}=infsup(0,1/2);
q=cap_quadratic(8); pi2=intval('pi')^2;
reactions={intval([-1,2;-2,-1]),intval([-2,1;1,-2]), ...
    [intval(1)/2,intval(2);intval(0),intval(-1)]};
for scenario=1:3
    args{13}=reactions{scenario};
    if scenario<=2
        a0=intval(1); c0=2/intval('pi');
    else
        a0=1-2/pi2;
        delta2=1-(-1+sqrt(intval(8)))/pi2;
        delta3=1-(-1+sqrt(intval(13)))/pi2;
        c0=(1/sqrt(delta2)+1/sqrt(delta3))/intval('pi');
    end
    a1=intval(3)/2+4/pi2; args{6}=[a0,a1];
    if scenario==2, c1=a1*sqrt(a1)*c0^2/(4*sqrt(a0));
    else, c1=a1^2*c0^2/a0; end
    ca=2+2/intval('sqrt3');
    c2=c1*(7+4*intval('log2')*a1/a0+ca);
    c3=ca*a1/(2*exp(intval(1))*a0);
    for weights={[1;0],[0;1],[1;2],[1;-2]}
        u=intval(weights{1}); reaction=args{13}'*u; k=args{12}(:).*u;
        normA=sqrt(4*sum(k.^2)-intval(2)/3*sum(k.*reaction)+sum(reaction.^2)/30);
        projection=2*norm(u)*args{10}^2/pi2;
        Q=2/sqrt(intval(2)); initial=intval(3)/4; target=intval(5)/4;
        time=args{3}*(Q*args{1}/2+c3*(initial+Q*args{1}))*normA;
        space=args{10}^2*c2*(initial+Q*args{1})*normA ...
            +((initial+Q*args{1})*ca+target)*projection;
        [err,b]=iv_coupled_heat_discretization_bound([u(1)*q;u(2)*q],initial,target,args);
        compare_upper(err,time+space);
        compare_upper(b.time,time); compare_upper(b.space,space);
    end
end
end

function testActualControlBoundsDetermineAmplitude(~)
args=init_coupled_heat_fem(8,4,intval('0.001'),false);
q=cap_quadratic(8); p=[q;2*q];
reference=iv_coupled_heat_discretization_bound(p,0,0,args);
% Cell 9 is legacy data: it must not override the true box constraints.
for stale=[0,1000]
    args{9}=intval(stale);
    actual=iv_coupled_heat_discretization_bound(p,0,0,args);
    assert(isequal(actual,reference));
end
args{14}=intval([-3,1]);
actual=iv_coupled_heat_discretization_bound(p,0,0,args);
assert(abs(actual-1.5*reference)<1e-9);
end

function testProjectionErrorDoesNotScaleWithControlMeasure(~)
q=cap_quadratic(8);
args=init_heat_fem(8,4,intval('0.01'),false); args{9}=intval(0);
first=iv_heat_discretization_bound(q,0,1,args,1);
args{13}=infsup(0,1);
second=iv_heat_discretization_bound(q,0,1,args,1);
assert(isequal(first,second)); compare_upper(first,2*args{10}^2/intval('pi')^2);
coupled=init_coupled_heat_fem(8,4,intval('0.001'),false);
coupled{14}=intval([0,0]);
first=iv_coupled_heat_discretization_bound([q;2*q],0,1,coupled);
coupled{15}=infsup(0,1);
second=iv_coupled_heat_discretization_bound([q;2*q],0,1,coupled);
assert(isequal(first,second));
compare_upper(first,2*sqrt(intval(5))*coupled{10}^2/intval('pi')^2);
end

function testSymmetryRequiresEqualPointCoefficients(~)
args=init_coupled_heat_fem(8,4,intval('0.001'),false);
q=cap_quadratic(8);
off_diagonal={intval(1),intval(2),infsup(.9,1.1)};
other={intval(1),intval(-2),infsup(.9,1.1)};
for scenario=1:3
    args{13}=[intval(-1),off_diagonal{scenario};other{scenario},intval(-1)];
    text=evalc('err=iv_coupled_heat_discretization_bound([q;2*q],0,0,args,true);');
    assert(isfinite(err) && err>=0);
    if scenario==1
        assert(contains(text,'Self-adjoint case'));
    else
        assert(contains(text,'self-adjointness is not certified'));
    end
end
end

function testVerboseIsOptionalAndDoesNotChangeBounds(~)
q=cap_quadratic(8);
scalar=init_heat_fem(8,4,intval('0.01'),false);
silent=evalc('[a,b]=iv_heat_discretization_bound(q,0,0,scalar,1);');
shown=evalc('[c,d]=iv_heat_discretization_bound(q,0,0,scalar,1,true);');
assert(isempty(strtrim(silent)) && contains(shown,'Time error bound'));
assert(isequal(a,c) && isequal(b,d));
coupled=init_coupled_heat_fem(8,4,intval('0.001'),false);
silent=evalc('[a,b]=iv_coupled_heat_discretization_bound([q;2*q],0,0,coupled);');
shown=evalc('[c,d]=iv_coupled_heat_discretization_bound([q;2*q],0,0,coupled,true);');
assert(isempty(strtrim(silent)) && contains(shown,'Time error bound'));
assert(isequal(a,c) && isequal(b,d));
end

function testNonpositiveRegularityDenominatorIsRejected(~)
args=init_coupled_heat_fem(8,4,intval('0.001'),false);
args{13}=100*intval(eye(2)); q=cap_quadratic(8);
cap_assert_rejected(@()iv_coupled_heat_discretization_bound([q;q],0,0,args),'regularity');
end

function compare_upper(actual,expected)
assert(isfinite(actual) && actual>=inf(expected));
assert(abs(actual-sup(expected))<1e-8);
end
