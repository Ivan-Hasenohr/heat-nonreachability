function tests = test_cap_input_contracts
tests = functiontests(localfunctions);
end

function testFloatingArgumentsAreRejected(~)
[~,scalar]=init_heat_fem(8,4,intval('0.01'),false);
[~,coupled]=init_coupled_heat_fem(8,4,intval('0.001'),false);
cap_assert_rejected(@()zero_certificate(scalar,false),'verified args');
cap_assert_rejected(@()zero_certificate(coupled,true),'verified args');
end

function testInconsistentStepsAndMatrixSizesAreRejected(~)
for coupled=[false,true]
    if coupled, args=init_coupled_heat_fem(8,4,intval('0.001'),false);
    else, args=init_heat_fem(8,4,intval('0.01'),false); end
    bad=args; bad{3}=2*args{3};
    cap_assert_rejected(@()zero_certificate(bad,coupled),'Dt must enclose');
    bad=args; bad{10}=2*args{10};
    cap_assert_rejected(@()zero_certificate(bad,coupled),'h must enclose');
    bad=args; bad{4}=args{4}(1:end-1,1:end-1);
    cap_assert_rejected(@()zero_certificate(bad,coupled),'inconsistent size');
    % Wider outward enclosures of the correct steps remain admissible.
    args{3}=infsup(inf(args{3})*.999,sup(args{3})*1.001);
    args{10}=infsup(inf(args{10})*.999,sup(args{10})*1.001);
    zero_certificate(args,coupled);
end
end

function testInvalidWitnessesAreRejected(~)
p=ones(7,1); args=init_heat_fem(8,4,intval('0.01'),false);
cap_assert_rejected(@()certify_heat_nonreachability(p',0,intval(0),args),'column');
p(1)=NaN;
cap_assert_rejected(@()certify_heat_nonreachability(p,0,intval(0),args),'finite');
p(1)=Inf;
cap_assert_rejected(@()certify_heat_nonreachability(p,0,intval(0),args),'finite');
end

function testScalarInitialStateMustBeExactlyZero(~)
p=zeros(7,1); args=init_heat_fem(8,4,intval('0.01'),false);
for value={1,realmin,midrad(0,1e-12),NaN,Inf}
    initial=value{1};
    cap_assert_rejected(@()certify_heat_nonreachability(p,initial,intval(0),args),'zero initial state');
end
end

function testCoupledRadiusAndComponentContracts(~)
p=zeros(14,1); args=init_coupled_heat_fem(8,4,intval('0.001'),false);
zero=intval([0;0]);
cap_assert_rejected(@()certify_coupled_heat_nonreachability(p,zero,zero,args,2),'target_radius');
for value={intval(-1),midrad(0,.01),intval(Inf),intval(NaN),intval([1,2]),.01}
    radius=value{1};
    cap_assert_rejected(@()certify_coupled_heat_nonreachability(p,zero,zero,args,2,radius),'target_radius');
end
cap_assert_rejected(@()certify_coupled_heat_nonreachability(p,intval([0;0;0]),zero,args), ...
    'equally sized column components');
cap_assert_rejected(@()certify_coupled_heat_nonreachability(p,zero,intval([0,0]),args), ...
    'equally sized column components');
% A zero-radius ball is a valid singleton target.
[ok,J]=certify_coupled_heat_nonreachability(p,zero,zero,args,2,intval(0));
assert(~ok && inf(J)==0 && sup(J)==0);
end

function testInitializersDoNotBuildFloatingTablesUnlessRequested(~)
for initializer={@init_heat_fem,@init_coupled_heat_fem}
    f=initializer{1};
    [args,floating]=f(8,4,intval('0.01'));
    assert(isscalar(args{8}) && isnan(args{8}));
    assert(isscalar(floating{8}) && isnan(floating{8}));
    [~,floating]=f(8,4,intval('0.01'),true);
    dimension=size(floating{4},1);
    assert(isequal(size(floating{8}),[dimension,dimension,5]));
    assert(isequal(floating{8}(:,:,1),eye(dimension)));
    cap_assert_rejected(@()f(8,4,.01,false),'final time');
end
end

function zero_certificate(args,coupled)
if coupled
    [ok,J]=certify_coupled_heat_nonreachability(zeros(14,1),intval([0;0]),intval([0;0]),args);
else
    [ok,J]=certify_heat_nonreachability(zeros(7,1),0,intval(0),args);
end
assert(~ok && inf(J)==0 && sup(J)==0);
end
