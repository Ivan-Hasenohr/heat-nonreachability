function tests = test_cap_certificates
tests = functiontests(localfunctions);
end

function testScalarTargetRepresentations(~)
[p,base,norm_squared,pairing]=cap_hat_fixture();
args=init_heat_fem(8,8,intval('0.01'),false); args{9}=intval(0);
frequency=intval('pi');
sine_pair=16*(1-cos(frequency/8))/frequency^2;
targets={intval(p),base,intval(1)};
expected_pair={pairing,norm_squared,sine_pair};
for representation=1:3
    [ok,J,support,actual_pair,err,b]=certify_heat_nonreachability(p,0,targets{representation},args);
    cap_assert_close(support,intval(0));
    cap_assert_close(actual_pair,expected_pair{representation});
    cap_assert_close(b.discrete_enclosure,-expected_pair{representation});
    if representation<=2
        assert(err==0); % Pairing evaluated directly on the spline witness.
    else
        expected=sqrt(intval(712704)/97)*args{10}^2/(sqrt(intval(2))*frequency^2);
        assert(abs(err-sup(expected))<1e-9);
    end
    assert(ok); check_metadata(ok,J,b);
end
end

function testCoupledTargetRepresentations(~)
[p,base,norm_squared,pairing]=cap_hat_fixture(); args=no_control_coupled();
frequency=intval('pi'); sine_pair=16*(1-cos(frequency/8))/frequency^2;
targets={intval(p),base,intval(1)}; expected={pairing,norm_squared,sine_pair};
radius=intval(1)/100;
for representation=1:3
    target=targets{representation};
    [ok,J,b]=certify_coupled_heat_nonreachability([p;2*p], ...
        intval([0;0]),[2*target;3*target],args,2,radius);
    cap_assert_close(b.discrete_enclosure,radius*sqrt(5*norm_squared)-8*expected{representation});
    assert(ok); check_metadata(ok,J,b);
    [ok,J,b]=certify_coupled_heat_nonreachability([p;0*p], ...
        intval([0;0]),[target;0*target],args,1);
    cap_assert_close(b.discrete_enclosure,-expected{representation});
    assert(ok); check_metadata(ok,J,b);
end
end

function testBallRadiusUsesContinuousJointNorm(~)
[p,base,norm_squared]=cap_hat_fixture(); args=no_control_coupled();
target=[2*base;3*base]; zero=intval([0;0]);
[~,J1,b1]=certify_coupled_heat_nonreachability([p;2*p],zero,target,args,2,intval(1)/100);
[~,J2,b2]=certify_coupled_heat_nonreachability([p;2*p],zero,target,args,2,intval(1)/20);
cap_assert_close(J2-J1,(intval(1)/20-intval(1)/100)*sqrt(5*norm_squared));
assert(isequal(b1.discretization,b2.discretization));
% An interval radius must also cover the endpoint-radius evaluations.
[~,J]=certify_coupled_heat_nonreachability([p;2*p],zero,target,args,2,infsup(.01,.05));
assert(inf(J)<=sup(J1) && sup(J)>=inf(J2));
end

function testInitialStateSignAndRepresentations(~)
[p,base,norm_squared]=cap_hat_fixture(); args=no_control_coupled();
mass=intval(1)/48*toeplitz([4,1,zeros(1,5)]);
stiffness=intval(8)*toeplitz([2,-1,zeros(1,5)]);
R=args{13}; k=args{12}; dt=args{3};
left=[mass+dt*(k(1)*stiffness-R(1,1)*mass),-dt*R(1,2)*mass; ...
    -dt*R(2,1)*mass,mass+dt*(k(2)*stiffness-R(2,2)*mass)];
resolvent=verifylss(left,blkdiag(mass,mass));
targets={intval([p;3*p]),[base;3*base],intval([1;3])};
weights={intval([0;0;1;4;1;0;0])/48, ...
    intval([12;-48;277;1462;277;-48;12])/15520, ...
    16*(1-cos(intval('pi')/8))/intval('pi')^2 ...
        *sin(intval('pi')*intval((1:7)')/8)};
initial_norms={sqrt(intval(10)/12),sqrt(10*norm_squared),sqrt(intval(5))};
for choice=[1,2]
    v=intval([p;2*p]); if choice==1, v(8:end)=0; end
    witness=mid(v);
    for step=1:args{2}, v=resolvent'*v; end
    for representation=1:3
        [ok,J,b]=certify_coupled_heat_nonreachability(witness,targets{representation}, ...
            intval([0;0]),args,choice,intval(0));
        expected=weights{representation}'*(v(1:7)+3*v(8:14));
        cap_assert_close(b.discrete_enclosure,expected);
        % The zero target contributes no projection error. Error scales with ||y0||.
        spline=base; if choice==2, spline=[base;2*base]; else, spline=[base;0*base]; end
        expected_error=iv_coupled_heat_discretization_bound(spline,initial_norms{representation},0,args);
        assert(abs(b.discretization-expected_error)<1e-8);
        assert(~ok); check_metadata(ok,J,b);
    end
end
end

function testNegativeDiscreteValueCanBeInconclusive(~)
[p,~,~,pairing]=cap_hat_fixture();
args=init_heat_fem(8,8,intval('0.01'),false);
support=iv_heat_positive_control_support(p,args);
amplitude=sup(support)/inf(pairing)+.001;
[ok,J,~,~,~,b]=certify_heat_nonreachability(p,0,amplitude*intval(p),args);
assert(sup(b.discrete_enclosure)<0 && sup(J)>0 && ~ok);
check_metadata(ok,J,b);
end

function testZeroPositiveAndInfiniteFinalValuesAreNotCertificates(~)
[p,base]=cap_hat_fixture();
scalar=init_heat_fem(8,8,intval('0.01'),false); scalar{9}=intval(0);
[ok,J]=certify_heat_nonreachability(0*p,0,intval(0),scalar);
assert(~ok && inf(J)==0 && sup(J)==0);
[ok,J]=certify_heat_nonreachability(p,0,-intval(p),scalar);
assert(~ok && inf(J)>0);
% The second component of a partial-target witness must vanish exactly.
vector=[p;0*p]; vector(8)=1e-12;
[ok,J]=certify_coupled_heat_nonreachability(vector,intval([0;0]), ...
    [base;0*base],no_control_coupled(),1);
assert(~ok && isinf(sup(J)));
end

function testMetadataSurvivesSaveLoad(~)
[p,base]=cap_hat_fixture();
[ok,J,bounds]=certify_coupled_heat_nonreachability([p;0*p], ...
    intval([0;0]),[base;0*base],no_control_coupled());
file=[tempname,'.mat']; cleanup=onCleanup(@()delete(file)); %#ok<NASGU>
save(file,'ok','J','bounds'); loaded=load(file);
assert(isequal(ok,loaded.ok) && isequaln(J,loaded.J) && isequaln(bounds,loaded.bounds));
check_metadata(loaded.ok,loaded.J,loaded.bounds);
end

function args = no_control_coupled()
args=init_coupled_heat_fem(8,4,intval('0.001'),false);
args{14}=intval([0,0]);
end

function check_metadata(ok,J,b)
% Bookkeeping checks complement, rather than replace, the analytic oracles.
assert(isequaln(J,b.final_enclosure));
assert(isequaln(b.roundoff,rad(b.discrete_enclosure)));
assert(isequaln(J,b.discrete_enclosure+midrad(0,b.discretization)));
assert(b.time>=0 && b.space>=0 && b.discretization>=0);
assert(logical(ok)==(isfinite(inf(J)) && isfinite(sup(J)) && sup(J)<0));
end
