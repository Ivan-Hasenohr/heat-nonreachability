function tests = test_cap_control_supports
tests = functiontests(localfunctions);
end

function testAlignedAndPartialControlCells(~)
% Dyadic endpoints give exact test geometry, including at x=0 and x=1.
domains = {[0,1],[1/8,7/8],[3/32,29/32],[0,3/8;3/4,1]};
fields = {[0;0;0;1;0;0;0],ones(7,1),-[1;2;1;2;1;2;1], ...
    [1;-2;3;-1;2;-3;1]};
for i = 1:length(domains)
    intervals = domains{i};
    for j = 1:length(fields)
        p = fields{j}; args = stationary_args(8,3,intervals);
        [positive,negative,square] = spatial_integrals(p,intervals);
        check_upper(iv_heat_positive_control_support(p,args),positive);
        check_upper(iv_heat_l2_control_support(p,args),sqrt(square));
        % Only the second component is controlled, with bounds [-1,2].
        check_upper(iv_coupled_box_control_support([10*p;p],couple(args)), ...
            2*positive-negative);
    end
end
end

function testPowerOfTwoAndBlockTransitions(~)
% A=-I and Dt=1/Nt give a nonconstant exact geometric time factor.
p = [0;1;0];
for kind = 1:3
    if kind==2, cap = 8192; else, cap = 2048; end
    for steps = [cap-1,cap,cap+1]
        args = stationary_args(4,steps,[0,1]);
        args{4} = -intval(eye(3));
        ratio = 1/(1+args{3});
        time_factor = args{3}*ratio*(1-ratio^steps)/(1-ratio);
        if kind==1
            actual = iv_heat_positive_control_support(p,args);
            expected = time_factor/4;
        elseif kind==2
            actual = iv_heat_l2_control_support(p,args);
            expected = time_factor/sqrt(intval(6));
        else
            args = couple(args); args{4} = -intval(eye(6));
            actual = iv_coupled_box_control_support([0*p;p],args);
            expected = time_factor/2;
        end
        check_upper(actual,expected,1e-7);
    end
end
end

function testControlledComponentAndAdjointCoupling(~)
args = couple(stationary_args(4,3,[0,1])); p = [0;1;0];
I = intval(eye(3)); Z = intval(zeros(3));
args{4} = [Z,I;Z,Z];
% (I-Dt*A')^(-q)*[p;0]=[p;q*Dt*p].
% Sum_q 2*Dt*(q*Dt)*integral(p), with integral(p)=1/4.
expected = args{3}^2*args{2}*(args{2}+1)/4;
check_upper(iv_coupled_box_control_support([p;0*p],args),expected);
end

function testIntervalWitnessesBoundSampledFields(~)
args = stationary_args(8,3,[3/32,29/32]);
center = [1;-2;0;1;-1;2;-1]; radius = ones(7,1)/8;
box = midrad(center,radius);
positive = iv_heat_positive_control_support(box,args);
l2 = iv_heat_l2_control_support(box,args);
box_support = iv_coupled_box_control_support([0*box;box],couple(args));
assert(isfinite(sup(positive)) && isfinite(sup(l2)) && isfinite(sup(box_support)));
for sample = [center-radius,center,center+radius]
    [pos,neg,square] = spatial_integrals(sample,[3/32,29/32]);
    assert(sup(positive)>=inf(pos) && sup(l2)>=inf(sqrt(square)));
    assert(sup(box_support)>=inf(2*pos-neg));
end
end

function testOneCellControlIsRejected(~)
args = stationary_args(8,3,[1/64,3/64]); p = ones(7,1);
cap_assert_rejected(@()iv_heat_positive_control_support(p,args),'mesh cell');
cap_assert_rejected(@()iv_heat_l2_control_support(p,args),'mesh cell');
cap_assert_rejected(@()iv_coupled_box_control_support([p;p],couple(args)),'mesh cell');
end

function args = stationary_args(N,steps,intervals)
% Synthetic kernel fixture: A=0 isolates spatial integration from propagation.
args = cell(1,13);
args{2}=steps; args{3}=intval(1)/steps; args{4}=intval(zeros(N-1));
args{9}=intval(1); args{10}=intval(1)/N; args{11}=N;
args{13}=infsup(intervals(:,1)',intervals(:,2)');
end

function result = couple(args)
result=args; N=args{11};
result{4}=intval(zeros(2*N-2));
result{14}=intval([-1,2]); result{15}=args{13};
end

function check_upper(actual,expected,tolerance)
% Mixed-sign cells may deliberately have a loose lower bound in the kernel.
if nargin<3, tolerance=1e-9; end
assert(isfinite(inf(actual)) && isfinite(sup(actual)));
assert(inf(actual)<=sup(expected) && sup(actual)>=inf(expected));
assert(abs(sup(actual)-sup(expected))<tolerance);
end

function [positive,negative,square] = spatial_integrals(p,intervals)
% Independent cell-by-cell oracle, without production masks or block logic.
N=length(p)+1; values=intval([0;p;0]);
positive=intval(0); negative=intval(0); square=intval(0);
for cell_index=1:N
    for region=1:size(intervals,1)
        left=max((cell_index-1)/N,intervals(region,1));
        right=min(cell_index/N,intervals(region,2));
        if right<=left, continue; end
        t_left=intval(left)*N-(cell_index-1);
        t_right=intval(right)*N-(cell_index-1);
        a=(1-t_left)*values(cell_index)+t_left*values(cell_index+1);
        b=(1-t_right)*values(cell_index)+t_right*values(cell_index+1);
        width=intval(right)-intval(left);
        if inf(a)>=0 && inf(b)>=0
            pos=width*(a+b)/2;
        elseif sup(a)<=0 && sup(b)<=0
            pos=intval(0);
        else
            pos=width*max(a,b)^2/(2*abs(a-b));
        end
        positive=positive+pos;
        negative=negative+width*(a+b)/2-pos;
        square=square+width*(a^2+a*b+b^2)/3;
    end
end
end
