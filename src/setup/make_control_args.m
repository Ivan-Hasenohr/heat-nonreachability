function args = make_control_args(varargin)
%MAKE_CONTROL_ARGS Store arguments in a row cell array without conversion.
% This function performs no validation and no verified computation.
% Each initializer defines the meaning and type of its arguments.
%
% Common layout:
%  1: T          - final time.
%  2: Nt         - number of time steps.
%  3: Dt         - time step, T/Nt.
%  4: A          - discrete state generator, size n-by-n.
%  5: B          - nodal control matrix, size n-by-m.
%                 Verified control supports use the control intervals.
%  6:            - scalar heat: sectorial angle alpha;
%                 coupled heat: [a_0,a_1], coercivity and continuity bounds.
%  7: norm_L_T   - control-map norm; -1 means not computed.
%  8:            - implicit-Euler propagation table, or NaN if not built.
%  9: M          - control-size bound used in error estimates.
% 10: h          - spatial mesh size.
% 11: Nx         - number of spatial cells.
%
% Scalar heat:
% 12: kappa      - diffusion coefficient.
% 13:            - control intervals.
%
% Coupled heat:
% 12:            - two diffusion coefficients.
% 13:            - 2-by-2 reaction matrix.
% 14:            - pointwise lower and upper control bounds.
% 15:            - control intervals.
%
% Values are stored unchanged. Wrapping floating data in this container
% does not turn them into verified data.

    args = {};
    args(1:nargin) = varargin;