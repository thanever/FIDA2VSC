function [x, f, eflag, output, lambda] = qps_mosek(H, c, A, l, u, xmin, xmax, x0, opt)
%QPS_MOSEK  Quadratic Program Solver based on MOSEK.
%   [X, F, EXITFLAG, OUTPUT, LAMBDA] = ...
%       QPS_MOSEK(H, C, A, L, U, XMIN, XMAX, X0, OPT)
%   A wrapper function providing a MATPOWER standardized interface for using
%   MOSEKOPT to solve the following QP (quadratic programming) problem:
%
%       min 1/2 X'*H*X + C'*X
%        X
%
%   subject to
%
%       L <= A*X <= U       (linear constraints)
%       XMIN <= X <= XMAX   (variable bounds)
%
%   Inputs (all optional except H, C, A and L):
%       H : matrix (possibly sparse) of quadratic cost coefficients
%       C : vector of linear cost coefficients
%       A, L, U : define the optional linear constraints. Default
%           values for the elements of L and U are -Inf and Inf,
%           respectively.
%       XMIN, XMAX : optional lower and upper bounds on the
%           X variables, defaults are -Inf and Inf, respectively.
%       X0 : optional starting value of optimization vector X
%       OPT : optional options structure with the following fields,
%           all of which are also optional (default values shown in
%           parentheses)
%           verbose (0) - controls level of progress output displayed
%               0 = no progress output
%               1 = some progress output
%               2 = verbose progress output
%           max_it (0) - maximum number of iterations allowed
%               0 = use algorithm default
%           mosek_opt - options struct for MOSEK, values in
%               verbose and max_it override these options
%       PROBLEM : The inputs can alternatively be supplied in a single
%           PROBLEM struct with fields corresponding to the input arguments
%           described above: H, c, A, l, u, xmin, xmax, x0, opt
%
%   Outputs:
%       X : solution vector
%       F : final objective function value
%       EXITFLAG : exit flag
%             1 = success
%             0 = terminated at maximum number of iterations
%            -1 = primal or dual infeasible
%           < 0 = the negative of the MOSEK return code
%       OUTPUT : output struct with the following fields:
%           r - MOSEK return code
%           res - MOSEK result struct
%       LAMBDA : struct containing the Langrange and Kuhn-Tucker
%           multipliers on the constraints, with fields:
%           mu_l - lower (left-hand) limit on linear constraints
%           mu_u - upper (right-hand) limit on linear constraints
%           lower - lower bound on optimization variables
%           upper - upper bound on optimization variables
%
%   Note the calling syntax is almost identical to that of QUADPROG
%   from MathWorks' Optimization Toolbox. The main difference is that
%   the linear constraints are specified with A, L, U instead of
%   A, B, Aeq, Beq.
%
%   Calling syntax options:
%       [x, f, exitflag, output, lambda] = ...
%           qps_mosek(H, c, A, l, u, xmin, xmax, x0, opt)
%
%       x = qps_mosek(H, c, A, l, u)
%       x = qps_mosek(H, c, A, l, u, xmin, xmax)
%       x = qps_mosek(H, c, A, l, u, xmin, xmax, x0)
%       x = qps_mosek(H, c, A, l, u, xmin, xmax, x0, opt)
%       x = qps_mosek(problem), where problem is a struct with fields:
%                       H, c, A, l, u, xmin, xmax, x0, opt
%                       all fields except 'c', 'A' and 'l' or 'u' are optional
%       x = qps_mosek(...)
%       [x, f] = qps_mosek(...)
%       [x, f, exitflag] = qps_mosek(...)
%       [x, f, exitflag, output] = qps_mosek(...)
%       [x, f, exitflag, output, lambda] = qps_mosek(...)
%
%   Example: (problem from from http://www.jmu.edu/docs/sasdoc/sashtml/iml/chap8/sect12.htm)
%       H = [   1003.1  4.3     6.3     5.9;
%               4.3     2.2     2.1     3.9;
%               6.3     2.1     3.5     4.8;
%               5.9     3.9     4.8     10  ];
%       c = zeros(4,1);
%       A = [   1       1       1       1;
%               0.17    0.11    0.10    0.18    ];
%       l = [1; 0.10];
%       u = [1; Inf];
%       xmin = zeros(4,1);
%       x0 = [1; 0; 0; 1];
%       opt = struct('verbose', 2);
%       [x, f, s, out, lambda] = qps_mosek(H, c, A, l, u, xmin, [], x0, opt);
%
%   See also MOSEKOPT.

%   MATPOWER
%   $Id: qps_mosek.m,v 1.5 2011/09/09 15:26:08 cvs Exp $
%   by Ray Zimmerman, PSERC Cornell
%   Copyright (c) 2010 by Power System Engineering Research Center (PSERC)
%
%   This file is part of MATPOWER.
%   See http://www.pserc.cornell.edu/matpower/ for more info.
%
%   MATPOWER is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published
%   by the Free Software Foundation, either version 3 of the License,
%   or (at your option) any later version.
%
%   MATPOWER is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
%   GNU General Public License for more details.
%
%   You should have received a copy of the GNU General Public License
%   along with MATPOWER. If not, see <http://www.gnu.org/licenses/>.
%
%   Additional permission under GNU GPL version 3 section 7
%
%   If you modify MATPOWER, or any covered work, to interface with
%   other modules (such as MATLAB code and MEX-files) available in a
%   MATLAB(R) or comparable environment containing parts covered
%   under other licensing terms, the licensors of MATPOWER grant
%   you additional permission to convey the resulting work.

%% check for Optimization Toolbox
% if ~have_fcn('mosek')
%     error('qps_mosek: requires MOSEK');
% end

%%----- input argument handling  -----
%% gather inputs
if nargin == 1 && isstruct(H)       %% problem struct
    p = H;
else                                %% individual args
    p = struct('H', H, 'c', c, 'A', A, 'l', l, 'u', u);
    if nargin > 5
        p.xmin = xmin;
        if nargin > 6
            p.xmax = xmax;
            if nargin > 7
                p.x0 = x0;
                if nargin > 8
                    p.opt = opt;
                end
            end
        end
    end
end

%% define nx, set default values for H and c
if ~isfield(p, 'H') || isempty(p.H) || ~any(any(p.H))
    if (~isfield(p, 'A') || isempty(p.A)) && ...
            (~isfield(p, 'xmin') || isempty(p.xmin)) && ...
            (~isfield(p, 'xmax') || isempty(p.xmax))
        error('qps_mosek: LP problem must include constraints or variable bounds');
    else
        if isfield(p, 'A') && ~isempty(p.A)
            nx = size(p.A, 2);
        elseif isfield(p, 'xmin') && ~isempty(p.xmin)
            nx = length(p.xmin);
        else    % if isfield(p, 'xmax') && ~isempty(p.xmax)
            nx = length(p.xmax);
        end
    end
    p.H = sparse(nx, nx);
    qp = 0;
else
    nx = size(p.H, 1);
    qp = 1;
end
if ~isfield(p, 'c') || isempty(p.c)
    p.c = zeros(nx, 1);
end
if ~isfield(p, 'x0') || isempty(p.x0)
    p.x0 = zeros(nx, 1);
end

%% default options
if ~isfield(p, 'opt')
    p.opt = [];
end
if ~isempty(p.opt) && isfield(p.opt, 'verbose') && ~isempty(p.opt.verbose)
    verbose = p.opt.verbose;
else
    verbose = 0;
end
if ~isempty(p.opt) && isfield(p.opt, 'max_it') && ~isempty(p.opt.max_it)
    max_it = p.opt.max_it;
else
    max_it = 0;
end
if ~isempty(p.opt) && isfield(p.opt, 'mosek_opt') && ~isempty(p.opt.mosek_opt)
    mosek_opt = mosek_options(p.opt.mosek_opt);
else
    mosek_opt = mosek_options;
end
if max_it
    mosek_opt.MSK_IPAR_INTPNT_MAX_ITERATIONS = max_it;
end
if qp
    mosek_opt.MSK_IPAR_OPTIMIZER = 0;   %% default solver only for QP
end

%% set up problem struct for MOSEK
prob.c = p.c;
if qp
   [prob.qosubi, prob.qosubj, prob.qoval] = find(tril(sparse(p.H)));
end
if isfield(p, 'A') && ~isempty(p.A)
    prob.a = sparse(p.A);
end
if isfield(p, 'l') && ~isempty(p.A)
    prob.blc = p.l;
end
if isfield(p, 'u') && ~isempty(p.A)
    prob.buc = p.u;
end
if isfield(p, 'xmin') && ~isempty(p.xmin)
    prob.blx = p.xmin;
end
if isfield(p, 'xmax') && ~isempty(p.xmax)
    prob.bux = p.xmax;
end

%% A is not allowed to be empty
if ~isfield(prob, 'a') || isempty(prob.a)
    unconstrained = 1;
    prob.a = sparse(1, 1, 1, 1, nx);
    prob.blc = -Inf;
    prob.buc =  Inf;
else
    unconstrained = 0;
end

%%-----  run optimization  -----
if verbose
    methods = {
        'default',
        'interior point',
        '<default>',
        '<default>',
        'primal simplex',
        'dual simplex',
        'primal dual simplex',
        'automatic simplex',
        '<default>',
        '<default>',
        'concurrent'
    };
    if isempty(H) || ~any(any(H))
        lpqp = 'LP';
    else
        lpqp = 'QP';
    end
    % (this code is also in mpver.m)
    % MOSEK Version 6.0.0.93 (Build date: 2010-10-26 13:03:27)
    % MOSEK Version 6.0.0.106 (Build date: 2011-3-17 10:46:54)
%    pat = 'Version (\.*\d)+.*Build date: (\d\d\d\d-\d\d-\d\d)';
    pat = 'Version (\.*\d)+.*Build date: (\d+-\d+-\d+)';
    [s,e,tE,m,t] = regexp(evalc('mosekopt'), pat);
    if isempty(t)
        vn = '<unknown>';
    else
        vn = t{1}{1};
    end
    fprintf('MOSEK Version %s -- %s %s solver\n', ...
            vn, methods{mosek_opt.MSK_IPAR_OPTIMIZER+1}, lpqp);
end
cmd = sprintf('minimize echo(%d)', verbose);
[r, res] = mosekopt(cmd, prob, mosek_opt);

%%-----  repackage results  -----
if isfield(res, 'sol')
    if isfield(res.sol, 'bas')
        sol = res.sol.bas;
    else
        sol = res.sol.itr;
    end
    x = sol.xx;
else
    sol = [];
    x = [];
end

%%-----  process return codes  -----
if isfield(res, 'symbcon')
    sc = res.symbcon;
else    
    [r2, res2] = mosekopt('symbcon echo(0)');
    sc = res2.symbcon;
end
eflag = -r;
msg = '';
switch (r)
    case sc.MSK_RES_OK
        if ~isempty(sol)
%            if sol.solsta == sc.MSK_SOL_STA_OPTIMAL
            if strcmp(sol.solsta, 'OPTIMAL')
                msg = 'The solution is optimal.';
                eflag = 1;
            else
                eflag = -1;
%                 if sol.prosta == sc.MSK_PRO_STA_PRIM_INFEAS
                if strcmp(sol.prosta, 'PRIMAL_INFEASIBLE')
                    msg = 'The problem is primal infeasible.';
%                 elseif sol.prosta == sc.MSK_PRO_STA_DUAL_INFEAS
                elseif strcmp(sol.prosta, 'DUAL_INFEASIBLE')
                    msg = 'The problem is dual infeasible.';
                else
                    msg = sol.solsta;
                end
            end
        end
    case sc.MSK_RES_TRM_MAX_ITERATIONS
        eflag = 0;
        msg = 'The optimizer terminated at the maximum number of iterations.';
    otherwise
        if isfield(res, 'rmsg') && isfield(res, 'rcodestr')
            msg = sprintf('%s : %s', res.rcodestr, res.rmsg);
        else
            msg = sprintf('MOSEK return code = %d', r);
        end
end

if (verbose || r == 1001) && ~isempty(msg)  %% always alert user if license is expired
    fprintf('%s\n', msg);
end

%%-----  repackage results  -----
if nargout > 1
    if r == 0
        f = p.c' * x;
        if ~isempty(p.H)
            f = 0.5 * x' * p.H * x + f;
        end
    else
        f = [];
    end
    if nargout > 3
        output.r = r;
        output.res = res;
        if nargout > 4
            if isfield(res, 'sol')
                lambda.lower = sol.slx;
                lambda.upper = sol.sux;
                lambda.mu_l  = sol.slc;
                lambda.mu_u  = sol.suc;
                if unconstrained
                    lambda.mu_l  = [];
                    lambda.mu_u  = [];
                end
            else
                lambda = [];
            end
        end
    end
end
