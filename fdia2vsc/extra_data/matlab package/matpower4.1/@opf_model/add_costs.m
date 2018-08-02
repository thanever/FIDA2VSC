function om = add_costs(om, name, cp, varsets)
%ADD_COSTS  Adds a set of user costs to the model.
%   OM = ADD_COSTS(OM, NAME, CP, VARSETS);
%
%   Adds a named block of user-defined costs to the model. Each set is
%   defined by the CP struct described below. All user-defined sets of
%   costs are combined together into a single set of cost parameters in
%   a single CP struct by BULD_COST_PARAMS. This full aggregate set of
%   cost parameters can be retreived from the model by GET_COST_PARAMS.
%
%   Examples:
%     cp1 = struct('N', N1, 'Cw', Cw1);
%     cp2 = struct('N', N2, 'Cw', Cw2, 'H', H, 'dd', dd, ...
%                   'rh', rh, 'kk', kk, 'mm', mm);
%     om = add_costs(om, 'usr1', cp1, {'Pg', 'Qg', 'z'});
%     om = add_costs(om, 'usr2', cp2, {'Vm', 'Pg', 'Qg', 'z'});
%
%   Let x refer to the vector formed by combining the specified VARSETS,
%   and f_u(x, CP) be the cost at x corresponding to the cost parameters
%   contained in CP, where CP is a struct with the following fields:
%       N      - nw x nx sparse matrix
%       Cw     - nw x 1 vector
%       H      - nw x nw sparse matrix (optional, all zeros by default)
%       dd, mm - nw x 1 vectors (optional, all ones by default)
%       rh, kk - nw x 1 vectors (optional, all zeros by default)
%
%   These parameters are used as follows to compute f_u(x, CP)
%
%       R  = N*x - rh
%
%               /  kk(i),  R(i) < -kk(i)
%       K(i) = <   0,     -kk(i) <= R(i) <= kk(i)
%               \ -kk(i),  R(i) > kk(i)
%
%       RR = R + K
%
%       U(i) =  /  0, -kk(i) <= R(i) <= kk(i)
%               \  1, otherwise
%
%       DDL(i) = /  1, dd(i) = 1
%                \  0, otherwise
%
%       DDQ(i) = /  1, dd(i) = 2
%                \  0, otherwise
%
%       Dl = diag(mm) * diag(U) * diag(DDL)
%       Dq = diag(mm) * diag(U) * diag(DDQ)
%
%       w = (Dl + Dq * diag(RR)) * RR
%
%       f_u(x, CP) = 1/2 * w'*H*w + Cw'*w
%
%   See also OPF_MODEL, BUILD_COST_PARAMS, GET_COST_PARAMS, COMPUTE_COST.

%   MATPOWER
%   $Id: add_costs.m,v 1.7 2010/04/26 19:45:25 ray Exp $
%   by Ray Zimmerman, PSERC Cornell
%   Copyright (c) 2008-2010 by Power System Engineering Research Center (PSERC)
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

%% prevent duplicate named cost sets
if isfield(om.cost.idx.N, name)
    error('@opf_model/add_costs: cost set named ''%s'' already exists', name);
end

if nargin < 4
    varsets = {};
end
if isempty(varsets)
    varsets = om.var.order;
end
[nw, nx] = size(cp.N);

%% check sizes
nv = 0;
for k = 1:length(varsets)
    nv = nv + om.var.idx.N.(varsets{k});
end
if nx ~= nv
    if nw == 0
        cp.N = sparse(nw, nx);
    else
        error('@opf_model/add_costs: number of columns in N (%d x %d) does not match\nnumber of variables (%d)\n', nw, nx, nv);
    end
end
if size(cp.Cw, 1) ~= nw
    error('@opf_model/add_costs: number of rows of Cw (%d x %d) and N (%d x %d) must match\n', size(cp.Cw), nw, nx);
end
if isfield(cp, 'H') && (size(cp.H, 1) ~= nw || size(cp.H, 2) ~= nw)
    error('@opf_model/add_costs: both dimensions of H (%d x %d) must match the number of rows in N (%d x %d)\n', size(cp.H), nw, nx);
end
if isfield(cp, 'dd') && size(cp.dd, 1) ~= nw
    error('@opf_model/add_costs: number of rows of dd (%d x %d) and N (%d x %d) must match\n', size(cp.dd), nw, nx);
end
if isfield(cp, 'rh') && size(cp.rh, 1) ~= nw
    error('@opf_model/add_costs: number of rows of rh (%d x %d) and N (%d x %d) must match\n', size(cp.rh), nw, nx);
end
if isfield(cp, 'kk') && size(cp.kk, 1) ~= nw
    error('@opf_model/add_costs: number of rows of kk (%d x %d) and N (%d x %d) must match\n', size(cp.kk), nw, nx);
end
if isfield(cp, 'mm') && size(cp.mm, 1) ~= nw
    error('@opf_model/add_costs: number of rows of mm (%d x %d) and N (%d x %d) must match\n', size(cp.mm), nw, nx);
end

%% add info about this user cost set
om.cost.idx.i1.(name)  = om.cost.N + 1;     %% starting index
om.cost.idx.iN.(name)  = om.cost.N + nw;    %% ending index
om.cost.idx.N.(name)   = nw;                %% number of costs (nw)
om.cost.data.N.(name)  = cp.N;
om.cost.data.Cw.(name) = cp.Cw;
om.cost.data.vs.(name) = varsets;
if isfield(cp, 'H')
    om.cost.data.H.(name)  = cp.H;
end
if isfield(cp, 'dd')
    om.cost.data.dd.(name) = cp.dd;
end
if isfield(cp, 'rh')
    om.cost.data.rh.(name) = cp.rh;
end
if isfield(cp, 'kk')
    om.cost.data.kk.(name) = cp.kk;
end
if isfield(cp, 'mm')
    om.cost.data.mm.(name) = cp.mm;
end

%% update number of vars and var sets
om.cost.N  = om.cost.idx.iN.(name);
om.cost.NS = om.cost.NS + 1;

%% put name in ordered list of var sets
om.cost.order{om.cost.NS} = name;
