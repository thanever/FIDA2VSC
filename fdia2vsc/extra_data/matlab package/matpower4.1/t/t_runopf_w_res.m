function t_runopf_w_res(quiet)
%T_RUNOPF_W_RES  Tests RUNOPF_W_RES and the associated callbacks.

%   MATPOWER
%   $Id: t_runopf_w_res.m,v 1.11 2011/03/23 18:09:12 cvs Exp $
%   by Ray Zimmerman, PSERC Cornell
%   Copyright (c) 2009-2010 by Power System Engineering Research Center (PSERC)
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

if nargin < 1
    quiet = 0;
end

t_begin(46, quiet);

if quiet
    verbose = 0;
else
    verbose = 0;
end

casefile = 't_case30_userfcns';
mpopt = mpoption('OPF_VIOLATION', 1e-6, 'PDIPM_GRADTOL', 1e-8, ...
        'PDIPM_COMPTOL', 1e-8, 'PDIPM_COSTTOL', 1e-9);
mpopt = mpoption(mpopt, 'OUT_ALL', 0, 'VERBOSE', verbose, 'OPF_ALG', 560);

%% define named indices into data matrices
[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
    MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
    QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;

s6 = warning('query', 'MATLAB:nearlySingularMatrixUMFPACK');
warning('off', 'MATLAB:nearlySingularMatrixUMFPACK');

t = 'runopf_w_res(''t_case30_userfcns'') : ';
r = runopf_w_res(casefile, mpopt);
t_is(r.reserves.R, [25; 15; 0; 0; 19.3906; 0.6094], 4, [t 'R']);
t_is(r.reserves.prc, [2; 2; 2; 2; 5.5; 5.5], 6, [t 'prc']);
t_is(r.reserves.mu.l, [0; 0; 1; 2; 0; 0], 7, [t 'mu.l']);
t_is(r.reserves.mu.u, [0.1; 0; 0; 0; 0; 0], 7, [t 'mu.u']);
t_is(r.reserves.mu.Pmax, [0; 0; 0; 0; 0.5; 0], 7, [t 'mu.Pmax']);
mpc = loadcase(casefile);
t_is(r.reserves.cost, mpc.reserves.cost, 12, [t 'cost']);
t_is(r.reserves.qty, mpc.reserves.qty, 12, [t 'qty']);
t_is(r.reserves.totalcost, 177.8047, 4, [t 'totalcost']);

t = 'gen 5 no reserves : ';
mpc = loadcase(casefile);
mpc.reserves.zones(:, 5) = 0;
mpc.reserves.cost(5) = [];
mpc.reserves.qty(5) = [];
r = runopf_w_res(mpc, mpopt);
t_is(r.reserves.R, [25; 15; 0; 0; 0; 20], 4, [t 'R']);
t_is(r.reserves.prc, [2; 2; 2; 2; 0; 5.5], 6, [t 'prc']);
t_is(r.reserves.mu.l, [0; 0; 1; 2; 0; 0], 7, [t 'mu.l']);
t_is(r.reserves.mu.u, [0.1; 0; 0; 0; 0; 0], 6, [t 'mu.u']);
t_is(r.reserves.mu.Pmax, [0; 0; 0; 0; 0; 0], 7, [t 'mu.Pmax']);
t_is(r.reserves.cost, mpc.reserves.cost, 12, [t 'cost']);
t_is(r.reserves.qty, mpc.reserves.qty, 12, [t 'qty']);
t_is(r.reserves.totalcost, 187.5, 4, [t 'totalcost']);

t = 'extra offline gen : ';
mpc = loadcase(casefile);
idx = [1:3 5 4:6];
mpc.gen = mpc.gen(idx, :);
mpc.gencost = mpc.gencost(idx, :);
mpc.reserves.zones = mpc.reserves.zones(:, idx);
mpc.reserves.cost = mpc.reserves.cost(idx);
mpc.reserves.qty = mpc.reserves.qty(idx);
mpc.gen(4, GEN_STATUS) = 0;
r = runopf_w_res(mpc, mpopt);
t_is(r.reserves.R, [25; 15; 0; 0; 0; 19.3906; 0.6094], 4, [t 'R']);
t_is(r.reserves.prc, [2; 2; 2; 5.5; 2; 5.5; 5.5], 6, [t 'prc']);
t_is(r.reserves.mu.l, [0; 0; 1; 0; 2; 0; 0], 7, [t 'mu.l']);
t_is(r.reserves.mu.u, [0.1; 0; 0; 0; 0; 0; 0], 7, [t 'mu.u']);
t_is(r.reserves.mu.Pmax, [0; 0; 0; 0; 0; 0.5; 0], 7, [t 'mu.Pmax']);
t_is(r.reserves.cost, mpc.reserves.cost, 12, [t 'cost']);
t_is(r.reserves.qty, mpc.reserves.qty, 12, [t 'qty']);
t_is(r.reserves.totalcost, 177.8047, 4, [t 'totalcost']);

t = 'both extra & gen 6 no res : ';
mpc = loadcase(casefile);
idx = [1:3 5 4:6];
mpc.gen = mpc.gen(idx, :);
mpc.gencost = mpc.gencost(idx, :);
mpc.reserves.zones = mpc.reserves.zones(:, idx);
mpc.reserves.cost = mpc.reserves.cost(idx);
mpc.reserves.qty = mpc.reserves.qty(idx);
mpc.gen(4, GEN_STATUS) = 0;
mpc.reserves.zones(:, 6) = 0;
mpc.reserves.cost(6) = [];
mpc.reserves.qty(6) = [];
r = runopf_w_res(mpc, mpopt);
t_is(r.reserves.R, [25; 15; 0; 0; 0; 0; 20], 4, [t 'R']);
t_is(r.reserves.prc, [2; 2; 2; 5.5; 2; 0; 5.5], 6, [t 'prc']);
t_is(r.reserves.mu.l, [0; 0; 1; 0; 2; 0; 0], 7, [t 'mu.l']);
t_is(r.reserves.mu.u, [0.1; 0; 0; 0; 0; 0; 0], 6, [t 'mu.u']);
t_is(r.reserves.mu.Pmax, [0; 0; 0; 0; 0; 0; 0], 7, [t 'mu.Pmax']);
t_is(r.reserves.cost, mpc.reserves.cost, 12, [t 'cost']);
t_is(r.reserves.qty, mpc.reserves.qty, 12, [t 'qty']);
t_is(r.reserves.totalcost, 187.5, 4, [t 'totalcost']);

t = 'no qty (Rmax) : ';
mpc = loadcase(casefile);
mpc.reserves = rmfield(mpc.reserves, 'qty');
r = runopf_w_res(mpc, mpopt);
t_is(r.reserves.R, [39.3826; 0.6174; 0; 0; 19.3818; 0.6182], 4, [t 'R']);
t_is(r.reserves.prc, [2; 2; 2; 2; 5.5; 5.5], 5, [t 'prc']);
t_is(r.reserves.mu.l, [0; 0; 1; 2; 0; 0], 5, [t 'mu.l']);
t_is(r.reserves.mu.u, [0; 0; 0; 0; 0; 0], 7, [t 'mu.u']);
t_is(r.reserves.mu.Pmax, [0.1; 0; 0; 0; 0.5; 0], 5, [t 'mu.Pmax']);
t_is(r.reserves.cost, mpc.reserves.cost, 12, [t 'cost']);
t_is(r.reserves.totalcost, 176.3708, 4, [t 'totalcost']);

t = 'RAMP_10, no qty (Rmax) : ';
mpc = loadcase(casefile);
mpc.reserves = rmfield(mpc.reserves, 'qty');
mpc.gen(1, RAMP_10) = 25;
r = runopf_w_res(mpc, mpopt);
t_is(r.reserves.R, [25; 15; 0; 0; 19.3906; 0.6094], 4, [t 'R']);
t_is(r.reserves.prc, [2; 2; 2; 2; 5.5; 5.5], 6, [t 'prc']);
t_is(r.reserves.mu.l, [0; 0; 1; 2; 0; 0], 7, [t 'mu.l']);
t_is(r.reserves.mu.u, [0.1; 0; 0; 0; 0; 0], 7, [t 'mu.u']);
t_is(r.reserves.mu.Pmax, [0; 0; 0; 0; 0.5; 0], 7, [t 'mu.Pmax']);
t_is(r.reserves.cost, mpc.reserves.cost, 12, [t 'cost']);
t_is(r.reserves.totalcost, 177.8047, 4, [t 'totalcost']);

warning(s6.state, 'MATLAB:nearlySingularMatrixUMFPACK');

t_end;
