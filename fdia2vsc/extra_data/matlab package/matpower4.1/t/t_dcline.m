function t_dcline(quiet)
%T_DCLINE  Tests for DC line extension in TOGGLE_DCLINE.

%   MATPOWER
%   $Id: t_dcline.m,v 1.1 2011/12/08 20:34:20 cvs Exp $
%   by Ray Zimmerman, PSERC Cornell
%   Copyright (c) 2011 by Power System Engineering Research Center (PSERC)
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

num_tests = 50;

t_begin(num_tests, quiet);

%% define named indices into data matrices
[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
    VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
    MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
    QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;
[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, ...
    TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
    ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = idx_brch;
[PW_LINEAR, POLYNOMIAL, MODEL, STARTUP, SHUTDOWN, NCOST, COST] = idx_cost;
c = idx_dcline;

casefile = 't_case9_dcline';
if quiet
    verbose = 0;
else
    verbose = 0;
end
if have_fcn('octave')
    s1 = warning('query', 'Octave:load-file-in-path');
    warning('off', 'Octave:load-file-in-path');
end

t0 = '';
mpopt = mpoption('OPF_VIOLATION', 1e-6, 'PDIPM_GRADTOL', 1e-8, ...
        'PDIPM_COMPTOL', 1e-8, 'PDIPM_COSTTOL', 1e-9);
mpopt = mpoption(mpopt, 'OPF_ALG', 560, 'OPF_ALG_DC', 200);
mpopt = mpoption(mpopt, 'OUT_ALL', 0, 'VERBOSE', verbose);

%% set up indices
ib_data     = [1:BUS_AREA BASE_KV:VMIN];
ib_voltage  = [VM VA];
ib_lam      = [LAM_P LAM_Q];
ib_mu       = [MU_VMAX MU_VMIN];
ig_data     = [GEN_BUS QMAX QMIN MBASE:APF];
ig_disp     = [PG QG VG];
ig_mu       = (MU_PMAX:MU_QMIN);
ibr_data    = (1:ANGMAX);
ibr_flow    = (PF:QT);
ibr_mu      = [MU_SF MU_ST];
ibr_angmu   = [MU_ANGMIN MU_ANGMAX];

%% load case
mpc0 = loadcase(casefile);
mpc0 = rmfield(mpc0, 'dclinecost');
mpc = mpc0;
mpc = toggle_dcline(mpc, 'on');
mpc = toggle_dcline(mpc, 'off');
ndc = size(mpc.dcline, 1);

%% run AC OPF w/o DC lines
t = [t0 'AC OPF (no DC lines) : '];
[r0, success] = runopf(mpc0, mpopt);
t_ok(success, [t 'success']);
[r,  success] = runopf(mpc, mpopt);
t_ok(success, [t 'success']);
t_is(r.f, r0.f, 8, [t 'f']);
t_is(   r.bus(:,ib_data   ),    r0.bus(:,ib_data   ), 10, [t 'bus data']);
t_is(   r.bus(:,ib_voltage),    r0.bus(:,ib_voltage),  3, [t 'bus voltage']);
t_is(   r.bus(:,ib_lam    ),    r0.bus(:,ib_lam    ),  3, [t 'bus lambda']);
t_is(   r.bus(:,ib_mu     ),    r0.bus(:,ib_mu     ),  2, [t 'bus mu']);
t_is(   r.gen(:,ig_data   ),    r0.gen(:,ig_data   ), 10, [t 'gen data']);
t_is(   r.gen(:,ig_disp   ),    r0.gen(:,ig_disp   ),  3, [t 'gen dispatch']);
t_is(   r.gen(:,ig_mu     ),    r0.gen(:,ig_mu     ),  3, [t 'gen mu']);
t_is(r.branch(:,ibr_data  ), r0.branch(:,ibr_data  ), 10, [t 'branch data']);
t_is(r.branch(:,ibr_flow  ), r0.branch(:,ibr_flow  ),  3, [t 'branch flow']);
t_is(r.branch(:,ibr_mu    ), r0.branch(:,ibr_mu    ),  2, [t 'branch mu']);

t = [t0 'AC PF (no DC lines) : '];
mpc1 = struct('baseMVA', [], 'bus', [], 'branch', [], 'gencost', [], 'dcline', []);
[mpc1.baseMVA, mpc1.bus, mpc1.gen, mpc1.branch, mpc1.gencost, mpc1.dcline] = ...
    deal(r.baseMVA, r.bus(:, 1:VMIN), r.gen(:, 1:APF), ...
        r.branch(:, 1:ANGMAX), r.gencost, r.dcline(:, 1:c.LOSS1));
mpc1.bus(:, VM) = 1;
mpc1.bus(:, VA) = 0;
[rp, success] = runpf(mpc1, mpopt);
t_ok(success, [t 'success']);
t_is(   rp.bus(:,ib_voltage),    r.bus(:,ib_voltage),  3, [t 'bus voltage']);
t_is(   rp.gen(:,ig_disp   ),    r.gen(:,ig_disp   ),  3, [t 'gen dispatch']);
t_is(rp.branch(:,ibr_flow  ), r.branch(:,ibr_flow  ),  3, [t 'branch flow']);

%% run with DC lines
t = [t0 'AC OPF (with DC lines) : '];
mpc = toggle_dcline(mpc, 'on');
[r, success] = runopf(mpc, mpopt);
t_ok(success, [t 'success']);
expected = [
	10	8.9	-10	10	1.0674	1.0935;
	2.2776	2.2776	0	0	1.0818	1.0665;
	0	0	0	0	1.0000	1.0000;
	10	9.5	0.0563	-10	1.0778	1.0665;
];
t_is(r.dcline(:, c.PF:c.VT), expected, 4, [t 'P Q V']);
expected = [
	0	0.8490	0.6165	0	0	0.2938;
	0	0	0	0.4290	0.0739	0;
	0	0	0	0	0	0;
	0	7.2209	0	0	0.0739	0;
];
t_is(r.dcline(:, c.MU_PMIN:c.MU_QMAXT), expected, 3, [t 'mu']);

t = [t0 'AC PF (with DC lines) : '];
mpc1 = struct('baseMVA', [], 'bus', [], 'branch', [], 'gencost', [], 'dcline', []);
[mpc1.baseMVA, mpc1.bus, mpc1.gen, mpc1.branch, mpc1.gencost, mpc1.dcline] = ...
    deal(r.baseMVA, r.bus(:, 1:VMIN), r.gen(:, 1:APF), ...
        r.branch(:, 1:ANGMAX), r.gencost, r.dcline(:, 1:c.LOSS1));
mpc1 = toggle_dcline(mpc1, 'on');
mpc1.bus(:, VM) = 1;
mpc1.bus(:, VA) = 0;
[rp, success] = runpf(mpc1, mpopt);
t_ok(success, [t 'success']);
t_is(   rp.bus(:,ib_voltage),    r.bus(:,ib_voltage), 3, [t 'bus voltage']);
%t_is(   rp.gen(:,ig_disp   ),    r.gen(:,ig_disp   ), 3, [t 'gen dispatch']);
t_is(   rp.gen(1:2,ig_disp ),    r.gen(1:2,ig_disp ), 3, [t 'gen dispatch']);
t_is(   rp.gen(3,PG        ),    r.gen(3,PG        ), 3, [t 'gen dispatch']);
t_is(   rp.gen(3,QG)+rp.dcline(1,c.QF), r.gen(3,QG)+r.dcline(1,c.QF), 3, [t 'gen dispatch']);
t_is(rp.branch(:,ibr_flow  ), r.branch(:,ibr_flow  ), 3, [t 'branch flow']);

%% add appropriate P and Q injections and check angles and generation when running PF
t = [t0 'AC PF (with equivalent injections) : '];
mpc1 = struct('baseMVA', [], 'bus', [], 'branch', [], 'gencost', [], 'dcline', []);
[mpc1.baseMVA, mpc1.bus, mpc1.gen, mpc1.branch, mpc1.gencost, mpc1.dcline] = ...
    deal(r.baseMVA, r.bus(:, 1:VMIN), r.gen(:, 1:APF), ...
        r.branch(:, 1:ANGMAX), r.gencost, r.dcline(:, 1:c.LOSS1));
mpc1.bus(:, VM) = 1;
mpc1.bus(:, VA) = 0;
for k = 1:ndc
    if mpc1.dcline(k, c.BR_STATUS)
        ff = find(mpc1.bus(:, BUS_I) == mpc1.dcline(k, c.F_BUS));
        tt = find(mpc1.bus(:, BUS_I) == mpc1.dcline(k, c.T_BUS));
        mpc1.bus(ff, PD) = mpc1.bus(ff, PD) + r.dcline(k, c.PF);
        mpc1.bus(ff, QD) = mpc1.bus(ff, QD) - r.dcline(k, c.QF);
        mpc1.bus(tt, PD) = mpc1.bus(tt, PD) - r.dcline(k, c.PT);
        mpc1.bus(tt, QD) = mpc1.bus(tt, QD) - r.dcline(k, c.QT);
        mpc1.bus(ff, VM) = r.dcline(k, c.VF);
        mpc1.bus(tt, VM) = r.dcline(k, c.VT);
        mpc1.bus(ff, BUS_TYPE) = PV;
        mpc1.bus(tt, BUS_TYPE) = PV;
    end
end
[rp, success] = runpf(mpc1, mpopt);
t_ok(success, [t 'success']);
t_is(   rp.bus(:,ib_voltage),    r.bus(:,ib_voltage),  3, [t 'bus voltage']);
t_is(   rp.gen(:,ig_disp   ),    r.gen(:,ig_disp   ),  3, [t 'gen dispatch']);
t_is(rp.branch(:,ibr_flow  ), r.branch(:,ibr_flow  ),  3, [t 'branch flow']);

%% test DC OPF
t = [t0 'DC OPF (with DC lines) : '];
mpc = mpc0;
mpc.gen(1, PMIN) = 10;
mpc.branch(5, RATE_A) = 100;
mpc = toggle_dcline(mpc, 'on');
[r, success] = rundcopf(mpc, mpopt);
t_ok(success, [t 'success']);
expected = [
	10	8.9	0	0	1.01	1;
	2	2	0	0	1	1;
	0	0	0	0	1	1;
	10	9.5	0	0	1	0.98;
];
t_is(r.dcline(:, c.PF:c.VT), expected, 4, [t 'P Q V']);
expected = [
	0	1.8602	0	0	0	0;
	1.8507	0	0	0	0	0;
	0	0	0	0	0	0;
	0	0.2681	0	0	0	0;
];
t_is(r.dcline(:, c.MU_PMIN:c.MU_QMAXT), expected, 3, [t 'mu']);

t = [t0 'DC PF (with DC lines) : '];
mpc1 = struct('baseMVA', [], 'bus', [], 'branch', [], 'gencost', [], 'dcline', []);
[mpc1.baseMVA, mpc1.bus, mpc1.gen, mpc1.branch, mpc1.gencost, mpc1.dcline] = ...
    deal(r.baseMVA, r.bus(:, 1:VMIN), r.gen(:, 1:APF), ...
        r.branch(:, 1:ANGMAX), r.gencost, r.dcline(:, 1:c.LOSS1));
mpc1 = toggle_dcline(mpc1, 'on');
mpc1.bus(:, VA) = 0;
[rp, success] = rundcpf(mpc1, mpopt);
t_ok(success, [t 'success']);
t_is(   rp.bus(:,ib_voltage),    r.bus(:,ib_voltage), 3, [t 'bus voltage']);
t_is(   rp.gen(:,ig_disp   ),    r.gen(:,ig_disp   ), 3, [t 'gen dispatch']);
t_is(rp.branch(:,ibr_flow  ), r.branch(:,ibr_flow  ), 3, [t 'branch flow']);

%% add appropriate P injections and check angles and generation when running PF
t = [t0 'DC PF (with equivalent injections) : '];
mpc1 = struct('baseMVA', [], 'bus', [], 'branch', [], 'gencost', [], 'dcline', []);
[mpc1.baseMVA, mpc1.bus, mpc1.gen, mpc1.branch, mpc1.gencost, mpc1.dcline] = ...
    deal(r.baseMVA, r.bus(:, 1:VMIN), r.gen(:, 1:APF), ...
        r.branch(:, 1:ANGMAX), r.gencost, r.dcline(:, 1:c.LOSS1));
mpc1.bus(:, VA) = 0;
for k = 1:ndc
    if mpc1.dcline(k, c.BR_STATUS)
        ff = find(mpc1.bus(:, BUS_I) == mpc1.dcline(k, c.F_BUS));
        tt = find(mpc1.bus(:, BUS_I) == mpc1.dcline(k, c.T_BUS));
        mpc1.bus(ff, PD) = mpc1.bus(ff, PD) + r.dcline(k, c.PF);
        mpc1.bus(tt, PD) = mpc1.bus(tt, PD) - r.dcline(k, c.PT);
        mpc1.bus(ff, BUS_TYPE) = PV;
        mpc1.bus(tt, BUS_TYPE) = PV;
    end
end
[rp, success] = rundcpf(mpc1, mpopt);
t_ok(success, [t 'success']);
t_is(   rp.bus(:,ib_voltage),    r.bus(:,ib_voltage),  3, [t 'bus voltage']);
t_is(   rp.gen(:,ig_disp   ),    r.gen(:,ig_disp   ),  3, [t 'gen dispatch']);
t_is(rp.branch(:,ibr_flow  ), r.branch(:,ibr_flow  ),  3, [t 'branch flow']);

%% run with DC lines
t = [t0 'AC OPF (with DC lines + poly cost) : '];
mpc = loadcase(casefile);
mpc = toggle_dcline(mpc, 'on');
[r, success] = runopf(mpc, mpopt);
t_ok(success, [t 'success']);
expected1 = [
	10	8.9	-10	10	1.0663	1.0936;
	7.8429	7.8429	0	0	1.0809	1.0667;
	0	0	0	0	1.0000	1.0000;
	6.0549	5.7522	-0.5897	-10	1.0778	1.0667;
];
t_is(r.dcline(:, c.PF:c.VT), expected1, 4, [t 'P Q V']);
expected2 = [
	0	0.7605	0.6226	0	0	0.2980;
	0	0	0	0.4275	0.0792	0;
	0	0	0	0	0	0;
	0	0	0	0	0.0792	0;
];
t_is(r.dcline(:, c.MU_PMIN:c.MU_QMAXT), expected2, 3, [t 'mu']);

mpc.dclinecost(4, 1:8) = [2 0 0 4 0 0 7.3 0];
[r, success] = runopf(mpc, mpopt);
t_ok(success, [t 'success']);
t_is(r.dcline(:, c.PF:c.VT), expected1, 4, [t 'P Q V']);
t_is(r.dcline(:, c.MU_PMIN:c.MU_QMAXT), expected2, 3, [t 'mu']);

t = [t0 'AC OPF (with DC lines + pwl cost) : '];
mpc.dclinecost(4, 1:8) = [1 0 0 2 0 0 10 73];
[r, success] = runopf(mpc, mpopt);
t_ok(success, [t 'success']);
t_is(r.dcline(:, c.PF:c.VT), expected1, 4, [t 'P Q V']);
t_is(r.dcline(:, c.MU_PMIN:c.MU_QMAXT), expected2, 3, [t 'mu']);

if have_fcn('octave')
    warning(s1.state, 'Octave:load-file-in-path');
end

t_end;
