function fname_out = savecase(fname, varargin)
%SAVECASE  Saves a MATPOWER case file, given a filename and the data.
%   SAVECASE(FNAME, CASESTRUCT)
%   SAVECASE(FNAME, CASESTRUCT, VERSION)
%   SAVECASE(FNAME, BASEMVA, BUS, GEN, BRANCH)
%   SAVECASE(FNAME, BASEMVA, BUS, GEN, BRANCH, GENCOST)
%   SAVECASE(FNAME, BASEMVA, BUS, GEN, BRANCH, AREAS, GENCOST)
%   SAVECASE(FNAME, COMMENT, CASESTRUCT)
%   SAVECASE(FNAME, COMMENT, CASESTRUCT, VERSION)
%   SAVECASE(FNAME, COMMENT, BASEMVA, BUS, GEN, BRANCH)
%   SAVECASE(FNAME, COMMENT, BASEMVA, BUS, GEN, BRANCH, GENCOST)
%   SAVECASE(FNAME, COMMENT, BASEMVA, BUS, GEN, BRANCH, AREAS, GENCOST)
%
%   FNAME = SAVECASE(FNAME, ...)
%
%   Writes a MATPOWER case file, given a filename and data struct or list of
%   data matrices. The FNAME parameter is the name of the file to be created or
%   overwritten. If FNAME ends with '.mat' it saves the case as a MAT-file
%   otherwise it saves it as an M-file. Optionally returns the filename,
%   with extension added if necessary. The optional COMMENT argument is
%   either string (single line comment) or a cell array of strings which
%   are inserted as comments. When using a MATPOWER case struct, if the
%   optional VERSION argument is '1' it will modify the data matrices to
%   version 1 format before saving.

%   MATPOWER
%   $Id: savecase.m,v 1.36 2011/09/28 20:10:23 cvs Exp $
%   by Carlos E. Murillo-Sanchez, PSERC Cornell & Universidad Autonoma de Manizales
%   and Ray Zimmerman, PSERC Cornell
%   Copyright (c) 1996-2010 by Power System Engineering Research Center (PSERC)
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

%% define named indices into bus, gen, branch matrices
[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
    VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
    MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
    QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;
[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, ...
    TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
    ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = idx_brch;
[AREA_I, PRICE_REF_BUS] = idx_area;
[PW_LINEAR, POLYNOMIAL, MODEL, STARTUP, SHUTDOWN, NCOST, COST] = idx_cost;

%% default arguments
if ischar(varargin{1}) || iscell(varargin{1})
    comment = varargin{1};
    [args{1:(length(varargin)-1)}] = deal(varargin{2:end});
else
    comment = '';
    args = varargin;
end
mpc_ver = '2';              %% default MATPOWER case file version
if isstruct(args{1})        %% 1st real argument is a struct
    mpc = args{1};
    if length(args) > 1
        mpc.version = args{2};
        mpc_ver = mpc.version;
    end
    baseMVA = mpc.baseMVA;
    bus     = mpc.bus;
    gen     = mpc.gen;
    branch  = mpc.branch;
    if isfield(mpc, 'areas')
        areas   = mpc.areas;
    end
    if isfield(mpc, 'gencost')
        gencost = mpc.gencost;
    end
else                        %% 1st real argument is NOT a struct
    baseMVA = args{1};
    bus     = args{2};
    gen     = args{3};
    branch  = args{4};
    mpc.baseMVA = baseMVA;
    mpc.bus     = bus;
    mpc.gen     = gen;
    mpc.branch  = branch;
    if length(args) == 5
        gencost = args{5};
        mpc.gencost = gencost;
    end
    if length(args) == 6
        areas   = args{5};
        gencost = args{6};
        mpc.areas   = areas;
        mpc.gencost = gencost;
    end
end

%% modifications for version 1 format
if strcmp(mpc_ver, '1')
    %% remove extra columns of gen
    if size(gen, 2) >= MU_QMIN
        gen = gen(:, [1:PMIN, MU_PMAX:MU_QMIN]);
    else
        gen = gen(:, 1:PMIN);
    end
    %% use the version 1 values for column names
    shift = MU_PMAX - PMIN - 1;
    tmp = num2cell([MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN] - shift);
    [MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN] = deal(tmp{:});
    
    %% remove extra columns of branch
    if size(branch, 2) >= MU_ST
        branch = branch(:, [1:BR_STATUS, PF:MU_ST]);
    elseif size(branch, 2) >= QT
        branch = branch(:, [1:BR_STATUS, PF:QT]);
    else
        branch = branch(:, 1:BR_STATUS);
    end
    %% use the version 1 values for column names
    shift = PF - BR_STATUS - 1;
    tmp = num2cell([PF, QF, PT, QT, MU_SF, MU_ST] - shift);
    [PF, QF, PT, QT, MU_SF, MU_ST] = deal(tmp{:});
end

%% verify valid filename
l = length(fname);
rootname = [];
if l > 2
    if strcmp(fname(l-1:l), '.m')
        rootname = fname(1:l-2);
        extension = '.m';
    elseif l > 4
        if strcmp(fname(l-3:l), '.mat')
            rootname = fname(1:l-4);
            extension = '.mat';
        end
    end
end
if isempty(rootname)
    rootname = fname;
    extension = '.m';
    fname = [rootname, extension];
end

%% open and write the file
if strcmp(extension, '.mat')        %% MAT-file
    vflag = '';
    if str2double(version('-release')) > 13
        vflag = ' -V6';
    end
    if strcmp(mpc_ver, '1')
        if exist('areas', 'var') && exist('gencost', 'var')
            cmd = sprintf('save %s baseMVA bus gen branch areas gencost%s;', rootname, vflag); 
        else
            cmd = sprintf('save %s baseMVA bus gen branch%s;', rootname, vflag); 
        end
    else
        cmd = sprintf('save %s mpc%s;', rootname, vflag); 
    end
    eval(cmd);
else                                %% M-file
    %% open file
    [fd, msg] = fopen(fname, 'wt');     %% print it to an M-file
    if fd == -1
        error(['savecase: ', msg]);
    end
    
    %% function header, etc.
    if strcmp(mpc_ver, '1')
        if exist('areas', 'var') && exist('gencost', 'var') && ~isempty(gencost)
            fprintf(fd, 'function [baseMVA, bus, gen, branch, areas, gencost] = %s\n', rootname);
        else
            fprintf(fd, 'function [baseMVA, bus, gen, branch] = %s\n', rootname);
        end
        prefix = '';
    else
        fprintf(fd, 'function mpc = %s\n', rootname);
        prefix = 'mpc.';
    end
    if ~isempty(comment)
        if ischar(comment)
            fprintf(fd, '%%%s\n', comment);
        elseif iscell(comment)
            for k = 1:length(comment)
                fprintf(fd, '%%%s\n', comment{k});
            end
        end
    end
    fprintf(fd, '\n%%%% MATPOWER Case Format : Version %s\n', mpc_ver);    
    if ~strcmp(mpc_ver, '1')
        fprintf(fd, 'mpc.version = ''%s'';\n', mpc_ver);    
    end
    fprintf(fd, '\n%%%%-----  Power Flow Data  -----%%%%\n');
    fprintf(fd, '%%%% system MVA base\n');
    fprintf(fd, '%sbaseMVA = %.9g;\n', prefix, baseMVA);
    
    %% bus data
    ncols = size(bus, 2);
    fprintf(fd, '\n%%%% bus data\n');
    fprintf(fd, '%%\tbus_i\ttype\tPd\tQd\tGs\tBs\tarea\tVm\tVa\tbaseKV\tzone\tVmax\tVmin');
    if ncols >= MU_VMIN             %% opf SOLVED, save with lambda's & mu's
        fprintf(fd, '\tlam_P\tlam_Q\tmu_Vmax\tmu_Vmin');
    end
    fprintf(fd, '\n%sbus = [\n', prefix);
    if ncols < MU_VMIN              %% opf NOT SOLVED, save without lambda's & mu's
        fprintf(fd, '\t%d\t%d\t%.9g\t%.9g\t%.9g\t%.9g\t%d\t%.9g\t%.9g\t%.9g\t%d\t%.9g\t%.9g;\n', bus(:, 1:VMIN).');
    else                            %% opf SOLVED, save with lambda's & mu's
        fprintf(fd, '\t%d\t%d\t%.9g\t%.9g\t%.9g\t%.9g\t%d\t%.9g\t%.9g\t%.9g\t%d\t%.9g\t%.9g\t%.4f\t%.4f\t%.4f\t%.4f;\n', bus(:, 1:MU_VMIN).');
    end
    fprintf(fd, '];\n');
    
    %% generator data
    ncols = size(gen, 2);
    fprintf(fd, '\n%%%% generator data\n');
    fprintf(fd, '%%\tbus\tPg\tQg\tQmax\tQmin\tVg\tmBase\tstatus\tPmax\tPmin');
    if ~strcmp(mpc_ver, '1')
        fprintf(fd, '\tPc1\tPc2\tQc1min\tQc1max\tQc2min\tQc2max\tramp_agc\tramp_10\tramp_30\tramp_q\tapf');
    end
    if ncols >= MU_QMIN             %% opf SOLVED, save with mu's
        fprintf(fd, '\tmu_Pmax\tmu_Pmin\tmu_Qmax\tmu_Qmin');
    end
    fprintf(fd, '\n%sgen = [\n', prefix);
    if ncols < MU_QMIN              %% opf NOT SOLVED, save without mu's
        if strcmp(mpc_ver, '1')
            fprintf(fd, '\t%d\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%d\t%.9g\t%.9g;\n', gen(:, 1:PMIN).');
        else
            fprintf(fd, '\t%d\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%d\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g;\n', gen(:, 1:APF).');
        end
    else
        if strcmp(mpc_ver, '1')
            fprintf(fd, '\t%d\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%d\t%.9g\t%.9g\t%.4f\t%.4f\t%.4f\t%.4f;\n', gen(:, 1:MU_QMIN).');
        else
            fprintf(fd, '\t%d\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%d\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%.4f\t%.4f\t%.4f\t%.4f;\n', gen(:, 1:MU_QMIN).');
        end
    end
    fprintf(fd, '];\n');
    
    %% branch data
    ncols = size(branch, 2);
    fprintf(fd, '\n%%%% branch data\n');
    fprintf(fd, '%%\tfbus\ttbus\tr\tx\tb\trateA\trateB\trateC\tratio\tangle\tstatus');
    if ~strcmp(mpc_ver, '1')
        fprintf(fd, '\tangmin\tangmax');
    end
    if ncols >= QT                  %% power flow SOLVED, save with line flows
        fprintf(fd, '\tPf\tQf\tPt\tQt');
    end
    if ncols >= MU_ST               %% opf SOLVED, save with mu's
        fprintf(fd, '\tmu_Sf\tmu_St');
        if ~strcmp(mpc_ver, '1')
            fprintf(fd, '\tmu_angmin\tmu_angmax');
        end
    end
    fprintf(fd, '\n%sbranch = [\n', prefix);
    if ncols < QT                   %% power flow NOT SOLVED, save without line flows or mu's
        if strcmp(mpc_ver, '1')
            fprintf(fd, '\t%d\t%d\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%d;\n', branch(:, 1:BR_STATUS).');
        else
            fprintf(fd, '\t%d\t%d\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%d\t%.9g\t%.9g;\n', branch(:, 1:ANGMAX).');
        end
    elseif ncols < MU_ST            %% power flow SOLVED, save with line flows but without mu's
        if strcmp(mpc_ver, '1')
            fprintf(fd, '\t%d\t%d\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%d\t%.4f\t%.4f\t%.4f\t%.4f;\n', branch(:, 1:QT).');
        else
            fprintf(fd, '\t%d\t%d\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%d\t%.9g\t%.9g\t%.4f\t%.4f\t%.4f\t%.4f;\n', branch(:, 1:QT).');
       end
    else                            %% opf SOLVED, save with lineflows & mu's
        if strcmp(mpc_ver, '1')
            fprintf(fd, '\t%d\t%d\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%d\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f;\n', branch(:, 1:MU_ST).');
        else
            fprintf(fd, '\t%d\t%d\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%d\t%.9g\t%.9g\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f;\n', branch(:, 1:MU_ANGMAX).');
        end
    end
    fprintf(fd, '];\n');
    
    %% OPF data
    if (exist('areas', 'var') && ~isempty(areas)) || ...
        (exist('gencost', 'var') && ~isempty(gencost))
        fprintf(fd, '\n%%%%-----  OPF Data  -----%%%%');
    end
    if exist('areas', 'var') && ~isempty(areas)
        %% area data
        fprintf(fd, '\n%%%% area data\n');
        fprintf(fd, '%%\tarea\trefbus\n');
        fprintf(fd, '%sareas = [\n', prefix);
        if ~isempty(areas)
            fprintf(fd, '\t%d\t%d;\n', areas(:, 1:PRICE_REF_BUS).');
        end
        fprintf(fd, '];\n');
    end
    if exist('gencost', 'var') && ~isempty(gencost)
        %% generator cost data
        fprintf(fd, '\n%%%% generator cost data\n');
        fprintf(fd, '%%\t1\tstartup\tshutdown\tn\tx1\ty1\t...\txn\tyn\n');
        fprintf(fd, '%%\t2\tstartup\tshutdown\tn\tc(n-1)\t...\tc0\n');
        fprintf(fd, '%sgencost = [\n', prefix);
        if ~isempty(gencost)
            n1 = 2 * max(gencost(gencost(:, MODEL) == PW_LINEAR,  NCOST));
            n2 =     max(gencost(gencost(:, MODEL) == POLYNOMIAL, NCOST));
            n = max([n1; n2]);
            if size(gencost, 2) < n + 4
                error('savecase: gencost data claims it has more columns than it does');
            end
            template = '\t%d\t%.9g\t%.9g\t%d';
            for i = 1:n
                template = [template, '\t%.9g'];
            end
            template = [template, ';\n'];
            fprintf(fd, template, gencost.');
        end
        fprintf(fd, '];\n');
    end
    
    %% generalized OPF user data
    if (isfield(mpc, 'A') && ~isempty(mpc.A)) || ...
            (isfield(mpc, 'N') && ~isempty(mpc.N))
        fprintf(fd, '\n%%%%-----  Generalized OPF User Data  -----%%%%');
    end

    %% user constraints
    if isfield(mpc, 'A') && ~isempty(mpc.A)
        %% A
        fprintf(fd, '\n%%%% user constraints\n');
        print_sparse(fd, sprintf('%sA', prefix), mpc.A);
        if isfield(mpc, 'l') && ~isempty(mpc.l) && ...
                isfield(mpc, 'u') && ~isempty(mpc.u)
            fprintf(fd, 'lu = [\n');
            fprintf(fd, '\t%.9g\t%.9g;\n', [mpc.l mpc.u].');
            fprintf(fd, '];\n');
            fprintf(fd, '%sl = lu(:, 1);\n', prefix);
            fprintf(fd, '%su = lu(:, 2);\n\n', prefix);
        elseif isfield(mpc, 'l') && ~isempty(mpc.l)
            fprintf(fd, '%sl = [\n', prefix);
            fprintf(fd, '\t%.9g;\n', mpc.l);
            fprintf(fd, '];\n\n');
        elseif isfield(mpc, 'u') && ~isempty(mpc.u)
            fprintf(fd, '%su = [\n', prefix);
            fprintf(fd, '\t%.9g;\n', mpc.u);
            fprintf(fd, '];\n');
        end
    end

    %% user costs
    if isfield(mpc, 'N') && ~isempty(mpc.N)
        fprintf(fd, '\n%%%% user costs\n');
        print_sparse(fd, sprintf('%sN', prefix), mpc.N);
        if isfield(mpc, 'H') && ~isempty(mpc.H)
            print_sparse(fd, sprintf('%sH', prefix), mpc.H);
        end
        if isfield(mpc, 'fparm') && ~isempty(mpc.fparm)
            fprintf(fd, 'Cw_fparm = [\n');
            fprintf(fd, '\t%.9g\t%d\t%.9g\t%.9g\t%.9g;\n', [mpc.Cw mpc.fparm].');
            fprintf(fd, '];\n');
            fprintf(fd, '%sCw    = Cw_fparm(:, 1);\n', prefix);
            fprintf(fd, '%sfparm = Cw_fparm(:, 2:5);\n', prefix);
        else
            fprintf(fd, '%sCw = [\n', prefix);
            fprintf(fd, '\t%.9g;\n', mpc.Cw);
            fprintf(fd, '];\n');
        end
    end

    %% user vars
    if isfield(mpc, 'z0') || isfield(mpc, 'zl') || isfield(mpc, 'zu')
        fprintf(fd, '\n%%%% user vars\n');
    end
    if isfield(mpc, 'z0') && ~isempty(mpc.z0)
        fprintf(fd, '%sz0 = [\n', prefix);
        fprintf(fd, '\t%.9g;\n', mpc.z0);
        fprintf(fd, '];\n');
    end
    if isfield(mpc, 'zl') && ~isempty(mpc.zl)
        fprintf(fd, '%szl = [\n', prefix);
        fprintf(fd, '\t%.9g;\n', mpc.zl);
        fprintf(fd, '];\n');
    end
    if isfield(mpc, 'zu') && ~isempty(mpc.zu)
        fprintf(fd, '%szu = [\n', prefix);
        fprintf(fd, '\t%.9g;\n', mpc.zu);
        fprintf(fd, '];\n');
    end

    %% execute userfcn callbacks for 'savecase' stage
    if isfield(mpc, 'userfcn')
        run_userfcn(mpc.userfcn, 'savecase', mpc, fd, prefix);
    end

    %% close file
    if fd ~= 1
        fclose(fd);
    end
end

if nargout > 0
    fname_out = fname;
end



function print_sparse(fd, varname, A)

[i, j, s] = find(A);
[m, n] = size(A);

if isempty(s)
    fprintf(fd, '%s = sparse(%d, %d);\n', varname, m, n);
else
    fprintf(fd, 'ijs = [\n');
    if m == 1           %% i, j, s are row vectors
        fprintf(fd, '\t%d\t%d\t%.9g;\n', [i; j; s]);
    else                %% i, j, s are column vectors
        fprintf(fd, '\t%d\t%d\t%.9g;\n', [i j s].');
    end
    fprintf(fd, '];\n');
    fprintf(fd, '%s = sparse(ijs(:, 1), ijs(:, 2), ijs(:, 3), %d, %d);\n', varname, m, n);
end
