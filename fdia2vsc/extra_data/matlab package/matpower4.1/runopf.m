function [MVAbase, bus, gen, gencost, branch, f, success, et] = ...
                runopf(casedata, mpopt, fname, solvedcase)
%RUNOPF  Runs an optimal power flow.
%   [RESULTS, SUCCESS] = RUNOPF(CASEDATA, MPOPT, FNAME, SOLVEDCASE)
%
%   Runs an optimal power flow (AC OPF by default), optionally returning
%   a RESULTS struct and SUCCESS flag.
%
%   Inputs (all are optional):
%       CASEDATA : either a MATPOWER case struct or a string containing
%           the name of the file with the case data (default is 'case9')
%           (see also CASEFORMAT and LOADCASE)
%       MPOPT : MATPOWER options vector to override default options
%           can be used to specify the solution algorithm, output options
%           termination tolerances, and more (see also MPOPTION).
%       FNAME : name of a file to which the pretty-printed output will
%           be appended
%       SOLVEDCASE : name of file to which the solved case will be saved
%           in MATPOWER case format (M-file will be assumed unless the
%           specified name ends with '.mat')
%
%   Outputs (all are optional):
%       RESULTS : results struct, with the following fields:
%           (all fields from the input MATPOWER case, i.e. bus, branch,
%               gen, etc., but with solved voltages, power flows, etc.)
%           order - info used in external <-> internal data conversion
%           et - elapsed time in seconds
%           success - success flag, 1 = succeeded, 0 = failed
%           (additional OPF fields, see OPF for details)
%       SUCCESS : the success flag can additionally be returned as
%           a second output argument
%
%   Calling syntax options:
%       results = runopf;
%       results = runopf(casedata);
%       results = runopf(casedata, mpopt);
%       results = runopf(casedata, mpopt, fname);
%       results = runopf(casedata, mpopt, fname, solvedcase);
%       [results, success] = runopf(...);
%
%       Alternatively, for compatibility with previous versions of MATPOWER,
%       some of the results can be returned as individual output arguments:
%
%       [baseMVA, bus, gen, gencost, branch, f, success, et] = runopf(...);
%
%   Example:
%       results = runopf('case30');
%
%   See also RUNDCOPF, RUNUOPF.

%   MATPOWER
%   $Id: runopf.m,v 1.17 2010/04/26 19:45:26 ray Exp $
%   by Ray Zimmerman, PSERC Cornell
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

%%-----  initialize  -----
%% default arguments
if nargin < 4
    solvedcase = '';                %% don't save solved case
    if nargin < 3
        fname = '';                 %% don't print results to a file
        if nargin < 2
            mpopt = mpoption;       %% use default options
            if nargin < 1
                casedata = 'case9'; %% default data file is 'case9.m'
            end
        end
    end
end

%%-----  run the optimal power flow  -----
[r, success] = opf(casedata, mpopt);

%%-----  output results  -----
if fname
    [fd, msg] = fopen(fname, 'at');
    if fd == -1
        error(msg);
    else
        printpf(r, fd, mpopt);
        fclose(fd);
    end
end
printpf(r, 1, mpopt);

%% save solved case
if solvedcase
    savecase(solvedcase, r);
end

if nargout == 1 || nargout == 2
    MVAbase = r;
    bus = success;
elseif nargout > 2
    [MVAbase, bus, gen, gencost, branch, f, et] = ...
        deal(r.baseMVA, r.bus, r.gen, r.gencost, r.branch, r.f, r.et);
% else  %% don't define MVAbase, so it doesn't print anything
end
