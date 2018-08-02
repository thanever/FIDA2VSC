function mpc = case_mc


%% bus measurement configuration
% 0 - no measurement, 1 - with measurement, 2 - with vitural measurement(also with very small sigma)
%	bus_i Mv Sigma_Mv Mp Sigma_Mp Mq Sigma_Mq
mpc.bus_mc = [
	 1	  1     1e-3  1    1e-3    1  1e-3;
	 2	  1     1e-3  1    1e-3    1  1e-3;
	 3	  1     1e-3  1    1e-3    1  1e-3;
	 4	  1     1e-3  1    1e-3    1  1e-3;
	 5	  1     1e-3  1    1e-3    1  1e-3;
	 6	  1     1e-3  1    1e-3    1  1e-3;
	 7	  1     1e-3  2    1e-3    2  1e-3;
	 8	  1     1e-3  1    1e-3    1  1e-3;
	 9	  1     1e-3  1    1e-3    1  1e-3;
	 10	  1     1e-3  1    1e-3    1  1e-3;
	 11	  1     1e-3  1    1e-3    1  1e-3;
	 12	  1     1e-3  1    1e-3    1  1e-3;
	 13	  1     1e-3  1    1e-3    1  1e-3;
	 14	  1     1e-3  1    1e-3    1  1e-3;
];

%% branch measurement configuration
%	fbus tbus Mfp, Sigma_Mfp, Mfq, Sigma_Mfq, Mtp, Sigma_Mtp, Mtq, Sigma_Mtq
mpc.branch_mc = [
	1	2	    1     1e-3     1      1e-3     1      1e-3     1     1e-3;
	1	5	    1     1e-3     1      1e-3     1      1e-3     1     1e-3;
	2	3	    1     1e-3     1      1e-3     1      1e-3     1     1e-3;
	2	4	    1     1e-3     1      1e-3     1      1e-3     1     1e-3;
	2	5	    1     1e-3     1      1e-3     1      1e-3     1     1e-3;
	3	4	    1     1e-3     1      1e-3     1      1e-3     1     1e-3;
	4	5	    1     1e-3     1      1e-3     1      1e-3     1     1e-3;
	4	7	    1     1e-3     1      1e-3     1      1e-3     1     1e-3;
	4	9	    1     1e-3     1      1e-3     1      1e-3     1     1e-3;
	5	6	    1     1e-3     1      1e-3     1      1e-3     1     1e-3;
	6	11	    1     1e-3     1      1e-3     1      1e-3     1     1e-3;
	6	12	    1     1e-3     1      1e-3     1      1e-3     1     1e-3;
	6	13	    1     1e-3     1      1e-3     1      1e-3     1     1e-3;
	7	8	    1     1e-3     1      1e-3     1      1e-3     1     1e-3;
	7	9	    1     1e-3     1      1e-3     1      1e-3     1     1e-3;
	9	10	    1     1e-3     1      1e-3     1      1e-3     1     1e-3;
	9	14	    1     1e-3     1      1e-3     1      1e-3     1     1e-3;
	10	11	    1     1e-3     1      1e-3     1      1e-3     1     1e-3;
	12	13	    1     1e-3     1      1e-3     1      1e-3     1     1e-3;
	13	14	    1     1e-3     1      1e-3     1      1e-3     1     1e-3;
];

%% bus data
%   busdc_i busac_i Mps Sigma_Mps Mqs Sigma_Mqs Mpc Sigma_Mpc Mqc Sigma_Mqc Mpcinj Sigma_Mpcing Mudc Sigma_Mudc Midc Sigma_Midc Mm Sigma_Mm Mtheta Sigma_Mtheta
mpc.busdc_mc = [
    1       6       1     1e-3     1     1e-3    1     1e-3    1     1e-3      2       1e-3       1     1e-3      1     1e-3     1    1e-3     1      1e-3;
    2       4       1     1e-3     1     1e-3    1     1e-3    1     1e-3      2       1e-3       1     1e-3      1     1e-3     1    1e-3     1      1e-3; 
];