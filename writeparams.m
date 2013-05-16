function writeparams(afilename,a0)

% if adding a parameter: - add paramname and paramunit in this script
%                        - add parameter to all ss_ global lists
%                        - define parameter in _main
%                        - add parameter to a0 vector in _main
%                        - add functionality of parameter to all scripts
%                        - add parameter to global list in setparams
%                        - set global parameter from pmatrix in setparams

paramnames = {'r_';'h_';'sigma_';'k_';'c';'b';'beta_';'eps_';'mu_';...
    'dh_';'K_';'R_';'capon';'hsaton';'Pdim1';'Ldim1';'x0';'chi_';'Qstep';...
    'Gamma_';'Nstep';'nrandon';'delta_'};
paramunits = {'day^(-1)';'\mul/cell/day';'day^(-1)';'cell/\mul';'--';...
    'sites';'sites';'sites';'cell/site';'day^(-1)';'cell/\mul';'cell/\mul';'--';'--';...
    'sites';'sites';'--';'--';'days';'cell/site/day';'days';'--';'day^{-1}'};
pmatrix = cell(size(a0,1),3); %index: row = parameter; 3 columns (name, number, units)
for i=1:size(a0,1)
    pmatrix{i,1} = paramnames{i};
    pmatrix{i,2} = num2str(a0(i));
    pmatrix{i,3} = paramunits{i};
end

if size(paramnames)~=size(a0)
    error('Wrong number of input parameters!')
end

% ensuring no overwrite
if isequal(exist(afilename,'file'),2)
    error('Cannot overwrite existing parameter file!');
end

cell2csv(afilename,pmatrix,0); % will DELETE what is currently in afilename

end