function writeparams(afilename,a0)

% if adding a parameter: - add paramname and paramunit in this script
%                        - define parameter in _main
%                        - add parameter to a0 vector in _main, _dy, and _seed
%                        - add functionality of parameter to all scripts

paramnames = {'r_';'h_';'sigma_';'k_';'b';'eps_';'mu_';...
    'dh_';'K_';'R_';'capon';'hsaton';'Pdim1';'Ldim1';'x0';'chi_';'Qstep';...
    'Gamma_';'nrandon';'delta_';'muton';'pinit';'w_'};
paramunits = {'day^(-1)';'\mul/cell/day';'day^(-1)';'cell/\mul';...
    'sites';'sites';'cell/site';'day^(-1)';'cell/\mul';'cell/\mul';'--';'--';...
    'sites';'sites';'--';'sites';'days';'cell/site/day';'--';'day^{-1}';'--';'cell/\mul';'--'};

if size(paramnames)~=size(a0)
    error('Wrong number of input parameters!')
end

pmatrix = cell(size(a0,1),3); %index: row = parameter; 3 columns (name, number, units)
for i=1:size(a0,1)
    pmatrix{i,1} = paramnames{i};
    pmatrix{i,2} = num2str(a0(i));
    pmatrix{i,3} = paramunits{i};
end

% ensuring no overwrite
if isequal(exist(afilename,'file'),2)
    error('Cannot overwrite existing parameter file!');
end

cell2csv(afilename,pmatrix,0); % will DELETE what is currently in afilename

end