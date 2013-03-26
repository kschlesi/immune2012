function writeparams(afilename,a0)

paramnames = {'r_';'h_';'sigma_';'de_';'k_';'f_';'c';'b';'beta_';'eps_';'mu_';...
    'dh_';'K_';'R_';'capon';'hsaton';'Pdim1';'Ldim1';'x0'};
paramunits = {'day^(-1)';'\mul/cell/day';'day^(-1)';'day^(-1)';'cell/\mul';'--';'--';...
    'sites';'sites';'sites';'cell/site';'day^(-1)';'cell/\mul';'cell/\mul';'--';'--';...
    'sites';'sites';'--'};
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