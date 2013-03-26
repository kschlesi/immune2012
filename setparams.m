function pmatrix = setparams(afilename)
% sets global variables to values of all parameters used in previous run
% which generated 'afilename' paramsfile.
% returns pmatrix, a size (nparams x 3) cell array with 
% [name(str), value(num), units(str)] in each parameter row.

global r_ h_ sigma_ de_ k_ f_ c b beta_ eps_ mu_ dh_ K_ R_ capon hsaton Pdim1 Ldim1 x0;

fid = fopen(afilename,'r');
C = textscan(fid,'%s %n %s','delimiter',',');
nparams = size(C{1,1},1);
pmatrix = cell(nparams,3);
for i=1:nparams
    pmatrix{i,1} = C{1,1}(i,:);
    pmatrix{i,2} = C{1,2}(i,:);
    pmatrix{i,3} = C{1,3}(i,:);
end

%pmatrix = nparams;

% setting global parameter values
r_ = pmatrix{1,2};
h_ = pmatrix{2,2};
sigma_ = pmatrix{3,2};
de_ = pmatrix{4,2};
k_ = pmatrix{5,2};
f_ = pmatrix{6,2};
c = pmatrix{7,2};
b = pmatrix{8,2};
beta_ = pmatrix{9,2};
eps_ = pmatrix{10,2};
mu_ = pmatrix{11,2};
dh_ = pmatrix{12,2};
K_ = pmatrix{13,2};
R_ = pmatrix{14,2};
capon = pmatrix{15,2};
hsaton = pmatrix{16,2};
Pdim1 = pmatrix{17,2};
Ldim1 = pmatrix{18,2};
x0 = pmatrix{19,2};

fclose(fid)

end