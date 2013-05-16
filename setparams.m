function pmatrix = setparams(afilename)
% sets global variables to values of all parameters used in previous run
% which generated 'afilename' paramsfile.
% returns pmatrix, a size (nparams x 3) cell array with 
% [name(str), value(num), units(str)] in each parameter row.

global r_ h_ sigma_ k_ c b beta_ eps_ mu_ dh_ K_ nrandon ;
global R_ chi_ capon Qstep hsaton Pdim1 Ldim1 x0 Gamma_ Nstep delta_ ;

fid = fopen(afilename,'r');
C = textscan(fid,'%s %n %s','delimiter',',');
nparams = size(C{1,1},1);
pmatrix = cell(nparams,3);
for i=1:nparams
    for j=1:3
        pmatrix{i,j} = C{1,j}(i,:);
    end
end

% setting global parameter values
r_ = pmatrix{1,2};
h_ = pmatrix{2,2};
sigma_ = pmatrix{3,2};
k_ = pmatrix{4,2};
c = pmatrix{5,2};
b = pmatrix{6,2};
beta_ = pmatrix{7,2};
eps_ = pmatrix{8,2};
mu_ = pmatrix{9,2};
dh_ = pmatrix{10,2};
K_ = pmatrix{11,2};
R_ = pmatrix{12,2};
capon = pmatrix{13,2};
hsaton = pmatrix{14,2};
Pdim1 = pmatrix{15,2};
Ldim1 = pmatrix{16,2};
x0 = pmatrix{17,2};
chi_ = pmatrix{18,2};
Qstep = pmatrix{19,2};
Gamma_ = pmatrix{20,2};
Nstep = pmatrix{21,2};
nrandon = pmatrix{22,2};
delta_ = pmatrix{23,2};

fclose(fid);

end