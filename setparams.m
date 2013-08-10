function pmatrix = setparams(afilename)
% sets global variables to values of all parameters used in previous run
% which generated 'afilename' paramsfile.
% returns pmatrix, a size (nparams x 3) cell array with 
% [name(str), value(num), units(str)] in each parameter row.

fid = fopen(afilename,'r');
C = textscan(fid,'%s %n %s','delimiter',',');
nparams = size(C{1,1},1);
pmatrix = cell(nparams,3);
for i=1:nparams
    for j=1:3
        pmatrix{i,j} = C{1,j}(i,:);
    end
end

fclose(fid);

end