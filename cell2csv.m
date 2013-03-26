function cell2csv(filename,cellArray,toappend,delimiter)
% Writes cell array content into a file.
% 
% CELL2CSV(filename,cellArray,toappend,delimiter)
%
% filename     = Name of the file to save. [ i.e. 'text.csv' ]
% cellarray   = Name of the Cell Array with data
% toappend   = integer (1 if appending, 0 if overwriting) 
% delimiter = seperating sign, normally:',' (default)
%
% by Sylvain Fiedler, KA, 2004
% modified by Rob Kohr, Rutgers, 2005 - changed to english and fixed delimiter
if nargin<4
    delimiter = ',';
end

permission = 'a';
if ~toappend
    permission = 'w';
end

fid = fopen(filename,permission);
for i=1:size(cellArray,1)
    for j=1:size(cellArray,2)

        var = eval('cellArray{i,j}');

        if size(var,1) == 0
            var = '';
        end

        if isnumeric(var) == 1
            var = num2str(var);
        end

        fprintf(fid,'%s',var);
        
        if j ~= size(cellArray,2)
            fprintf(fid,'%s',delimiter);
        end
    end
    fprintf(fid,'\n');
end
fclose(fid);

end