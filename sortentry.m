function sortmat = sortentry(inmat,dir,which)
% given a 2D matrix 'inmat,' this function treats either rows or columns as
% linked entries (indicated by 'dir') and sorts the entries by the value of
% a given field 'which.' Default: rows are entries; 'which' names a column
% by whose values to sort

if strcmp(dir,'row')
    inmat = inmat';
end

[sortmat,ix] = sort(inmat,1);
sortmat = zeros(size(sortmat));
for i=1:size(inmat,1)
    sortmat(i,:) = inmat(ix(i,which),:);
end

if strcmp(dir,'row')
    sortmat = sortmat';
end

end