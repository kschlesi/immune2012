% concatenate two continuous saved runs
function concatenate = concat(run1,run2)

Plin1 = csvread(run1);
Plin2in = csvread(run2);
Plin2 = Plin2in(2:end,:);
Plin3 = [Plin1;Plin2];
dlmwrite(run1,Plin3);
delete(run2);

concatenate=11;
