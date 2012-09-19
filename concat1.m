% concatenate two continuous saved runs
function concatenate = concat(run1,run2,run3)

Plin1 = csvread(run1);
Plin2in = csvread(run2);
Plin3in = csvread(run3);
Plin2 = Plin2in(2:end,:);
Plin3 = Plin3in(2:end,:);
Plin4 = [Plin1;Plin2;Plin3];
dlmwrite(run1,Plin4);

delete(run2);
delete(run3);

concatenate=1;
