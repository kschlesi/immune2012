% GILLESPIE ALGORITHM IMPLEMENTATION

rxnt = 5;    % duration of reaction to simulate
tinit = 0;   % time to start

dim1 = 10;                  % size of physical space
X0 = 100.*ones(dim1,1);     % initial prey population
Y0 = 100.*ones(dim1,1);     % initial predator population

rxno = 4;
alpha = 0.04;   % reaction rates
beta = 0.0005;
gamma = 0.1;
delta = 0.0005;

% Gillespie algorithm loop
t = tinit;
X = X0;
Y = Y0;
Xtot = sum(X0,1);
Ytot = sum(Y0,1);

% estimate initial savematrix size & initialise
smsize = floor(rxnt/((alpha*(Xtot/dim1)+(Ytot/dim1)*((Xtot/dim1)*(beta+delta)+gamma))));
savematrix = zeros(smsize,1+size(X0,1)+size(Y0,1));
savematrix(1,:) = [t,X',Y'];
smi = 1;

while t<(tinit+rxnt) && (sum(X,1)+sum(Y,1))
    
    % draw next reaction time
        a1 = alpha.*X;
        a2 = beta.*X.*Y;
        a3 = gamma.*Y;
        a4 = delta.*X.*Y;
        a0 = sum(a1 + a2 + a3 + a4,1);
        tnext = exprnd(1/a0);
        
    % draw next reaction site, compute omegas
        totpop = sum(X,1)+sum(Y,1);
        rxn = totpop;
        while (rxn==totpop && rxn)
            rxn = rand*totpop;
        end
        i = find(cumsum(sum([X,Y],2))>rxn,1,'first');
        %disp([rxn,i]);
        
    % draw next reaction and update
        rxn = a0;
        while (rxn==a0 && rxn)
            rxn = rand*a0;
        end
        
        if 0<=rxn && rxn<a1(i)
            X = X+1;
        end
        if a1(i)<=rxn && rxn<(a1(i)+a2(i))
            X = X-1;
        end
        if (a1(i)+a2(i))<=rxn && rxn<(a1(i)+a2(i)+a3(i))
            Y = Y-1;
        end
        if (a1(i)+a2(i)+a3(i))<=rxn && rxn<a0
            Y = Y+1;
        end
        
        t = t + tnext;
        disp(t);
        
        savematrix(smi+1,:) = [t,X',Y'];
        smi = smi+1;   % number of rows written in savematrix
        if (smi==smsize)
           savematrix = [savematrix;zeros(max(10,rxnt*10^-3),1+size(X0,1)+size(Y0,1))]; 
        end
        
end
cutrow = find(~sum(savematrix,2),1,'first');
if cutrow
    savematrix = savematrix(1:cutrow-1,:);
end
dlmwrite('/Users/kimberly/Documents/MATLAB/gill1.txt',savematrix);
result = csvread('/Users/kimberly/Documents/MATLAB/gill1.txt');
figure
plot(result(:,1),result(:,2),'r',result(:,1),result(:,3),'b');
%figure
%plot(result(:,1))
        
        
