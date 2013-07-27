%% markov_pi

N = 10000; % number of steps
d = 0.3;   % throwing range
clubhouse = [1 1];  % initial location

Nxy = [clubhouse;zeros(N,2)];
deltaxy = -d + (2*d).*rand(N,2);
tic
for i=1:N
    Nxy(i+1,:) = Nxy(i,:)+deltaxy(i,:).*(abs(Nxy(i,:)+deltaxy(i,:))<1);
end
Nhits = sum(sum(Nxy.^2,2)<1)-(sum(clubhouse.^2)>1);
toc
disp(Nhits*4/N); % estimate of pi

%% markov_pi_periodic

N = 10000; % number of steps
d = 0.3;   % throwing range
clubhouse = [1 1];  % initial location

Nxy = [clubhouse;zeros(N,2)];
deltaxy = -d + (2*d).*rand(N,2);
tic
for i=1:N
    Nxy(i+1,:) = Nxy(i,:)+ deltaxy(i,:)...
            - 2.*((Nxy(i,:)+deltaxy(i,:))>1)...
            + 2.*((Nxy(i,:)+deltaxy(i,:))<-1);
end
Nhits = sum(sum(Nxy.^2,2)<1)-(sum(clubhouse.^2)>1);
toc
disp(Nhits*4/N); % estimate of pi