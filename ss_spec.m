%ss_spec: makes a series of plots at times given in "frametimes"
%to show time-slices of the speciation process

clear

global r_ h_ sigma_ k_ c dh_ K_ R_ capon hsaton nrandon mrates Gamma_ delta_ muton beta_;
global b eps_ mu_ Pdim1 Ldim1 x0 chi_ Qstep Nstep tgone ntgone gammas1D lambdas1D;

%frametimes = (380:15:440);   % row vector of times (in days) at which to plot
frametimes = 470;
%frametimes = (85:60:265);
ismovie = 0;
runnum = 1;
basecode = 'pldyn';
Pdim1 = 400;
Ldim1 = 400;

y_min = 1;
y_max = 23*10^4;
x_min = 30;
x_max = 90;
% x_min = 50;
% x_max = 70;
%x_max = min(Pdim1,Ldim1);

datapath = ['/Users/kimberly/Google Drive/immunedata/PL13/' basecode '/'];
bfilename = [datapath 'b' basecode num2str(runnum) '.txt'];
tfilename = [datapath 't' basecode num2str(runnum) '.txt'];
Pfilename = [datapath 'P' basecode num2str(runnum) '.txt'];
Lfilename = [datapath 'L' basecode num2str(runnum) '.txt'];

params = setparams(bfilename);
days = params{end,2};

timevec = csvread(tfilename);
nframes = size(frametimes,2);
vid = moviein(nframes,gcf);
for i=1:nframes
    if(frametimes(i)<days)
        tin = find(timevec>frametimes(i));
        t0index = tin(1)-1;
    else
        t0index = size(timevec,1);
    end
    P = csvread(Pfilename,t0index-1,0,[t0index-1,0,t0index-1,Pdim1-1]);
    L = csvread(Lfilename,t0index-1,0,[t0index-1,0,t0index-1,Ldim1-1]);
    %y_max = max(y_max,max(max(P),max(L)));
    disp(frametimes(i));
    disp(P(52:64));

if ~ismovie
    figure
end
semilogy((1:1:400),P,(1:1:400),L);xlabel('location in shape space (site)');
ylabel('population density (cells/\mul)');
legend('Pathogen','Lymphocytes');
axis([x_min x_max y_min y_max]);
set(gca,'NextPlot','replacechildren');
title(['P and L distributions at t = ' num2str(frametimes(i)) ' days']);
    if ismovie
        vid(i) = getframe(gcf);
    end
end

if ismovie
    movie(gcf,vid,4,4,[x_min y_min 0 0]);
end
    