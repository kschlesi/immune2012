%ss_spec: makes a series of plots at times given in "frametimes"
%to show time-slices of the speciation process

global r_ h_ sigma_ k_ c dh_ K_ R_ capon hsaton nrandon mrates Gamma_ delta_ muton beta_;
global b eps_ mu_ Pdim1 Ldim1 x0 chi_ Qstep Nstep tgone ntgone gammas1D lambdas1D;

frametimes = (100:5:150);   % row vector of times (in days) at which to plot
%frametimes = [100,110];
ismovie = 1;
runnum = 1;
basecode = 'pldyn';
Pdim1 = 400;
Ldim1 = 400;

y_min = 0;
y_max = 0;
x_min = 0;
x_max = min(Pdim1,Ldim1);

datapath = ['/Users/kimberly/Google Drive/immunedata/PL/' basecode '/']; %KONG datapath
%datapath = ['C:\Users\Kimberly\Google Drive\immunedata\PL\' basecode '\']; %laptop datapath
bfilename = [datapath 'b' basecode num2str(runnum) '.txt'];
tfilename = [datapath 't' basecode num2str(runnum) '.txt'];
Pfilename = [datapath 'P' basecode num2str(runnum) '.txt'];
Lfilename = [datapath 'L' basecode num2str(runnum) '.txt'];

params = setparams(bfilename);
days = params{end,2};

timevec = csvread(tfilename);
nframes = size(frametimes,2);
vid = moviein(nframes); % should include gcf argument??
for i=1:nframes
    if(frametimes(i)<days)
        tin = find(timevec>frametimes(i));
        t0index = tin(1)-1;
    else
        t0index = size(timevec,1);
    end
    P = transpose(csvread(Pfilename,t0index-1,0,[t0index-1,0,t0index-1,Pdim1-1]));
    L = transpose(csvread(Lfilename,t0index-1,0,[t0index-1,0,t0index-1,Ldim1-1]));
    y_max = max(y_max,max(max(P),max(L)));
    
figure;
hold on;
hold all;
plot((1:1:400),P);
plot((1:1:400),L);
title(['P and L distributions at t = ' num2str(frametimes(i)) ' days']);
xlabel('location in shape space (site)');
ylabel('population (cells/\mul)');
legend('Pathogen','Lymphocytes');
axis([x_min x_max y_min y_max]);
    if ismovie
        vid(i) = getframe;
        if i == nframes
            mhandle = get(gcf,'CurrentAxes');
        end
    end
hold off; 
end

if ismovie
    movie(mhandle,vid,1,1,[x_min y_min 0 0]);
end
    