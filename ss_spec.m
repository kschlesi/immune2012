%ss_spec: makes a speciation series plot

global r_ h_ sigma_ k_ c dh_ K_ R_ capon hsaton nrandon mrates Gamma_ delta_ muton beta_;
global b eps_ mu_ Pdim1 Ldim1 x0 chi_ Qstep Nstep tgone ntgone gammas1D lambdas1D;

frametimes = [325;350;375];
runnum = 1;
basecode = 'pldyn';
Pdim1 = 400;
Ldim1 = 400;

datapath = ['/Users/kimberly/Google Drive/immunedata/PL/' basecode '/']; %KONG datapath
%datapath = ['C:\Users\Kimberly\Google Drive\immunedata\PL\' basecode '\']; %laptop datapath

bfilename = [datapath 'b' basecode num2str(runnum) '.txt'];
tfilename = [datapath 't' basecode num2str(runnum) '.txt'];
Pfilename = [datapath 'P' basecode num2str(runnum) '.txt'];
Lfilename = [datapath 'L' basecode num2str(runnum) '.txt'];

params = setparams(bfilename);
days = params{end,2};

nframes = size(frametimes,1);
for i=1:nframes
    timevec = csvread(tfilename);
    if(frametimes(i)<days)
        tin = find(timevec>frametimes(i));
        t0index = tin(1)-1;
    else
        t0index = size(timevec,1);
    end
    P = transpose(csvread(Pfilename,t0index-1,0,[t0index-1,0,t0index-1,Pdim1-1]));
    L = transpose(csvread(Lfilename,t0index-1,0,[t0index-1,0,t0index-1,Ldim1-1]));

    figure    % plot of P0 and L0 distributions at frametimes(i)
hold on
hold all
plot((1:1:400),P)
plot((1:1:400),L)
title(['P and L distributions at t = ' num2str(frametimes(i)) ' days'])
xlabel('location in shape space (site)')
ylabel('population (cells/\mul)')
legend('Pathogen','Lymphocytes')
hold off
    
end


    
    

    