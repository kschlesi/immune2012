clear

days = 4444;
runnum = 1;
basecode = 'pldyn';
datapath = ['/Users/kimberly/Google Drive/immunedata/PL/' basecode '/'];
bfilename = [datapath 'b' basecode num2str(runnum) '.txt'];
tfilename = [datapath 't' basecode num2str(runnum) '.txt'];
Pfilename = [datapath 'P' basecode num2str(runnum) '.txt'];
Lfilename = [datapath 'L' basecode num2str(runnum) '.txt'];

% set parameters
params = setparams(bfilename);
for i=1:size(params,1)
    if ~strcmp(char(params{i,1}),'days')
        eval([char(params{i,1}) ' = ' num2str(params{i,2}) ';']);
        eval([char(params{i,1}) 'units = char(params{i,3});']);
    end
end
clear params;
% 
% Lflux = Gamma_;
% disp(size(Lflux));
% dL = Lflux + ones(Ldim1,1);
% disp(dL(1:3));
% 
% mplot = randn(5);
% disp(mplot);
% mplot = mplot.*(mplot<1).*(mplot>-1);
% disp(mplot);

% make a matrix of size (Pdim1,Pdim1) whose entry (i,j) = i-j
no = 10;
imat = repmat((1:no)',1,no); % each index names its own row
jmat = repmat(1:no,no,1); % each index names its own column
%disp(imat);
%disp(jmat);
%disp(abs(imat-jmat));

iminj = abs(repmat((1:Pdim1)',1,Pdim1)-repmat(1:Pdim1,Pdim1,1));

%Pdim1 = 5;
%mrates = Qmatrix(Pdim1,chi_);
% P = (1:Pdim1)';
% 
% tic
% Pmat = repmat(P,1,Pdim1);
% dmut1 = sum(Pmat.*mrates,1)';
% clear Pmat;
% toc
% 
% tic
% dmut2 = zeros(Pdim1,1);
% for i=1:Pdim1
%     dmut2(i) = sum(P.*mrates(:,i));
% end
% toc

%disp(dmut1);
%disp(dmut2);

% chii = 2.5*10^-4;
% figure
% plot(sqrt(2/pi)*chii*(1:400).^(-2))
% figure
% loglog(sqrt(2/pi)*chii*(1:400).^(-2))

% name = 1.131;
% matchstring = char(regexp(num2str(name),'\D\w*','match'));
% str2num(matchstring(2:end))+1

% result1 = csvread(['/Users/kimberly/Google Drive/immunedata/PL13/'...
%             'phlind/result.txt']);
% result2 = csvread(['/Users/kimberly/Google Drive/immunedata/PL13/'...
%             'phlinc/result.txt']);
% result3 = csvread(['/Users/kimberly/Google Drive/immunedata/PL13/'...
%             'phline/result.txt']);
% result4 = csvread(['/Users/kimberly/Google Drive/immunedata/PL13/'...
%             'phlin/result.txt']);
% result5 = [23,9.5,10;24,8.5,8.75;25,7.5,7.75;26,6.5,6.75];
% tests5 = [23,10,0,9.0;23,9,1,9.1;23,8,1,9;23,9.5,1,9.2;23,9.75,1,9.3];
% tests6 = [24,10,0,10;24,9,0,10.1;24,8,1,10.2;24,8.5,1,10.3;24,8.75,0,10.4];
% tests7 = [25,8,0,11;24,7,1,11.1;24,7.5,1,11.2;24,7.75,0,11.3];
% tests8 = [26,6.5,1,12;26,6.75,0,12.1];
% tests1 = csvread(['/Users/kimberly/Google Drive/immunedata/PL13/'...
%             'phlind/tests.txt']);
% tests2 = csvread(['/Users/kimberly/Google Drive/immunedata/PL13/'...
%             'phlinc/tests.txt']);
% tests3 = csvread(['/Users/kimberly/Google Drive/immunedata/PL13/'...
%             'phline/tests.txt']);
% tests4 = csvread(['/Users/kimberly/Google Drive/immunedata/PL13/'...
%             'phlin/tests.txt']);
        
bchi_result = csvread(['/Users/kimberly/Google Drive/immunedata/PL13/'...
            'bchi_phaseline.txt']);
bchi_result1 = csvread(['/Users/kimberly/Google Drive/immunedata/PL13/'...
            'bchi_phaseline1.txt']);
bchi_tests = csvread(['/Users/kimberly/Google Drive/immunedata/PL13/'...
            'bchi_tests.txt']);
bchi_result1 = sortentry(bchi_result1,'col',1);
disp(bchi_result1);

windowsize = 0.25;
ymax = max(bchi_tests(:,2));
ymin = windowsize;
yaxis = ymin:windowsize:ymax;
xmin = min(bchi_tests(:,1));
xmax = max(bchi_tests(:,1));
xaxis = xmin:xmax;
%disp([xmin xmax ymin ymax]);

figure
escapes = bchi_tests(:,2).*(bchi_tests(:,3)==1);
chronics = bchi_tests(:,2).*(~bchi_tests(:,3));
clears = bchi_tests(:,2).*(bchi_tests(:,3)==-1);
plot(bchi_tests(:,1),escapes,'xr',bchi_tests(:,1),chronics,'*b',bchi_tests(:,1),clears,'.g')
axis([xmin xmax ymin ymax])

phasmat = zeros(size(yaxis,2),size(bchi_result,1));
chilist = yaxis';
for brow=1:size(bchi_result,1)
    phasmat(:,brow) = (chilist>=bchi_result(brow,3)) - (chilist<=bchi_result(brow,2));
end

figure
v = [0 1];
contourf(bchi_result(:,1),chilist,phasmat,v);

% plot 1/chi_ instead
overone_result = [bchi_result(:,1),1./bchi_result(:,3),1./bchi_result(:,2)];
overone_tests = [bchi_tests(:,1),1./bchi_tests(:,2),bchi_tests(:,3:4)];

figure
escapes = overone_tests(:,2).*overone_tests(:,3)==1;
chronics = overone_tests(:,2).*(~overone_tests(:,3));
clears = overone_tests(:,2).*overone_tests(:,3)==-1;
plot(overone_tests(:,1),escapes,'xr',...
     overone_tests(:,1),chronics,'*b',...
     overone_tests(:,1),clears,'og')
axis([xmin xmax 1/ymax 1/ymin])

chilist = 1./yaxis';
phasmat = zeros(size(chilist,1),size(overone_result,1));
for brow=1:size(overone_result,1)
    phasmat(:,brow) = (chilist>=overone_result(brow,3)) - (chilist<=overone_result(brow,2));
end

figure
v = [0 1];
contourf(bchi_result(:,1),chilist,phasmat,v);
