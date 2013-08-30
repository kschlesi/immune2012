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

name = 1.131;
matchstring = char(regexp(num2str(name),'\D\w*','match'));
str2num(matchstring(2:end))+1