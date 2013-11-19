%%%%% PHASEPLOTS 2 %%%%%
% in which plots are made of the phase diagrams that have appropriate
% colors and things.

% twovars = 'bchi';
% xvar = 'b';
% yvar = 'chi_';
% windowsize = 0.25; % y-resolution
% xresolution = 1;
twovars = 'Gammabeta';
xvar = 'one';
yvar = 'two';
windowsize = 0.25; % y-resolution
xresolution = 0.025e-4;

% first, we load the necessary files
orig_tests = csvread(['/Users/kimberly/Google Drive/immunedata/PL13/'...
            twovars '_tests.txt']);
    % .... and make a quick plot of all tests.
        ymax = max(orig_tests(:,2));
        ymin = windowsize;
        yaxis = ymin:windowsize:ymax;
        xmin = min(orig_tests(:,1));
        xmax = max(orig_tests(:,1));
        xaxis = xmin:xresolution:xmax;
        %disp([xmin xmax ymin ymax]);
        figure
        escapes = orig_tests(:,2).*(orig_tests(:,3)==1);
        chronics = orig_tests(:,2).*(~orig_tests(:,3));
        clears = orig_tests(:,2).*(orig_tests(:,3)==-1);
        plot(orig_tests(:,1).*(escapes>0),escapes,'xr',...
            orig_tests(:,1).*(chronics>0),chronics,'*b',...
            orig_tests(:,1).*(clears>0),clears,'.g')
        axis([xmin xmax ymin ymax])
        %axis([xmin 40 ymin 1000])
        legend('early escape','chronic infection','early clearance')
        
        
% load phaseline results from file        
% orig_result = csvread(['/Users/kimberly/Google Drive/immunedata/PL13/'...
%            twovars '_phaseline.txt']);
% orig_result1 = csvread(['/Users/kimberly/Google Drive/immunedata/PL13/'...
%            twovars '_phaseline1.txt']);
% % sort results individually...
% orig_result = sortentry(orig_result,'col',1);
% orig_result1 = sortentry(orig_result1,'col',1);

% better way of grabbing orig_results from orig_tests
%        (all info should already be in orig_tests)

% find places FOR EACH x-value where: (a) changes from 1 to 0
%                                  0R (b) changes from 0 to -1 ;
% put each of these pairs of y-points in aLINE or bLINE matrix
orig_tests = sortentry(orig_tests,'col',2);
orig_tests = sortentry(orig_tests,'col',1);
x_hist = histc(orig_tests(:,1),xaxis);
    % to convert x_hist index to x-value : 
    %           (index-1)*xresolution + xmin = x-value

% preallocate aLine, bLine
aLine = zeros(floor(size(orig_tests,1)/2),3);
bLine = zeros(size(aLine));
current_arow = 1;
current_brow = 1;
%current_crow = 1;

% for loop that searches each x-value for phaseline points & saves them to
% matrices aLine and bLine
cum_used_entries = 0;
for i=1:size(x_hist,1)
    if x_hist(i)
        x_chunk = orig_tests(1+cum_used_entries:x_hist(i)+cum_used_entries,1:3);
        grab_index = find(diff(x_chunk(:,3)));
        if grab_index 
        for j=1:size(grab_index,1)    
            if ( ~x_chunk(grab_index(j),3) && x_chunk(grab_index(j)+1,3)==1 ) ||...
               ( x_chunk(grab_index(j),3)==1 && ~x_chunk(grab_index(j)+1,3) )
                aLine(current_arow,:) = ...
                    [x_chunk(grab_index(j),1),x_chunk(grab_index(j),2),x_chunk(grab_index(j)+1,2)];
                current_arow = current_arow+1;
            else
                bLine(current_brow,:) = ...
                    [x_chunk(grab_index(j),1),x_chunk(grab_index(j),2),x_chunk(grab_index(j)+1,2)];
                current_brow = current_brow+1;
            end
        end
        end
        clear xchunk;
        clear grab_index;
    end
    cum_used_entries = cum_used_entries + x_hist(i);
end
% truncate and plot aLine, bLine
aLine = aLine(1:current_arow-1,:);
bLine = bLine(1:current_brow-1,:);
figure
plot(aLine(:,1),mean([aLine(:,2),aLine(:,3)],2),'go')
hold on
plot(bLine(:,1),mean([bLine(:,2),bLine(:,3)],2),'ro')



% preparing to fit a curve to the relevant phaseline points
orig_result = aLine;
mean_y = mean([orig_result(:,2),orig_result(:,3)],2);

% plot original guess of fit parameters
figure
plot(powerfun([10^5,3],orig_result(:,1)))

% interactively plot power-law fit
% nlintool(orig_result(:,1),mean_y,@powerfun,[10^5,3]);

% calculate power-law fit and error on fit
[bhat,resid,J] = nlinfit(orig_result(:,1),mean_y,@powerfun,[10^5,3]);
betaci = nlparci(bhat,resid,J);
disp(betaci(1,:));
disp(betaci(2,:));

% plot relevant test points (green), fit (blue), fit error lines (red)
figure
plot(orig_result(:,1),mean_y,'go')
hold on
plot(orig_result(:,1),powerfun(betaci(:,1),orig_result(:,1)),'r--')
plot(orig_result(:,1),powerfun(betaci(:,2),orig_result(:,1)),'r--')
plot(orig_result(:,1),powerfun(bhat,orig_result(:,1)),'b')


