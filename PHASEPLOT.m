%%%%% PHASEPLOTS %%%%%
% in which plots are made of the phase diagrams that have appropriate
% colors and things.

twovars = 'bchi';
xvar = 'b';
yvar = 'chi_';
windowsize = 0.25;

% first, we load the necessary files
orig_tests = csvread(['/Users/kimberly/Google Drive/immunedata/PL13/'...
            twovars '_tests.txt']);
    % .... and make a quick plot of all tests.
        ymax = max(orig_tests(:,2));
        ymin = windowsize;
        yaxis = ymin:windowsize:ymax;
        xmin = min(orig_tests(:,1));
        xmax = max(orig_tests(:,1));
        xaxis = xmin:xmax;
        %disp([xmin xmax ymin ymax]);
        figure
        escapes = orig_tests(:,2).*(orig_tests(:,3)==1);
        chronics = orig_tests(:,2).*(~orig_tests(:,3));
        clears = orig_tests(:,2).*(orig_tests(:,3)==-1);
        plot(orig_tests(:,1).*(escapes>0),escapes,'xr',...
        orig_tests(:,1).*(chronics>0),chronics,'*b',...
        orig_tests(:,1).*(clears>0),clears,'.g')
        %axis([xmin xmax ymin ymax])
        axis([xmin 40 ymin 1000])
        legend('early escape','chronic infection','early clearance')
        
orig_result = csvread(['/Users/kimberly/Google Drive/immunedata/PL13/'...
            twovars '_phaseline.txt']);
orig_result1 = csvread(['/Users/kimberly/Google Drive/immunedata/PL13/'...
            twovars '_phaseline1.txt']);
% sort results individually...
%orig_result = sortentry(orig_result,'col',1);
%orig_result1 = sortentry(orig_result1,'col',1);

% or, put all results together and sort them that way.
orig_result = [orig_result;orig_result1];
orig_result = sortentry(orig_result,'col',1);
%disp(orig_result);


% next, we define the separate regions between the phaselines.
% like so:
% FOR EACH x-axis number (listed in orig_tests(:,1))
%       determine how many y-axis crossings exist
% 

% finally, we create the plot.
