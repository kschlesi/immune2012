%%%%% PHASEPLOTS %%%%%
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
        
orig_result = csvread(['/Users/kimberly/Google Drive/immunedata/PL13/'...
            twovars '_phaseline.txt']);
% orig_result1 = csvread(['/Users/kimberly/Google Drive/immunedata/PL13/'...
%             twovars '_phaseline1.txt']);
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
%       determine how many y-axis crossings exist (xcounts)
% 

% finally, we create the plot.
% figure out how many xes there are.
xcounts = histc(orig_result(:,1),xaxis)';
phasmat = zeros(size(yaxis,2),size(xaxis,2)); % includes all tested points
for xval=1:size(xaxis,2)
    num_lines = 0;
    while num_lines<xcounts(xval)
        phasmat(:,xval) = phasmat(:,xval) + ...
            ((yaxis')>=orig_result(xval,3)) - ((yaxis')<=orig_result(xval,2));
        num_lines = num_lines+1;
    end
end

figure
v = [-2 -1 0 1 2 3];
contourf(xaxis,(yaxis'),phasmat);
