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
        %axis([xmin 40 ymin 1000])
        legend('early escape','chronic infection','early clearance')
        
% overone_tests = [orig_tests(:,1),1./orig_tests(:,2),orig_tests(:,3:4)];
% 
% figure
% escapes = overone_tests(:,2).*(overone_tests(:,3)==1);
% chronics = overone_tests(:,2).*(~overone_tests(:,3));
% clears = overone_tests(:,2).*(overone_tests(:,3)==-1);
% plot(overone_tests(:,1),escapes,'xr',overone_tests(:,1),chronics,'*b',...
%      overone_tests(:,1),clears,'.g')
% axis([xmin xmax 1/ymax 1/ymin])
% %axis([xmin 20 1/ymax 0.02])
        
% create a matrix with mesh as small as tests. 
meshmatrix = zeros(size(yaxis,2),size(xaxis,2));

% fill all test-spots in matrix with appropriate numbers (round to nearest window)
orig_tests(:,3) = orig_tests(:,3) + 2; % new key: escape = 3, chronic = 2, clear = 1;
for i=1:size(orig_tests,1)
    if (orig_tests(i,2)>=ymin)&&(orig_tests(i,1)>=xmin)
    meshy = floor((orig_tests(i,2)-ymin)/windowsize)+1;
    meshx = floor((orig_tests(i,1)-xmin)/xresolution)+1;
    meshmatrix(meshy,meshx) = orig_tests(i,3);
    end
end

% fill interstices (by xval)
[yzeros,xzeros] = find(meshmatrix~=0);
% for i=min(xzeros):max(xzeros)
%     spikes = yzeros(xzeros==i);
%     if spikes
%         meshmatrix(1:spikes(1)-1,i) = meshmatrix(spikes(1),i);
%         for j=1:size(spikes,1)-1
%             if meshmatrix(spikes(j),i)==meshmatrix(spikes(j+1),i)
%                 meshmatrix(spikes(j)+1:spikes(j+1)-1,i) = meshmatrix(spikes(j),i);
%             end
%         end
%         meshmatrix(spikes(end)+1:end,i) = meshmatrix(spikes(end),i);
%     else
%         meshmatrix(:,i) = meshmatrix(:,i-1);
%     end
% end

% contour this matrix.
vi = [3 2 1 0];
v = [3.5 2.5 1.5 0.5 -0.5];
figure
contourf(meshmatrix,vi)