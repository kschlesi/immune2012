%%%%% PHASEPLOTS 2 %%%%% FOR B MU PHASE DIAGRAMS
% in which plots are made of the phase diagrams that have appropriate
% colors and things.

twovars = 'cmeshtwoA';
windowsize = 1e-6; % y-resolution (MU, not chi_)
xresolution = 1;
% twovars = 'Gammabeta';
% xvar = 'one';
% yvar = 'two';
% windowsize = 0.25; % y-resolution
% xresolution = 0.025e-4;

% first, we load the necessary files
%orig_tests = csvread(['/Users/kimberly/Google Drive/immunedata/PL13/'...
%            twovars '/tests.txt']);
        orig_tests = csvread('tests.txt');        
        orig_tests(:,2) = 2.2e-4./orig_tests(:,2);
    % .... and make a quick plot of all tests.
        ymax = max(orig_tests(:,2));
        ymin = min(orig_tests(:,2));
        yaxis = ymin:windowsize:ymax;
        xmin = min(orig_tests(:,1));
        xmax = max(orig_tests(:,1));
        xaxis = xmin:xresolution:xmax;

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


% find places FOR EACH x-value where: (a) changes from 1 to 0
%                                  0R (b) changes from 0 to -1 ;
% put each of these pairs of y-points in aLINE or bLINE matrix
orig_tests = sortentry(orig_tests,'col',2);
orig_tests = sortentry(orig_tests,'col',1);
x_hist = histc(orig_tests(:,1),xaxis);
    % to convert x_hist index to x-value : 
    %           (index-1)*xresolution + xmin = x-value

% preallocate aLine, bLine
aLine = zeros(floor(size(orig_tests,1)/3),3);
bLine = zeros(size(aLine));
cLine = zeros(size(aLine));
dLine = zeros(size(aLine));
eLine = zeros(size(aLine));
current_arow = 1;
current_brow = 1;
current_crow = 1;
current_drow = 1;
current_erow = 1;

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
                if x_chunk(grab_index(j),1)==x_chunk(grab_index(j)+1,1) && x_chunk(grab_index(j),2)
                aLine(current_arow,:) = ...
                    [x_chunk(grab_index(j),1),x_chunk(grab_index(j),2),x_chunk(grab_index(j)+1,2)];
                current_arow = current_arow+1;
                end
            else
                if ( ~x_chunk(grab_index(j),3) && x_chunk(grab_index(j)+1,3)==-1 ) 
                if x_chunk(grab_index(j),1)==x_chunk(grab_index(j)+1,1) && x_chunk(grab_index(j),2)
                bLine(current_brow,:) = ...
                    [x_chunk(grab_index(j),1),x_chunk(grab_index(j),2),x_chunk(grab_index(j)+1,2)];
                current_brow = current_brow+1;
                end
                else
                if x_chunk(grab_index(j),1)==x_chunk(grab_index(j)+1,1) && x_chunk(grab_index(j),2)
                    cLine(current_crow,:) = ...
                        [x_chunk(grab_index(j),1),x_chunk(grab_index(j),2),x_chunk(grab_index(j)+1,2)];
                    current_crow = current_crow+1;
                end
                end
                
                if ( x_chunk(grab_index(j),3)==-1 && ~x_chunk(grab_index(j)+1,3) ) 
                if x_chunk(grab_index(j),1)==x_chunk(grab_index(j)+1,1) && x_chunk(grab_index(j),2)
                dLine(current_drow,:) = ...
                    [x_chunk(grab_index(j),1),x_chunk(grab_index(j),2),x_chunk(grab_index(j)+1,2)];
                current_drow = current_drow+1;
                end
                else
                if x_chunk(grab_index(j),1)==x_chunk(grab_index(j)+1,1) && x_chunk(grab_index(j),2)
                    eLine(current_erow,:) = ...
                        [x_chunk(grab_index(j),1),x_chunk(grab_index(j),2),x_chunk(grab_index(j)+1,2)];
                    current_erow = current_erow+1;
                end
                end
               
            end
        end
        end
        clear x_chunk;
        clear grab_index;
    end
    cum_used_entries = cum_used_entries + x_hist(i);
end
% truncate and plot aLine, bLine, &c
aLine = aLine(1:current_arow-1,:);
bLine = bLine(1:current_brow-1,:);
cLine = cLine(1:current_crow-1,:);
dLine = dLine(1:current_drow-1,:);
eLine = eLine(1:current_erow-1,:);
figure
plot(aLine(:,1),mean([aLine(:,2),aLine(:,3)],2),'go')
hold on
plot(bLine(:,1),mean([bLine(:,2),bLine(:,3)],2),'ro')
plot(cLine(:,1),mean([cLine(:,2),cLine(:,3)],2),'bo')
plot(dLine(:,1),mean([dLine(:,2),dLine(:,3)],2),'ko')
plot(eLine(:,1),mean([eLine(:,2),eLine(:,3)],2),'bo')


tomoves=(dLine(:,1)>35).*(dLine(:,2)<2e-5);
d1Line = dLine(tomoves>0,:);
d2Line = dLine(~tomoves,:);

figure
plot(aLine(:,1),mean([aLine(:,2),aLine(:,3)],2),'k-')
hold on
plot(bLine(:,1),mean([bLine(:,2),bLine(:,3)],2),'k-')
plot(d1Line(:,1),mean([d1Line(:,2),d1Line(:,3)],2),'k-')
plot(d2Line(:,1),mean([d2Line(:,2),d2Line(:,3)],2),'k-')


