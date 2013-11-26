%%%% phaseplot gammabeta

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
        plot(orig_tests(:,1).*(escapes>0),escapes,'or',...
            orig_tests(:,1).*(chronics>0),chronics,'ob',...
            orig_tests(:,1).*(clears>0),clears,'og')
        axis([xmin xmax ymin ymax])
        axis([xmin xmax ymin 100])
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
current_arow = 1;
current_brow = 1;
current_crow = 1;

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
                if ( ~x_chunk(grab_index(j),3) && x_chunk(grab_index(j)+1,3)==-1 ) ||...
               ( x_chunk(grab_index(j),3)==-1 && ~x_chunk(grab_index(j)+1,3) )
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
figure
plot(aLine(:,1),mean([aLine(:,2),aLine(:,3)],2),'g')
hold on
plot(bLine(:,1),mean([bLine(:,2),bLine(:,3)],2),'r')
plot(cLine(:,1),mean([cLine(:,2),cLine(:,3)],2),'b')

%%% now that we have the points, fit the lines!!!

%%%%%%%%%%%%%%%%%%aLine plot%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% preparing to fit a curve to the relevant phaseline points
keep_indx = find(aLine(:,2)~=max(aLine(:,2)));
orig_result1 = aLine(keep_indx,:);
mean_y1 = mean([orig_result1(:,2),orig_result1(:,3)],2);

fitrange1 = (orig_result1(1,1):xresolution:orig_result1(end,1));

% plot original guess of fit parameters
figure
%exponential(blue)
plot(fitrange1,expfun([10^2,3e5],(orig_result1(1,1):xresolution:orig_result1(end,1))),'b')
hold on
%powerlaw (red)
plot(fitrange1,powerfun([10^-14/4,3],(orig_result1(1,1):xresolution:orig_result1(end,1))),'r')

% calculate fit and error on fit
[bhat1,resid,J] = nlinfit(orig_result1(:,1),mean_y1,@powerfun,[10^-14/4,3]);
betaci1 = nlparci(bhat1,resid,J);
[bhat1_e,resid_e,J_e] = nlinfit(orig_result1(:,1),mean_y1,@expfun,[10^2,3e5]);
betaci1_e = nlparci(bhat1_e,resid_e,J_e);
%disp(betaci1(1,:));
%disp(betaci1(2,:));

% plot relevant test points (green), fit (blue), fit error lines (red)
figure
plot(orig_result1(:,1),mean_y1,'go')
hold on
%plot(orig_result1(:,1),powerfun(betaci1(:,1),orig_result1(:,1)),'r--')
%plot(orig_result1(:,1),powerfun(betaci1(:,2),orig_result1(:,1)),'r--')
plot(fitrange1,powerfun(bhat1,fitrange1),'r')
plot(fitrange1,expfun(bhat1_e,fitrange1),'b')

%%%%%%%%%%%%%%%%%%%%%%%%%bLine plot%%%%%%%%%%%%%%%%%%%%%%%%%

% % preparing to fit a curve to the relevant phaseline points
% orig_result2 = bLine;
% mean_y2 = mean([orig_result2(:,2),orig_result2(:,3)],2);
% 
% fitrange2 = (orig_result2(1,1):xresolution:orig_result2(end,1));
% 
% % plot original guess of fit parameters
% figure
% %exponential (blue)
% plot(fitrange2,expfun([10^3,3e5],fitrange2),'b')
% hold on 
% %powerlaw (red)
% plot(fitrange2,powerfun([10^-8,2],fitrange2),'r')
% 
% % calculate power-law fit and error on fit
% % ...or exponential??
% [bhat2,resid,J] = nlinfit(orig_result2(:,1),mean_y2,@powerfun,[10^-8,2]);
% betaci2 = nlparci(bhat2,resid,J);
% [bhat2_e,resid_e,J_e] = nlinfit(orig_result2(:,1),mean_y2,@expfun,[10^3,3e5]);
% betaci2_e = nlparci(bhat2_e,resid_e,J_e);
% %disp(betaci2(1,:));
% %disp(betaci2(2,:));
% 
% % plot relevant test points (green), fit (blue), fit error lines (red)
% figure
% plot(orig_result2(:,1),mean_y2,'ro')
% hold on
% %plot(orig_result2(:,1),expfun(betaci2(:,1),orig_result2(:,1)),'r--')
% %plot(orig_result2(:,1),expfun(betaci2(:,2),orig_result2(:,1)),'r--')
% plot(fitrange2,expfun(bhat2_e,fitrange2),'b')
% plot(fitrange2,powerfun(bhat2,fitrange2),'r')

%%%%%%%%%%%%%bLine plot (bottom half)%%%%%%%%%%%%%%

% preparing to fit a curve to the relevant phaseline points
keep_indx = find((bLine(:,2)~=max(bLine(:,2))).*(abs(bLine(:,3)-bLine(:,2)<1)));
orig_resultB = bLine(keep_indx,:);
orig_result3 = orig_resultB(find(orig_resultB(:,2)<30),:);
mean_y3 = mean([orig_result3(:,2),orig_result3(:,3)],2);

fitrange3 = (orig_result3(1,1):xresolution:orig_result3(end,1));

% plot original guess of fit parameters
figure
%exponential (blue)
plot(fitrange3,expfun([10^3,3e5],fitrange3),'b')
hold on 
%powerlaw (red)
plot(fitrange3,powerfun([10^-8,2],fitrange3),'r')

% calculate power-law fit and error on fit
% ...or exponential??
[bhat3,resid,J] = nlinfit(orig_result3(:,1),mean_y3,@powerfun,[10^-8,2]);
betaci3 = nlparci(bhat3,resid,J);
[bhat3_e,resid_e,J_e] = nlinfit(orig_result3(:,1),mean_y3,@expfun,[10^3,3e5]);
betaci3_e = nlparci(bhat3_e,resid_e,J_e);
% disp(betaci3(1,:));
% disp(betaci3(2,:));

% plot relevant test points (green), fit (blue), fit error lines (red)
figure
plot(orig_result3(:,1),mean_y3,'ro')
hold on
%plot(orig_result3(:,1),powerfun(betaci3(:,1),orig_result3(:,1)),'r--')
%plot(orig_result3(:,1),powerfun(betaci3(:,2),orig_result3(:,1)),'r--')
plot(fitrange3,powerfun(bhat3,fitrange3),'r')
plot(fitrange3,expfun(bhat3_e,fitrange3),'b')

%%%%%%%%%%%%%bLine plot (top half)%%%%%%%%%%%%%%

orig_result5 = orig_resultB(find(orig_resultB(:,2)>20),:);
mean_y5 = mean([orig_result5(:,2),orig_result5(:,3)],2);

fitrange5 = (orig_result5(1,1):xresolution:orig_result5(end,1));

% plot original guess of fit parameters
figure
%exponential (blue)
%plot(fitrange5,expfun([10^3,3e5],fitrange5),'b')
%hold on 
%powerlaw (red)
plot(fitrange5,logfun([1,20],fitrange5),'r')

% calculate power-law fit and error on fit
% ...or exponential??
[bhat5,resid,J] = nlinfit(orig_result5(:,1),mean_y5,@logfun,[10^-8,2]);
betaci5 = nlparci(bhat5,resid,J);
%[bhat5_e,resid_e,J_e] = nlinfit(orig_result5(:,1),mean_y5,@expfun,[10^3,3e5]);
%betaci5_e = nlparci(bhat5_e,resid_e,J_e);
% disp(betaci3(1,:));
% disp(betaci3(2,:));

% plot relevant test points (green), fit (blue), fit error lines (red)
figure
plot(orig_result5(:,1),mean_y5,'ro')
hold on
%plot(orig_result3(:,1),powerfun(betaci3(:,1),orig_result3(:,1)),'r--')
%plot(orig_result3(:,1),powerfun(betaci3(:,2),orig_result3(:,1)),'r--')
plot(fitrange5,logfun(bhat5,fitrange5),'r')
%plot(fitrange5,expfun(bhat5_e,fitrange5),'b')


%%%%%%%%%%%%%%%%%%%bLine plot (sideways)%%%%%%%%%%%%%%%%%%%%

% preparing to fit a curve to the relevant phaseline points
% sideways_y = bLine(:,1);
% sideways_x = mean([bLine(:,2),bLine(:,3)],2);
% 
% fitrange4 = (min(sideways_x):windowsize:max(sideways_x));
% 
% % % plot original guess of fit parameters
% % figure
% % %exponential (blue)
% % plot(fitrange4,expfun([10^3,3e5],fitrange4),'b')
% % hold on 
% % %powerlaw (red)
% % plot(fitrange4,powerfun([10^-8,2],fitrange4),'r')
% 
% % calculate power-law fit and error on fit
% % ...or exponential??
% 
% %[theFit,goFit] = fit(sideways_x,sideways_y,'poly2');
% %plot(theFit,sideways_x,sideways_y);
% 
% %figure
% %plot(sideways_x,sideways_y,'go')
% %hold on
% %theFit = polyfit(sideways_x,sideways_y,2);
% %fitPlotVals = polyval(theFit,fitrange4);
% %plot(fitPlotVals,fitrange4,'b')
% 
% 
% [bhat2,resid,J] = nlinfit(orig_result4(:,1),mean_y4,@powerfun,[10^-8,2]);
% betaci2 = nlparci(bhat2,resid,J);
% [bhat2_e,resid_e,J_e] = nlinfit(orig_result4(:,1),mean_y4,@expfun,[10^3,3e5]);
% betaci2_e = nlparci(bhat2_e,resid_e,J_e);
%disp(betaci2(1,:));
%disp(betaci2(2,:));

% % plot relevant test points (green), fit (blue), fit error lines (red)
% figure
% plot(orig_result4(:,1),mean_y4,'ro')
% hold on
% %plot(orig_result2(:,1),expfun(betaci2(:,1),orig_result2(:,1)),'r--')
% %plot(orig_result2(:,1),expfun(betaci2(:,2),orig_result2(:,1)),'r--')
% plot(fitrange4,expfun(bhat2_e,fitrange4),'b')
% plot(fitrange4,powerfun(bhat2,fitrange4),'r')


% truncate and plot aLine, bLine, &c
figure
plot(orig_result1(:,1),mean([orig_result1(:,2),orig_result1(:,3)],2),'b')
hold on
%plot(fitrange1,powerfun(bhat1,fitrange1),'k')
plot(orig_result3(:,1),mean([orig_result3(:,2),orig_result3(:,3)],2),'r')
%plot(fitrange3,powerfun(bhat3,fitrange3),'k')
plot(orig_result5(:,1),mean([orig_result5(:,2),orig_result5(:,3)],2),'r')
axis([xmin xmax ymin 80])
