% PHASELINE: this script systematically calls ss_seed with different values
% of two parameters of interest, to plot the line between escape and
% chronic control on the model's phase diagram.

clear

% starting values of parameters & first seedfile
chi00_lo = 0.0037;  % definite lower bound on chi_ for first (lowest) b-value
brange = (13);
jumptest = 0.00025;
windowsize00 = 0.25e-3;   % should be smaller than jumptest
xresolution = 1;
maxtests = 25; % max number of tests per b-value
minvalue00 = 0; % absolute minimum value of gamma00 to be tested
seednum = 12;
seedbasecode = 'qtune';

realnum = 2600;
realbasecode = 'BCoofill';
bseedfile = ['/Users/kimberly/Google Drive/immunedata/PL13/' seedbasecode...
    '/b' seedbasecode num2str(seednum) '.txt' ];
savefile = ['/Users/kimberly/Google Drive/immunedata/PL13/'...
            realbasecode '/tests.txt'];

for bb=brange
   
    % choose and name next run
    chi00_hi = max(chi00_lo+jumptest,minvalue00);
    isbound = 0;
    realnum = realnum + 100 - mod(realnum,100);
    while ~isbound % first, find a lower bound of chi_
        
        chi00_try = chi00_hi;
        
        % create new paramfile from bseedfile & specified changes
        params = setparams(bseedfile);
        for i=1:size(params,1)
            if ~strcmp(char(params{i,1}),'days')
                eval([char(params{i,1}) ' = ' num2str(params{i,2}) ';']);
                eval([char(params{i,1}) 'units = char(params{i,3});']);
            end
        end
        clear params;
        b0 = [r_;h_;sigma_;k_;bb;eps_;mu_;dh_;K_;R_;capon;hsaton;...
            Pdim1;Ldim1;x0;1/chi00_try;Gamma_;nrandon;delta_;spliton;pinit];
        temppath = ['/Users/kimberly/Google Drive/immunedata/PL13/'...
            realbasecode '/b' realbasecode '999.txt'];
        writeparams(temppath,b0,temppath);
        
        disp(['trying ' num2str(chi00_try) '...']);
        pause;
        
        % call ss_seed from seedfile run with new paramfile
        didescape = ss_seed([seedbasecode num2str(seednum)],...
           [realbasecode '999'],0,300,0.1,realbasecode,realnum,1);
       
        if ~didescape
           isbound = 1; 
           %save params, runnum, and result
           dlmwrite(['/Users/kimberly/Google Drive/immunedata/PL13/'...
            realbasecode '/tests.txt'],[bb,1/chi00_try,didescape,realnum],'-append');
        
        disp(didescape);
        disp(['Moving on to narrowing window: ' num2str(chi00_lo) ',' num2str(chi00_hi)]);
        pause;
        
        else
           chi00_lo = chi00_try;
           chi00_hi = chi00_hi+jumptest;
           disp(['Still finding window: ' num2str(chi00_lo) ',' num2str(chi00_hi-jumptest)]);
           pause;
        end
           realnum = realnum + 1;
        
        % emergency stop   
        if(mod(realnum,100)>maxtests)
            isbound = 1;
        end
           
    end
    
    while abs(chi00_lo-chi00_hi)>windowsize00; % now, narrow the window appropriately
        
        % determine chi_try
        chi00_try = (chi00_lo+chi00_hi)/2;
        
        % create new paramfile from bseedfile & specified changes
        params = setparams(bseedfile);
        for i=1:size(params,1)
            if ~strcmp(char(params{i,1}),'days')
                eval([char(params{i,1}) ' = ' num2str(params{i,2}) ';']);
                eval([char(params{i,1}) 'units = char(params{i,3});']);
            end
        end
        clear params;
        b0 = [r_;h_;sigma_;k_;bb;eps_;mu_;dh_;K_;R_;capon;hsaton;...
            Pdim1;Ldim1;x0;1/chi00_try;Gamma_;nrandon;delta_;spliton;pinit];
        temppath = ['/Users/kimberly/Google Drive/immunedata/PL13/'...
            realbasecode '/b' realbasecode '999.txt'];
        writeparams(temppath,b0,temppath);
        
        % SECOND: call ss_seed from seedfile run with new paramfile
        didescape = ss_seed([seedbasecode num2str(seednum)],...
           [realbasecode '999'],0,300,0.1,realbasecode,realnum,1);
        
        % THIRD: save params, runnum, and result
        dlmwrite(['/Users/kimberly/Google Drive/immunedata/PL13/'...
            realbasecode '/tests.txt'],[bb,1/chi00_try,didescape,realnum],'-append');
               
        % FOURTH: choose new runname, bseedfile, runnum, and values
        realnum = realnum + 1;
        if didescape
            chi00_lo = chi00_try;
        else
            chi00_hi = chi00_try;
        end
                    
        % emergency stop
        if(mod(realnum,100)>maxtests)
            windowsize00 = abs(chi00_hi-chi00_lo);
        end

    end
    
    % save result window at THAT particular value of bb
    dlmwrite(['/Users/kimberly/Google Drive/immunedata/PL13/'...
            realbasecode '/result.txt'],[bb,1/chi00_hi,1/chi00_lo],'-append');
    %chi_hi = chi_hi; % absolute upper bound for next bb value
    
end

mesh_tests_now = csvread(['/Users/kimberly/Google Drive/immunedata/PL13/'...
            realbasecode '/tests.txt']);
result_now = csvread(['/Users/kimberly/Google Drive/immunedata/PL13/'...
            realbasecode '/result.txt']);
dlmwrite(['/Users/kimberly/Google Drive/immunedata/PL13/'...
             'bchi_tests.txt'],mesh_tests_now,'-append');
dlmwrite(['/Users/kimberly/Google Drive/immunedata/PL13/'...
             'bchi_phaseline.txt'],result_now,'-append');
         
         % first, we load the necessary files
orig_tests = csvread(['/Users/kimberly/Google Drive/immunedata/PL13/'...
            'bchi_tests.txt']);
    % .... and make a quick plot of all tests.
        ymax = max(orig_tests(:,2));
        ymin = windowsize00;
        yaxis = ymin:windowsize00:ymax;
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
        
        overone_tests = [orig_tests(:,1),1./orig_tests(:,2),orig_tests(:,3:4)];

figure
escapes = overone_tests(:,2).*(overone_tests(:,3)==1);
chronics = overone_tests(:,2).*(~overone_tests(:,3));
clears = overone_tests(:,2).*(overone_tests(:,3)==-1);
plot(overone_tests(:,1),escapes,'xr',...
     overone_tests(:,1),chronics,'*b',...
     overone_tests(:,1),clears,'og')
%axis([xmin xmax 1/ymax 1/ymin])
axis([xmin 20 1/ymax 0.02])
        
        