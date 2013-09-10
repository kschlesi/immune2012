% PHASELINE: this script systematically calls ss_seed with different values
% of two parameters of interest, to plot the line between escape and
% chronic control on the model's phase diagram.

clear

% starting values of parameters & first seedfile
chi_hi = 800;  % definite upper bound on chi_ for first (lowest) b-value
brange = 5:9;
jumptest = 50;
windowsize = 0.25;   % should be smaller than jumptest
maxtests = 25; % max number of tests per b-value
seednum = 12;
seedbasecode = 'qtune';

realnum = 0;
realbasecode = 'phlinb';
bseedfile = ['/Users/kimberly/Google Drive/immunedata/PL13/' seedbasecode...
    '/b' seedbasecode num2str(seednum) '.txt' ];
savefile = ['/Users/kimberly/Google Drive/immunedata/PL13/'...
            realbasecode '/tests.txt'];

for bb=brange
   
    % choose and name next run
    chi_lo = chi_hi-jumptest;
    isbound = 0;
    realnum = realnum + 100 - mod(realnum,100);
    while ~isbound % first, find a lower bound of chi_
        
        chi_try = chi_lo;
        
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
            Pdim1;Ldim1;x0;chi_try;Gamma_;nrandon;delta_;spliton;pinit];
        temppath = ['/Users/kimberly/Google Drive/immunedata/PL13/'...
            realbasecode '/b' realbasecode '999.txt'];
        writeparams(temppath,b0,temppath);
        
        % call ss_seed from seedfile run with new paramfile
        didescape = ss_seed([seedbasecode num2str(seednum)],...
           [realbasecode '999'],0,300,0.1,realbasecode,realnum,1);
        if didescape
           isbound = 1; 
           %save params, runnum, and result
           dlmwrite(['/Users/kimberly/Google Drive/immunedata/PL13/'...
            realbasecode '/tests.txt'],[bb,chi_try,didescape,realnum],'-append');           
        else
           chi_hi = chi_try;
           chi_lo = chi_lo-jumptest;
        end
           realnum = realnum + 1;
        
        % emergency stop   
        if(mod(realnum,100)>maxtests)
            isbound = 1;
        end
           
    end
    
    while abs(chi_lo-chi_hi)>windowsize; % now, narrow the window appropriately
        
        % determine chi_try
        chi_try = (chi_lo+chi_hi)/2;
        
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
            Pdim1;Ldim1;x0;chi_try;Gamma_;nrandon;delta_;spliton;pinit];
        temppath = ['/Users/kimberly/Google Drive/immunedata/PL13/'...
            realbasecode '/b' realbasecode '999.txt'];
        writeparams(temppath,b0,temppath);
        
        % SECOND: call ss_seed from seedfile run with new paramfile
        didescape = ss_seed([seedbasecode num2str(seednum)],...
           [realbasecode '999'],0,300,0.1,realbasecode,realnum,1);
        
        % THIRD: save params, runnum, and result
        dlmwrite(['/Users/kimberly/Google Drive/immunedata/PL13/'...
            realbasecode '/tests.txt'],[bb,chi_try,didescape,realnum],'-append');
                
        % FOURTH: choose new runname, bseedfile, runnum, and values
        realnum = realnum + 1;
        if didescape
            chi_lo = chi_try;
        else
            chi_hi = chi_try;
        end
                    
        % emergency stop
        if(mod(realnum,100)>maxtests)
            windowsize = chi_hi-chi_lo;
        end

    end
    
    % save result window at THAT particular value of bb
    dlmwrite(['/Users/kimberly/Google Drive/immunedata/PL13/'...
            realbasecode '/result.txt'],[bb,chi_lo,chi_hi],'-append');
    %chi_hi = chi_hi; % absolute upper bound for next bb value
    
end