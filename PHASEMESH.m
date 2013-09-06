% PHASEMESH

brange = 15:5:35;
chi_range = 40:20:100;
% jumptest = 1;
% windowsize = 0.25;   % should be smaller than jumptest
% maxtests = 25; % max number of tests per b-value
seednum = 12;
seedbasecode = 'qtune';

realnum = 0;
realbasecode = 'phames';
bseedfile = ['/Users/kimberly/Google Drive/immunedata/PL13/' seedbasecode...
    '/b' seedbasecode num2str(seednum) '.txt' ];
savefile = ['/Users/kimberly/Google Drive/immunedata/PL13/'...
            realbasecode '/tests.txt'];
        
for bb=brange
    
    realnum = realnum + 100 - mod(realnum,100);
    for chi_try=chi_range
        
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
       
        %save params, runnum, and result
        dlmwrite(['/Users/kimberly/Google Drive/immunedata/PL13/'...
         realbasecode '/tests.txt'],[bb,chi_try,didescape,realnum],'-append');
        
        realnum = realnum+1;
     
    end
    
    
    
    
end
        
        