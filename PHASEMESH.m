% PHASEMESH

brange = 13;
chi_range = 200:100:1000;
%brange2 = 5:3:14;
%chi_range2 = 200:200:800;
% jumptest = 1;
% windowsize = 0.25;   % should be smaller than jumptest
% maxtests = 25; % max number of tests per b-value
seednum = 12;
seedbasecode = 'qtune';

realnum = 0;
realbasecode = 'phamesc';
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
        dlmwrite(savefile,[bb,chi_try,didescape,realnum],'-append');
        
        realnum = realnum+1;
     
    end   
    
    
end

% for bb=brange2
%     
%     realnum = realnum + 100 - mod(realnum,100);
%     for chi_try=chi_range2
%         
%         % create new paramfile from bseedfile & specified changes
%         params = setparams(bseedfile);
%         for i=1:size(params,1)
%             if ~strcmp(char(params{i,1}),'days')
%                 eval([char(params{i,1}) ' = ' num2str(params{i,2}) ';']);
%                 eval([char(params{i,1}) 'units = char(params{i,3});']);
%             end
%         end
%         clear params;
%         b0 = [r_;h_;sigma_;k_;bb;eps_;mu_;dh_;K_;R_;capon;hsaton;...
%             Pdim1;Ldim1;x0;chi_try;Gamma_;nrandon;delta_;spliton;pinit];
%         temppath = ['/Users/kimberly/Google Drive/immunedata/PL13/'...
%             realbasecode '/b' realbasecode '999.txt'];
%         writeparams(temppath,b0,temppath);
%         
%         % call ss_seed from seedfile run with new paramfile
%         didescape = ss_seed([seedbasecode num2str(seednum)],...
%            [realbasecode '999'],0,300,0.1,realbasecode,realnum,1);
%        
%         %save params, runnum, and result
%         dlmwrite(savefile,[bb,chi_try,didescape,realnum],'-append');
%         
%         realnum = realnum+1;
%      
%     end   
%     
%     
% end

mesh_tests_now = csvread(savefile);
mesh_tests1 = csvread(['/Users/kimberly/Google Drive/immunedata/PL13/'...
            'phamesa/tests.txt']);
mesh_tests2 = csvread(['/Users/kimberly/Google Drive/immunedata/PL13/'...
            'phames/tests.txt']); 
mesh_tests3 = csvread(['/Users/kimberly/Google Drive/immunedata/PL13/'...
            'phamesb/tests.txt']); 
mesh_tests4 = csvread(['/Users/kimberly/Google Drive/immunedata/PL13/'...
            'phlinb/tests.txt']); 
mesh_tests5 = csvread(['/Users/kimberly/Google Drive/immunedata/PL13/'...
            'phlina/tests.txt']); 
tests5 = [23,10,0,9.0;23,9,1,9.1;23,8,1,9;23,9.5,1,9.2;23,9.75,1,9.3];
tests6 = [24,10,0,10;24,9,0,10.1;24,8,1,10.2;24,8.5,1,10.3;24,8.75,0,10.4];
tests7 = [25,8,0,11;24,7,1,11.1;24,7.5,1,11.2;24,7.75,0,11.3];
tests8 = [26,6.5,1,12;26,6.75,0,12.1];
tests1 = csvread(['/Users/kimberly/Google Drive/immunedata/PL13/'...
            'phlind/tests.txt']);
tests2 = csvread(['/Users/kimberly/Google Drive/immunedata/PL13/'...
            'phlinc/tests.txt']);
tests3 = csvread(['/Users/kimberly/Google Drive/immunedata/PL13/'...
            'phline/tests.txt']);
tests4 = csvread(['/Users/kimberly/Google Drive/immunedata/PL13/'...
            'phlin/tests.txt']);
bchi_tests = [tests1;tests2;tests5;tests6;tests7;tests8;tests3;tests4]; 
mesh_tests = [mesh_tests1;mesh_tests2;mesh_tests3;mesh_tests4;mesh_tests5;...
   bchi_tests;mesh_tests_now];


figure
escapes = mesh_tests(:,2).*(mesh_tests(:,3)==1);
chronics = mesh_tests(:,2).*(~mesh_tests(:,3));
clears = mesh_tests(:,2).*(mesh_tests(:,3)==-1);
plot(mesh_tests(:,1),escapes,'xr',mesh_tests(:,1),chronics,'*b',...
    mesh_tests(:,1),clears,'.g')
%axis([0,50,5,1000])
        
        