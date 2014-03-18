% PHASEMESH for b v. mu diagram, plosruns, etc.
% fleshing out cmeshtwo

brange = [50,48:-2:30];
mu_range = 1e-5:1e-6:2.8e-5;
seednum = 507;
seedbasecode = 'cmeshtwoA';

realnum = 2000;
realbasecode = 'cmeshtwoA';
bseedfile = ['/Users/kimberly/Google Drive/immunedata/PL13/' seedbasecode...
    '/b' seedbasecode num2str(seednum) '.txt' ];
savefile = ['/Users/kimberly/Google Drive/immunedata/PL13/'...
            realbasecode '/tests.txt'];

tic        
for bb=brange
    
    realnum = realnum + 100 - mod(realnum,100);
    for mu_try=mu_range
        chi_try = 2.2e-4/mu_try;
        
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
toc

mesh_tests_now = csvread(savefile);
dlmwrite(['/Users/kimberly/Google Drive/immunedata/PL13/'...
             realbasecode '_tests.txt'],mesh_tests_now,'-append');
mesh_tests = csvread(['/Users/kimberly/Google Drive/immunedata/PL13/'...
             realbasecode '_tests.txt']);

%mesh_tests = mesh_tests_now;
figure
mesh_tests(:,2) = 2.2e-4./mesh_tests(:,2);
escapes = mesh_tests(:,2).*(mesh_tests(:,3)==1);
chronics = mesh_tests(:,2).*(~mesh_tests(:,3));
clears = mesh_tests(:,2).*(mesh_tests(:,3)==-1);
plot(mesh_tests(:,1).*(escapes>0),escapes,'xr',...
    mesh_tests(:,1).*(chronics>0),chronics,'*b',...
    mesh_tests(:,1).*(clears>0),clears,'.g')
%axis([0,50,5,1000])
        
        
