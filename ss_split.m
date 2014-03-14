%%% runs several ss_seeds in a split

seedbasecode = 'cmeshone';
seednum = 516;
realbasecode = 'cmeshextA';

savefile = ['/Users/kimberly/Google Drive/immunedata/PL13/'...
            realbasecode '/list.txt'];
realnum = seednum*100;
didesc = ss_seed([seedbasecode num2str(seednum)],...
            [seedbasecode num2str(seednum)],'end',500,0.1,realbasecode,realnum,1);
dlmwrite(savefile,[realnum,didesc],'-append');

tic
while didesc~=1
    
    didesc = ss_seed([realbasecode num2str(realnum)],...
            [realbasecode num2str(realnum)],'end',500,0.1,realbasecode,realnum+1,1);
    dlmwrite(savefile,[realnum+1,didesc],'-append');
    realnum = realnum+1;
        
end
toc