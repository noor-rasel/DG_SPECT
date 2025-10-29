clear all;
clc

pidall = {'BV010'};

for pidnum = 1:1
     
    pid = pidall{pidnum} ;   

    load (['/Users/wenyuanqi/Desktop/clinical/full/NormalDataBase/' pid '/' pid '_sino_resp_nor_full.mat'])

    for aa = 1:64
       for rr = 1:7

           temp = sino_resp(:,:,aa,rr);
           counts(aa,rr) = sum(temp(:));

       end
    end

    countsum = sum(counts,2);
    angletime_resp = counts./repmat(countsum,[1,7]);

    save (['/Users/wenyuanqi/Desktop/clinical/full/NormalDataBase/' pid '/' pid '_angletime_resp_nor_full.mat'],'angletime_resp');


    clear counts temp countsum;
end
