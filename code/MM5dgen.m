% load 5dNCAT_orgIm x_roi
% x=mean(x_roi,5);clear x_roi
% yin=zeros(size(x));
% for g=1:8
%     yin(:,:,:,g)=convn(x(:,:,:,g),ones(3,3,3)/27,'same');
% end
% clear x
% M1=zeros(64,64,32,3,8);
% for i=1:7
%     [vx,vy,vz]=motionele3d(yin(:,:,:,i:i+1),50,.1);
%     M1(:,:,:,:,i)=cat(4,vx,vy,vz);
% end
% [vx,vy,vz]=motionele3d(yin(:,:,:,[8 1]),50,.1);
% M1(:,:,:,:,8)=cat(4,vx,vy,vz);
% save 5dNCAT_orgM1 M1
%load 5dNCAT_orgM1 M1
load 5dNCAT_motion
for g=1:8
    MM=mf2matrix3d(M,g,1);%M only
    filename=['5dMM' num2str(g)];
    save(filename,'MM');
end

%%%%%%%%%%%%%%
%generate (I-M)'*(I-M)=Q for exact motion prior!
load 5dNCAT_orgM1 M1
Q=0;
for g=1:8
    MM=mf2matrix3d(M1,g,1);%(I-M);
    Q1=MM'*MM;
    Q=Q+Q1;
end
clear Q1 MM M1
%save 5dExactQ Q
%split Q
Q=Q';
for g=1:8
    Qs=Q(:,(g-1)*64*64*32+1:g*64*64*32);
    Qs=Qs';
    filename=['5dQs' num2str(g)];
    save(filename,'Qs');
end
%Elapsed time is 4566.138050 seconds.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Possible implementation for exact motion prior
%g=1, calculate Q=MM'*MM, then split it to 8 sub-matrix, store.
%g=2,..., calculate Q=MM'*MM, then split it to 8 sub-matrix, load
%sub-matrix and add new sub-matrix, store. Repeat untile 8 gate finished.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Row index longer than column index is more efficient on memory usage.


% tic,MMc=[];for n=1:8,m_name=['5dMM' num2str(n)];load(m_name);MMc=cat(1,MMc,MM);end,toc
% Elapsed time is 24.967374 seconds.
% >> whos
%   Name         Size                    Bytes  Class
% 
%   MM      131072x1048576            92274704  double array (sparse)
%   MMc    1048576x1048576           708837476  double array (sparse)
%   m_name       1x5                        10  char array
%   n            1x1                         8  double array
% 
% Grand total is 66060303 elements using 801112198 bytes
% >> tic,MMc*ones(1048576,1);toc
% Elapsed time is 0.616826 seconds.
% >> 0.6*64
% ans =   38.4000
% >> log2(1048576)
% ans =    20
% >> save 5dMMall MMc
% >> tic,speye(1048576)-MMc;toc
% Elapsed time is 7.008918 seconds.
% >> tic,load 5dMMall;toc
% Elapsed time is 9.744784 seconds.
% >> tic,MMc';toc
% Elapsed time is 3.910245 seconds.
% >> clear MMc
% >> clear ans
% >> tic,MM*MM';toc
% Elapsed time is 1.884669 seconds.
% >> 1.88*64
% ans =  120.3200
 