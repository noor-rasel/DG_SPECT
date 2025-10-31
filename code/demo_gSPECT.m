%SPECT recon demo - Mingwu Jin, UTA, Nov. 01, 2014.
%refer to the publication "MingwuJin_4D_PMBSept2009.pdf"

%% ideal recon: SIMIND simulation no/min noise, attenuation, and scatter
%%best recon can be achieved.
tc=8e6;%total counts
G=16;%# of gates
sub_num=16;%# of subsets
it_num=10;%# of iteration
OF_tag=0;% whether to calculate the objective function (0: no (save time); 1: yes)
blur=1;% include the distance-dependent blurring

load ncatd64_air sino
load roi
Im_ideal=zeros(64,64,64,G);
for g=1:G%each gate recon (10 iterations) takes about 50s on a i-5 computer. 
    Im_ideal(:,:,:,1)=mbsrem4dv2(sino(:,:,:,g),repmat(roi,[1,1,64]),sub_num,it_num,0,0,0);
end
Im_ideal=Im_ideal*tc/sum(Im_ideal(:));
%% Sinogram from more realistic simulation using SIMIND
load ncatd64 sino_p sino_s sino_L
%nact64.mat: SIMIND simulated projection data (sino_P: primary counts in 
%peak-energy window; sino_s: scatter in peak-energy window; sino_L: counts
%in the lower-energy window for scatter estimatoin)
sinop=random('poiss',sino_p);
sinos=random('poiss',sino_s);
sinoT=sinop+sinos;
sinoL=random('poiss',sino_L);

%% FBP
Im_rawfbp=fbp_3dM(sino_p);%recon for 16 gates.
%for one gate use: Im_rawfbp=fbp_3dM(sinoT(:,:,:,1));
Im_fbp=clinicfilt3d_yyjin(Im_rawfbp,2.4,.4);
Im_fbp=Im_fbp*tc/sum(Im_fbp(:));
%Im_fbp=Im_fbp(23:52,16:43,29:48,:); %region containing the heart

%% 3D iterative. Spatial regularization - modified BSREM

gbeta=0;% temproal regulization weight (0, no temproal smoothing)
%%MAP-S no attenuation (AC) and scatter correction (SC),sbeta=[0 5e-4 1e-3 1e-2 3e-2];
sbeta=1e-3; %spatial regulization weight (gbeta=0 & sbeta=0 <=> OSEM)
Im_maps=mbsrem4dv2(sinoT,repmat(roi,[1,1,64,16]),sub_num,it_num,OF_tag,...
    sbeta,gbeta,blur,0,0,0); 

%%MAP-S with AC no SC, sbeta=[0 1e-4 1e-3 5e-3 9e-3];
%Im_maps=mbsrem4dv2(sinoT,repmat(roi,[1,1,64,16]),sub_num,it_num,OF_tag,sbeta(n),gbeta,blur,1,0,0);

%%MAP-S with AC and SC,sbeta=[0 5e-4 2e-3 6e-3 1e-2];
% sinoL=reshape(sinoL,[64^2 64 16])*4;
% tew_scat=buttlpf(sinoL,0.2*.634,3);
% tew_scat=reshape(tew_scat,[64 64 64 16]);
% Im_maps=mbsrem4dv2(sinoT,repmat(roi,[1,1,64,16]),sub_num,it_num,OF_tag,sbeta(n),gbeta,blur,1,tew_scat,0);

Im_maps=Im_maps*8e6/sum(Im_maps(:));

%% motion estimation and generation of prediction matrix
MM4dgen(Im_maps); %some big files will be saved for 4D recon.

%% 4D iterative. Sptiotemporal regularization
gbeta=5e-3; %or 1e-2
sbeta=5e-4;
Im_mapst=mbsrem4dv2(sinoT,repmat(roi,[1,1,64,16]),sub_num,it_num,OF_tag,...
    sbeta,gbeta,blur,0,0,0);
Im_mapst=Im_mapst*tc/sum(Im_mapst(:));
%AC and SC can be applied similar to 3D recon.
%AC only: geta=3e-3 or 6e-3; sbeta=1e-4;
%AC and SC: geta=4e-3 or 7e-3; sbeta=5e-4;