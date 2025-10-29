%cho 2 channels

% load cho_result_fbp_140_g1f8_center
% load cho_result_st121_140_g1f8_center
% load cho_result_nAS140_g1f8_center
% load cho_result_A140_g1f8_center
% load cho_result_AS140_g1f8_center

% load cho_result_fbp_140_g1f8_centerlu
% load cho_result_st121_140_g1f8_centerlu
% load cho_result_nAS140_g1f8_centerlu
% load cho_result_A140_g1f8_centerlu
% load cho_result_AS140_g1f8_centerlu

load cho_result_fbp_140_g1f8_centerluN
load cho_result_st121_140_g1f8_centerluN
load cho_result_nAS140_g1f8_centerluN
load cho_result_A140_g1f8_centerluN
load cho_result_AS140_g1f8_centerluN

int_noise=(1/256)^2/12;ch_numb=2;
for s=1:S
    for t=1:T
        Az_nAs(s,t)=Azfinder(squeeze(ch_nAs_nor(1:ch_numb,:,s,t))',squeeze(ch_nAs_les(1:ch_numb,:,s,t))',int_noise);
    end
end
for s=1:S
    for t=1:T 
        Az_As(s,t)=Azfinder(squeeze(ch_As_nor(1:ch_numb,:,s,t))',squeeze(ch_As_les(1:ch_numb,:,s,t))',int_noise);
    end
end
for s=1:S
    for t=1:T        
        Az_A(s,t)=Azfinder(squeeze(ch_A_nor(1:ch_numb,:,s,t))',squeeze(ch_A_les(1:ch_numb,:,s,t))',int_noise);
    end
end
Az_st121=Azfinder(ch_st121_nor(1:ch_numb,:)',ch_st121_les(1:ch_numb,:)',int_noise);
Az_fbp=Azfinder(ch_fbp_nor(1:ch_numb,:)',ch_fbp_les(1:ch_numb,:)',int_noise);

az_fbp=Az_fbp;
az_st121=Az_st121;%0.5*(1+erf(snr_st121/2));
az_mapnAS=Az_nAs';
az_mapA=Az_A';%0.5*(1+erf(snr_mapA/2));
az_mapAS=Az_As';%0.5*(1+erf(snr_mapAS/2));

sbeta=[0 5e-4 1e-3 1e-2 3e-2];
figure,plot(sbeta,az_mapnAS(1,:),'*-',sbeta,az_mapnAS(2,:),'+-',sbeta,az_mapnAS(3,:),'o-')
hold on,plot(sbeta,repmat(az_fbp,1,5),'k-',sbeta,repmat(az_st121,1,5),'k--')
legend('\beta_t=0','\beta_t=1\times 10^{-2}','\beta_t=2\times 10^{-2}','FBP','ST 121')
title('NC'),xlabel('\beta_s'),ylabel('Az')

sbeta=[0 5e-5 1e-4 1e-3 3e-3];
figure,plot(sbeta,az_mapA(1,:),'*-',sbeta,az_mapA(2,:),'+-',sbeta,az_mapA(3,:),'o-')
hold on,plot(sbeta,repmat(az_fbp,1,5),'k-',sbeta,repmat(az_st121,1,5),'k--')
legend('\beta_t=0','\beta_t=1\times 10^{-3}','\beta_t=2\times 10^{-3}','FBP','ST 121')
title('AC'),xlabel('\beta_s'),ylabel('Az')

sbeta=[0 5e-4 1e-3 2e-3 4e-3];
figure,plot(sbeta,az_mapAS(1,:),'*-',sbeta,az_mapAS(2,:),'+-',sbeta,az_mapAS(3,:),'o-')
hold on,plot(sbeta,repmat(az_fbp,1,5),'k-',sbeta,repmat(az_st121,1,5),'k--')
legend('\beta_t=0','\beta_t=2\times 10^{-3}','\beta_t=3\times 10^{-3}','FBP','ST 121')
title('ASC'),xlabel('\beta_s'),ylabel('Az')