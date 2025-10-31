w_st121=[0.5 0.25 zeros(1,13) 0.25];
for n=1:16
    wm_st121(:,n)=circshift(w_st121',n-1);
end

noise_num=12;
filename={['D:\imagereconstruction\4D\data_16g_fbp\Im_fbp_n' num2str(noise_num) '.mat'];...
    ['D:\imagereconstruction\4D\data_16g_5st0\Im_ncat16g_nAS_s4t0_n' num2str(noise_num) '.mat'];...
    ['D:\imagereconstruction\4D\data_16g_5st0\Im_ncat16g_A_s4t0_n' num2str(noise_num) '.mat'];...
    ['D:\imagereconstruction\4D\data_16g_5st0\Im_ncat16g_AS_s4t0_n' num2str(noise_num) '.mat'];...
    ['D:\imagereconstruction\4D\data_16g_5st2\Im_ncat16g_nAS_s2t2_n' num2str(noise_num) '.mat'];...    
    ['D:\imagereconstruction\4D\data_16g_5st2\Im_ncat16g_A_s2t2_n' num2str(noise_num) '.mat'];...    
    ['D:\imagereconstruction\4D\data_16g_5st2\Im_ncat16g_AS_s2t2_n' num2str(noise_num) '.mat']};

load(filename{1});
Im_st121=reshape(Im_fbp,[30*28*20 16])*wm_st121;
Im_st121=reshape(Im_st121,[30 28 20 16]);
Im_part=single(Im_st121);
for n=1:7
    if n==1
        dsp(Im_part(:,:,4:12,1),1)
    else
        load(filename{n})
        dsp(Im_maps(:,:,4:12,1),1)
    end
end