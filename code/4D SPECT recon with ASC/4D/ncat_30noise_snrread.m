function ncat_30noise_snrread

%read SNRs
%(1)
%sbeta=[0.001:0.002:0.02]; no AC and SC
%sbeta=[0.0001:0.0002:0.002]; AC only
%sbeta=[0.0001:0.0002:0.002]; AC+SC
cat_name={'Im_les_nAS_s';'Im_ncatles_A_s';'Im_ncatles_AS_s'};
snr_name={'snr_s';'snr_sac';'snr_sasc'}; 
snr_maps=zeros(3,10,30);
for m=1:3
    %for n=1:10
        for k=1:30
            filename=[cat_name{m} num2str(10) 't0_n' num2str(k) '.mat'];
            load(filename,snr_name{m});
            snr_maps(m,:,k)=eval(snr_name{m});
        end
end
save ncatles_4D_maps10_snr snr_maps

%(2) snr for lesion
cat_name={'Im_ncatles_nAS_s';'Im_ncatles_A_s';'Im_ncatles_AS_s'};
snr_name1={'snr_s';'snr_sac';'snr_sasc'};
snr_name23={'snr_t';'snr_tac';'snr_tasc'};
snr_les=zeros(3,5,3,30);
t_index=[0 2 3];
for m=1:3
    for t=1:3
        for k=1:30
            filename=[cat_name{m} num2str(5) 't' num2str(t_index(t)) '_n' num2str(k) '.mat'];
            if t==1
                load(filename,snr_name1{m});
                snr_les(m,:,t,k)=eval(snr_name1{m});
            else
                load(filename,snr_name23{m});
                snr_les(m,:,t,k)=eval(snr_name23{m});
            end            
        end
    end
end
save ncatles_4D_snr snr_les
%(3) snr for normal spatial
cat_name={'./data_16g_5st0/Im_ncat16g_nAS_s';'./data_16g_5st0/Im_ncat16g_A_s';'./data_16g_5st0/Im_ncat16g_AS_s'};
snr_name1={'snr_s';'snr_sac';'snr_sasc'};
snr_name23={'snr_t';'snr_tac';'snr_tasc'};
snr_16g=zeros(3,5,1,30);
t_index=[0 2 3];
for m=1:3
    for t=1
        for k=1:30
            filename=[cat_name{m} num2str(5) 't' num2str(t_index(t)) '_n' num2str(k) '.mat'];
            if t==1
                load(filename,snr_name1{m});
                snr_16g(m,:,t,k)=eval(snr_name1{m});
            else
                load(filename,snr_name23{m});
                snr_16g(m,:,t,k)=eval(snr_name23{m});
            end            
        end
    end
end