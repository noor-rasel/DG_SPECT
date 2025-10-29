

%clincal
load myo_template
myo_ind=find(myo_temp>0);

p_slice=18;M=20;G=16;
b_slice=8;
load ncatlesktew_Im_idealUm
Im_ideal=Im_ideal*8e6/sum(Im_ideal(:));
Im_part=Im_ideal(23:52,16:43,29:48,:);
immax=max(Im_part(myo_ind));
Im_part(Im_part>immax)=immax;
Im_part=Im_part/immax;
for n=1:G
    temp=Im_part(:,:,:,n);
    temp=cardiac3drot(temp);
    id_short(:,:,n)=temp(:,:,p_slice);
end
Im_comb=zeros(M,M,G,b_slice);
Im_comb(:,:,:,1)=id_short(7:26,:,:);

w_st121=[0.5 0.25 zeros(1,13) 0.25];
for n=1:G
    wm_st121(:,n)=circshift(w_st121',n-1);
end

noise_num=1;
filename={['D:\imagereconstruction\4D\data_lesk_fbp\Im_fbp_n' num2str(noise_num) '.mat'];...
    ['D:\imagereconstruction\4D\data_lesk_5st0\Im_ncatlesk_nAS_s3t0_n' num2str(noise_num) '.mat'];...
    ['D:\imagereconstruction\4D\data_lesk_5st0\Im_ncatlesk_A_s4t0_n' num2str(noise_num) '.mat'];...
    ['D:\imagereconstruction\4D\data_lesk_5st0\Im_ncatlesk_AS_s4t0_n' num2str(noise_num) '.mat'];...
    ['D:\imagereconstruction\4D\data_lesk_5st2\Im_ncatlesk_nAS_s1t2_n' num2str(noise_num) '.mat'];...
    ['D:\imagereconstruction\4D\data_lesk_5st2\Im_ncatlesk_A_s1t2_n' num2str(noise_num) '.mat'];...
    ['D:\imagereconstruction\4D\data_lesk_5st2\Im_ncatlesk_AS_s1t2_n' num2str(noise_num) '.mat']};
for k=1:7
    load(filename{k})
    if k>1
        Im_part=Im_maps(:,:,:,:);
        immax=max(Im_part(myo_ind));
        Im_part(Im_part>immax)=immax;
        Im_part=Im_part/immax;
        for n=1:G
            temp=Im_part(:,:,:,n);
            temp=cardiac3drot(temp);
            id_short(:,:,n)=temp(:,:,p_slice);
        end
        Im_comb(:,:,:,k+1)=id_short(7:26,:,:);
    else
        Im_st121=reshape(Im_fbp,[30*28*20 16])*wm_st121;
        Im_st121=reshape(Im_st121,[30 28 20 16]);
        Im_part=Im_st121(:,:,:,:);
        immax=max(Im_part(myo_ind));
        Im_part(Im_part>immax)=immax;
        Im_part=Im_part/immax;
        for n=1:G
            temp=Im_part(:,:,:,n);
            temp=cardiac3drot(temp);
            id_short(:,:,n)=temp(:,:,p_slice);
        end
        Im_comb(:,:,:,k+1)=id_short(7:26,:,:);
    end
end

Im_combs=zeros((M+2)*2,(M+2)*4,G);
im_index=[1 2 3 6 4 7 5 8];
for g=1:G
    for n=1:2
        for k=1:4
            Im_combs((n-1)*(M+2)+1:n*M+(n-1)*2,(k-1)*(M+2)+1:k*M+(k-1)*2,g)=Im_comb(:,:,g,im_index((k-1)*2+n));
        end
    end
end

G=16
for g=1:16 
interp_im(:,:,g)=interp2(Im_combs(:,:,g),2);
end
interp_im=interp_im/max(interp_im(:));
for g=1:G
    figure('position',[100 100 800 400]),imagesc(interp_im(:,:,g),[0 1]),clinicalcolor,axis image%colormap(gray)
    axis off,im(g)=getframe;pause
    close
end
figure('position',[100 100 800 400]),movie(im,10)