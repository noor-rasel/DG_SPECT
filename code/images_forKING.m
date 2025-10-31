function images_forKING(con)

if nargin<1
    con=0;
end

filename={'ncat16gtew_Im_idealUm.mat';
    'D:\imagereconstruction\4D\data_16g_5st2\Im_ncat16g_nAS_s1t2_n1.mat';...
    'D:\imagereconstruction\4D\data_16g_5st2\Im_ncat16g_A_s1t2_n1.mat';...
    'D:\imagereconstruction\4D\data_16g_5st2\Im_ncat16g_AS_s1t2_n1.mat';...
    'D:\imagereconstruction\4D\optSC\Im_optSC_CF1_n2.mat';...
    'D:\imagereconstruction\4D\data_16g_fbp\Im_fbp_n1.mat'};
disp_name={'Ideal';'4D-NC';'4D-A';'4D-AS';'4D-AS2';'FBP';'FBP ST121'};
if con==0
    first_fr=5;last_fr=5+8;
    for nm_ind=1:6
        load(filename{nm_ind})
        if nm_ind==1
            Im_ideal=Im_ideal*8e6/sum(Im_ideal(:));
            Im_ideal=Im_ideal(23:52,16:43,29:48,1);
            Im_show=Im_ideal;
        elseif nm_ind==6
            Im_show=Im_fbp;
            w_st121=[0.5 0.25 zeros(1,13) 0.25];
            for g=1:16
                wm_st121(:,g)=circshift(w_st121',g-1);
            end
            Im_st121=reshape(Im_fbp,[30*28*20 16])*wm_st121;
            Im_st121=reshape(Im_st121,[30 28 20 16]);
        else
            Im_show=Im_maps;
        end
        Im_show=Im_show(:,:,first_fr:last_fr,1);
        s=size(Im_show);
        d=ceil(sqrt(s(3)));
        com_im=zeros(s(2)*d,s(1)*d);
        for n=1:s(3)
            [k,m]=ind2sub([d d],n);
            com_im((k-1)*s(2)+1:k*s(2),(m-1)*s(1)+1:m*s(1))=Im_show(:,:,n)';
        end
        if nm_ind==1
            temp=interp2(com_im,2);
            immin=min(temp(:));immax=max(temp(:));
        end
        figure,imagesc(interp2(com_im,2),[immin immax])
        axis equal,axis tight,axis off,title(disp_name{nm_ind});clinicalcolor,colorbar;
    end
    Im_show=Im_st121(:,:,first_fr:last_fr,1);
    s=size(Im_show);
    d=ceil(sqrt(s(3)));
    com_im=zeros(s(2)*d,s(1)*d);
    for n=1:s(3)
        [k,m]=ind2sub([d d],n);
        com_im((k-1)*s(2)+1:k*s(2),(m-1)*s(1)+1:m*s(1))=Im_show(:,:,n)';
    end
    figure,imagesc(interp2(com_im,2),[immin immax])
    axis equal,axis tight,axis off,title(disp_name{7});clinicalcolor,colorbar;

else
    first_fr=15;last_fr=23;
    for nm_ind=1:6
        load(filename{nm_ind})
        if nm_ind==1
            Im_ideal=Im_ideal*8e6/sum(Im_ideal(:));
            Im_ideal=Im_ideal(23:52,16:43,29:48,1);
            Im_show=Im_ideal;
        elseif nm_ind==6
            Im_show=Im_fbp;
            w_st121=[0.5 0.25 zeros(1,13) 0.25];
            for g=1:16
                wm_st121(:,g)=circshift(w_st121',g-1);
            end
            Im_st121=reshape(Im_fbp,[30*28*20 16])*wm_st121;
            Im_st121=reshape(Im_st121,[30 28 20 16]);
        else
            Im_show=Im_maps;
        end
        Im_show=cardiac3drot(Im_show(:,:,:,1));
        Im_show=Im_show(:,:,first_fr:last_fr);
        s=size(Im_show);
        d=ceil(sqrt(s(3)));
        com_im=zeros(s(1)*d,s(2)*d);
        for n=1:s(3)
            [k,m]=ind2sub([d d],n);
            com_im((k-1)*s(1)+1:k*s(1),(m-1)*s(2)+1:m*s(2))=Im_show(:,:,n);
        end
        if nm_ind==1
            temp=interp2(com_im,2);
            immin=min(temp(:));immax=max(temp(:));
        end
        figure,imagesc(interp2(com_im,2),[immin immax])
        axis equal,axis tight,axis off,title(disp_name{nm_ind});clinicalcolor,colorbar;
    end
    Im_show=cardiac3drot(Im_st121(:,:,:,1));
    Im_show=Im_show(:,:,first_fr:last_fr);
    s=size(Im_show);
    d=ceil(sqrt(s(3)));
    com_im=zeros(s(1)*d,s(2)*d);
    for n=1:s(3)
        [k,m]=ind2sub([d d],n);
        com_im((k-1)*s(1)+1:k*s(1),(m-1)*s(2)+1:m*s(2))=Im_show(:,:,n);
    end
    figure,imagesc(interp2(com_im,2),[immin immax])
    axis equal,axis tight,axis off,title(disp_name{7});clinicalcolor,colorbar;
end

%%%%%%%%%%%%%
%origin
for g=1
    filename=['C:\simind\ncat16g_act_' num2str(g) '.bin'];
    fid=fopen(filename,'rb');
    temp=fread(fid,64^3,'single');
    fclose(fid);
end
temp=reshape(temp,64,64,64);
temp=permute(temp,[2 1 3])*5e5/sum(temp(:));
dsp(temp,1,0)%min(temp(:)),max(temp(:))
axis off
temp=temp(16:43,23:52,29:48);
temp=temp(:,:,6-1-1:14-1-1);
dsp(temp,1,2)
% d=3;s(1)=28;s(2)=30;
% com_im=zeros(s(1)*d,s(2)*d);
% for n=1:9
%     [k,m]=ind2sub([d d],n);
%     com_im((k-1)*s(1)+1:k*s(1),(m-1)*s(2)+1:m*s(2))=temp(:,:,n);
% end
% figure,imagesc(interp2(com_im,2),[min(temp(:)) max(temp(:))])
    axis equal,axis tight,axis off,title('Origin');clinicalcolor,colorbar;