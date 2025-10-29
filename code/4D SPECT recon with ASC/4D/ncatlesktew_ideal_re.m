function ncatlestew_ideal_re(g)

g=g+1;

%(1) air projection
% load ncatlestew_air sino
% load roi
% Im_ideal=mbsrem4dv2(sino(:,:,:,g),repmat(roi,[1,1,64]),16,10,0,0,0);
% Im_ideal=Im_ideal*5e5/sum(Im_ideal(:));
% fn=['ncatlestew_Im_ideal_g' num2str(g) '.mat'];
% save(fn,'Im_ideal')

% compose 16 gate together
% temp=[];
% for g=1:16
%     fn=['ncatlestew_Im_ideal_g' num2str(g) '.mat'];
%     load(fn,'Im_ideal')
%     temp=cat(4,temp,Im_ideal);
% end
% Im_ideal=temp;
% save ncatlestew_Im_ideal Im_ideal

%(1.1) air projection + UMass blur
load ncatlesktew_air sino
load roi
Im_ideal=mbsrem4dv2(sino(:,:,:,g),repmat(roi,[1,1,64]),16,10,0,0,0);
Im_ideal=Im_ideal*5e5/sum(Im_ideal(:));
fn=['ncatlesktew_Im_ideal_g' num2str(g) '.mat'];
save(fn,'Im_ideal')

% compose 16 gate together
% temp=[];
% for g=1:16
%     fn=['ncatlesktew_Im_ideal_g' num2str(g) '.mat'];
%     load(fn,'Im_ideal')
%     temp=cat(4,temp,Im_ideal);
% end
% Im_ideal=temp;
% save ncatlesktew_Im_idealUm Im_ideal

%(2) analytic
% load weight64_mn
% load gbk64
% load roi
% sino=zeros(64,64,64);
%     filename=['./motion_frame/ncat16gles_act_' num2str(g) '.bin'];
%     fid=fopen(filename,'r');
%     temp=fread(fid,64^3,'single');
%     fclose(fid);
%     temp=reshape(temp,64,64,64);
%     temp=temp.*repmat(roi,[1 1 64]);    
%     for n=1:64
%         if n<33
%             sino(:,:,n)=proj3d_sa(temp,n,wp_vray((n-1)*64+1:n*64),...
%                 wp_ipxl((n-1)*64+1:n*64),wp_wgt((n-1)*64+1:n*64),1,gb_temp);
%         else
%             m=n-32;
%             sino(:,:,n)=proj3d_sa(temp,n,wp_vray((m-1)*64+1:m*64),...
%                 wp_ipxl((m-1)*64+1:m*64),wp_wgt((m-1)*64+1:m*64),1,gb_temp);
%         end
%     end
% Im_ideal=mbsrem4dv2(sino,repmat(roi,[1,1,64]),16,10,0,0,0);
% Im_ideal=Im_ideal*5e5/sum(Im_ideal(:));
% fn=['ncatlestew_Im_ideal_Ag' num2str(g) '.mat'];
% save(fn,'Im_ideal')

% compose 16 gate together
% temp=[];
% for g=1:16
%     fn=['ncatlestew_Im_ideal_Ag' num2str(g) '.mat'];
%     load(fn,'Im_ideal')
%     temp=cat(4,temp,Im_ideal);
% end
% Im_ideal=temp;
% save ncatlestew_Im_idealA Im_ideal
%%%%%% analytic seems better!