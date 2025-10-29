function ncat16gtew_ideal_re(g)

g=g+1;

%(1.1) air projection + UMass blur
load ncat16gtew_air sino
load roi
Im_ideal=mbsrem4dv2(sino(:,:,:,g),repmat(roi,[1,1,64]),16,10,0,0,0);
Im_ideal=Im_ideal*5e5/sum(Im_ideal(:));
fn=['ncatlestew_Im_ideal_g' num2str(g) '.mat'];
save(fn,'Im_ideal')

% compose 16 gate together
% temp=[];
% for g=1:16
%     fn=['ncat16gtew_Im_ideal_g' num2str(g) '.mat'];
%     load(fn,'Im_ideal')
%     temp=cat(4,temp,Im_ideal);
% end
% Im_ideal=temp;
% save ncat16gtew_Im_idealUm Im_ideal