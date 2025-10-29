%validate projection and backprojection pair
%conlusion: projection and backprojection routines are equivalent to H*x
%(projection) and H'*y (backprojection)
x=magic(5);%image
H=rand(5);%projection
g1=[1 2 3 2 1];
g2=[0 1 2 1 0];
g3=[0 0 1 0 0];
G=[g1;g1;g2;g2;g3];
%projection
%way1
pro=x.*H;%weighting
for n=1:5%blurring
    temp=conv(pro(n,:),G(n,:));
    pro(n,:)=temp(3:7);
end
pro=sum(pro);%summing
% %way2:wrong
% for n=1:5
%     temp=conv(H(n,:),G(n,:));
%     pro(n,:)=temp(3:7);
% end
% pro=pro.*x;
% sum(pro)
%system matrix
sys_mat=zeros(5,25);

temp=H(:,1:3).*G(:,3:5);
temp(5,5)=0;
sys_mat(1,:)=temp(:)';

temp=H(:,1:4).*G(:,2:5);
temp(5,5)=0;
sys_mat(2,:)=temp(:)';

temp=H.*G;
sys_mat(3,:)=temp(:);

temp=H(:,2:5).*G(:,1:4);
temp1=zeros(5);
temp1(:,2:5)=temp;
sys_mat(4,:)=temp1(:)';

temp=H(:,3:5).*G(:,1:3);
temp1=zeros(5);
temp1(:,3:5)=temp;
sys_mat(5,:)=temp1(:)';%sys_mat*x(:)==pro
%backprojection
bpro=repmat(pro,5,1);%spreading
for n=1:5%blurring
    temp=conv(bpro(n,:),G(n,:));
    bpro(n,:)=temp(3:7);
end
bpro=bpro.*H%weighting
%bpro==reshape(sys_mat'*pro',5,5)
