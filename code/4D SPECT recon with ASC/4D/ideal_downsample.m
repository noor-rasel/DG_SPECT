function B=ideal_downsample(temp)

order=21;
[f1,f2] = freqspace(order,'meshgrid');
Hd = ones(order);
r = sqrt(f1.^2 + f2.^2);
Hd((r>0.25)) = 0;
h = fwind1(Hd,hamming(order));
ham21_cf25=h(11,:)/sum(h(11,:));

ord=21;
% ham21_cf25=fir1(32,.25);%32-order Hamming
A=zeros(ord,1,1);
A(:,1,1)=ham21_cf25;
B=convn(temp,A,'same');
A=zeros(1,ord,1);
A(1,:,1)=ham21_cf25;
B=convn(B,A,'same');
A=zeros(1,1,ord);
A(1,1,:)=ham21_cf25;
B=convn(B,A,'same');
B=B(1:4:256,1:4:256,1:4:256);