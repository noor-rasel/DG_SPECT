function x=buttlpf(x,cf,order)
tot_c=sum(x(:));
[N,T,G]=size(x);
temp=zeros(sqrt(N));
for g=1:G
    for t=1:T
        temp(:)=x(:,t,g);
        temp=imbutt(temp,cf,order);
        x(:,t,g)=temp(:);
    end
end
x(find(x<0))=0;
x=x*tot_c/sum(x(:));