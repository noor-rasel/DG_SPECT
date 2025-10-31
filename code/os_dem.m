function [y,x,cf]=os_dem(sino,ini_x,sMatrix,pl,sub,iter)
%One head system is not suitable for dEM.
%3-head
[N,T]=size(ini_x);%N=4096, T=64
[M,h]=size(sino);%M=4096, h=3
B=64;
if h==3
    sp=round(64/3);%3 heads at uniformly distributed position
    h_index=repmat((1:64)',1,3);
    h_index(:,2)=circshift(h_index(:,2),-sp);
    h_index(:,3)=circshift(h_index(:,3),-2*sp);
else
    h_index=(1:64)';
end

fv=ones(1,N)*sMatrix;
H=cell(T,1);
min_sys=1;
for n=1:T
    temp=full(sMatrix((n-1)*B+1:n*B,:)');
    H{n}=sparse(temp);%each 'H{n}' is a transpose of sub_systemmatrix
    nz=find(temp>0);
    mm=min(temp(nz));
    if mm<min_sys
        min_sys=mm;
    end
end
clear sMatrix

%subset size and orders are sensitive for OS-dEM
if sub==1
    index=1:T;index=index';
    sub_com=T;
else%sub=8;
    sub_com=T/sub;
    s_sub=reshape(1:(sub_com/4)*sub,sub_com/4,sub);
    index=zeros(sub_com/4,4,sub);
    for n=1:sub
        index(:,:,n)=repmat(s_sub(:,n),1,4)+repmat((0:3)*(sub_com/4)*sub,sub_com/4,1);
    end
    n1=1:sub/2;n2=sub/2+1:sub;
    n3=zeros(sub,1);
    for n=1:sub/2
        n3(2*n-1:2*n)=[n1(n) n2(n)];
    end
    index=index(:,:,n3);
    index=reshape(index,sub_com,sub);
end
% sub_com=T/sub;
% index=reshape(1:64,sub_com,sub);

q=backward(ones(M,1),H,N,T,B);
x=ini_x;x(fv==0,:)=0;
y=invAdp2(x,pl);
g=forward(y,H,h_index);
cf=costfunction(g,sino);
for n=1:iter
    for s=1:sub
        y=zeros(N,T);
        qq=zeros(N,T);
        for i=1:h
            sino_sub=zeros(B,1);
            g_sub=zeros(B,1);            
            qq(:,index(:,s))=qq(:,index(:,s))+...
                q(:,h_index(index(:,s),i));
            for m=1:sub_com
                sino_sub=sino((index(m,s)-1)*64+1:index(m,s)*64,i);
                g_sub=g((index(m,s)-1)*64+1:index(m,s)*64,i);
                nz=find(g_sub>0);
                g_sub(nz)=sino_sub(nz)./g_sub(nz);
                y(:,index(m,s))=y(:,index(m,s))+...
                    H{h_index(index(m,s),i)}*g_sub;
            end
        end
        y=invAdp2(y,pl,'t');
        qq=invAdp2(qq,pl,'t');
        nz1=find(qq>0);
        y(nz1)=y(nz1)./qq(nz1);
        x(nz1)=x(nz1).*y(nz1);
        y=invAdp2(x,pl);        
        %k=s-1;if k==0,k=8;end
        k=mod(s,sub)+1;
        for i=1:h
            for m=1:sub_com
                g((index(m,s)-1)*64+1:index(m,s)*64,i)=...
                    y(:,index(m,s))'*H{h_index(index(m,s),i)};%for cost function
                g((index(m,k)-1)*64+1:index(m,k)*64,i)=...
                    y(:,index(m,k))'*H{h_index(index(m,k),i)};%for update
            end
        end
    end
    nz2=find(g>0);
    cf=[cf sum(sum(sino(nz2).*log(g(nz2))-g(nz2)))];    
    subplot(2,2,2),plot(cf);%,axis([30 105 3470 3510]);
    subplot(2,2,3),imagesc(reshape(x(:,1),64,64));
    subplot(2,2,4),imagesc(reshape(x(:,64),64,64));pause(0.1);
end

function cf=costfunction(g,sino)
%pernalized likelihood
nz=find(g>0);
cf=sum(sum(sino(nz).*log(g(nz))-g(nz)));
return

function y=forward(x,H,h_index)
B=size(H{1},2);
[T,h]=size(h_index);
y=zeros(T*B,h);
for m=1:h
    for n=1:T
        y((n-1)*B+1:n*B,m)=(x(:,n)'*H{h_index(n,m)})';
    end
end
return

function x=backward(y,H,N,T,B)
x=zeros(N,T);
for n=1:T
    x(:,n)=H{n}*y((n-1)*B+1:n*B);
end
return