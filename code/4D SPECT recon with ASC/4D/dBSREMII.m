function [y,x,cf]=dBSREMII(sino,ini_x,sMatrix,pl,sub,iter,sbeta,gbeta)
%OS version of dEM with coordinate descent.
[N,T]=size(ini_x);%N=4096, T=64
[M,h]=size(sino);%M=4096, h=3
B=64;%bin #, usually M/T/gates=64.
%split system matrix by projection angles(dynamic samples)
%'min_sys' will be used to determine the upper bound of pixel activity
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
%sub: # of subsets. sub_com: # of components in each subset. T: total views
%each subset contains 4 perpendicular views(each view has one or more projections)
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
%each row of 'index' is a subset index!

%Pre-loop parameters.
%Fixed scaling denominator
d=backward(ones(M,1),H,N,T,B);
%******************************************************
%d=invAdp2(d,pl,'t')/sub;
%***not right for dEM!!!!!!very important difference!!!
%******************************************************
%Upper bound of pixel activity
U=max(sino(:)/min_sys);%*2
x=ini_x;x(fv==0,:)=0;
y=invAdp2(x,pl);
g=forward(y,H,B,T);
cf=costfunction(g,sino);
alpha0=1;
for n=1:iter
    for s=1:sub        
        sino_sub=zeros(B,1);
        g_sub=zeros(B,1);
        y=zeros(N,T);
        for m=1:sub_com
            sino_sub=sino((index(m,s)-1)*64+1:index(m,s)*64);
            g_sub=g((index(m,s)-1)*64+1:index(m,s)*64);
            nz=find(g_sub>0);
            g_sub(nz)=sino_sub(nz)./g_sub(nz)-1;
            y(:,index(m,s))=H{index(m,s)}*g_sub;
        end
        y=invAdp2(y,pl,'t');%subset dEM gradient
        qq=zeros(N,T);qq(:,index(:,s))=d(:,index(:,s));
        qq=invAdp2(qq,pl,'t');
%         y=zeros(4096,64);
%         y=invAtran_sub(temp,index(:,s),pl,N,T); %A^-T operation
        %BSEM
        nz=find(y~=0);
        w=zeros(4096,64);
        t1=x(nz);t2=qq(nz);
        t3=zeros(length(nz),1);
        tt=find(t1<U/2);
        t3(tt)=t1(tt)./t2(tt);
        tt=find(t1>=U/2);
        t3(tt)=(U-t1(tt))./t2(tt);
        w(nz)=t3;
        y(nz)=w(nz).*y(nz);
        gamma=n*.5;%adjustable weighting parameter
        x(nz)=x(nz)+alpha0/(gamma+1)*y(nz);
        t1=x(nz);
        t1(t1<1e-10)=1e-10;
        t1(t1>(U-(1e-10)))=U-(1e-10);
        x(nz)=t1;
        y=invAdp2(x,pl);
        k=mod(s,8)+1;
        for m=1:sub_com            
            g((index(m,s)-1)*64+1:index(m,s)*64)=y(:,index(m,s))'*H{index(m,s)};%for cost function
            g((index(m,k)-1)*64+1:index(m,k)*64)=y(:,index(m,k))'*H{index(m,k)};%for update
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

function y=forward(x,H,B,T)
y=zeros(B*T,1);
for n=1:T
    y((n-1)*B+1:n*B)=x(:,n)'*H{n};
end
return

function x=backward(y,H,N,T,B)
x=zeros(N,T);
for n=1:T
    x(:,n)=H{n}*y((n-1)*B+1:n*B);
end
return