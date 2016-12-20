theta=zeros(1,56);
theta=0.1*ones(1,56);
theta(19:20)=0.5;
n=200;
X=zeros(n,30);
for i=1:n
    for k=1:30
        X(i,k)=binornd(1,0.5);
        if X(i,k)==0
            X(i,k)=-1;
        end;
    end;
    for j=1:4000
        temp1=1;
        for k=2:29
            temp1=temp1*exp(2*theta(2*k-3)*X(i,k));
        end;
        P1=temp1/(temp1+1);
        X(i,1)=binornd(1,P1);
        if X(i,1)==0
            X(i,1)=-1;
        end;
        for k=2:29
            temp=exp(2*(theta(2*k-3)*X(i,1)+theta(2*k-2)*X(i,30)));
            P=temp/(temp+1);
            X(i,k)=binornd(1,P);
            if X(i,k)==0
                X(i,k)=-1;
            end;
        end;
        temp2=1;
        for k=2:29
            temp2=temp2*exp(2*theta(2*k-2)*X(i,k));
        end;
        P2=temp2/(temp2+1);
        X(i,30)=binornd(1,P2);
        if X(i,30)==0
            X(i,30)=-1;
        end;
    end;
end;
A=zeros(29,29);
for k=1:29
    for l=1:29
        suumm=0;
        for i=1:n
            temp1=1;
            for m=2:29
                temp1=temp1*exp(2*X(i,1)*theta(2*m-3)*X(i,m));
            end;
            suumm=suumm+(2*X(i,1)*X(i,k+1))*(2*X(i,1)*X(i,l+1))*temp1/((temp1+1)^2*n);
        end;
        A(k,l)=suumm;
    end;
end;

        
AA=zeros(29,29);
for k=1:3
    for l=1:3
        suumm=0;
        for i=1:n
            ttemp=exp(2*X(i,4)*(theta(5)*X(i,1)+theta(6)*X(i,30)));
            suumm=suumm+(2*X(i,4)*X(i,k))*(2*X(i,4)*X(i,l))*ttemp/((ttemp+1)^2*n);
        end;
        AA(k,l)=suumm;
    end;
end;
for k=4:29
    for l=4:29
        suumm=0;
        for i=1:n
            ttemp=exp(2*X(i,4)*(theta(5)*X(i,1)+theta(6)*X(i,30)));
            suumm=suumm+(2*X(i,4)*X(i,k+1))*(2*X(i,4)*X(i,l+1))*ttemp/((ttemp+1)^2*n);
        end;
        AA(k,l)=suumm;
    end;
end;
B=A(1:28,1:28);
C=A(29,1:28)*(B^(-1));
d=eig(B);
answ1=d(1);
answ2=norm(C,Inf);
ff=zeros(n,29);
for i=1:n
    for j=1:29
        if 2*X(i,1)*X(i,j+1)>0
           ff(i,j)=1;
        else
           ff(i,j)=0;
        end;
    end;
end;
answw=sum(ff);
BB=AA([1,29],[1,29]);
CC=AA(2:28,[1,29])*(BB^(-1));
dd=eig(BB);
answw1=dd(1);
answw2=norm(CC,Inf);
AAA=zeros(29,29);
for k=1:9
    for l=1:9
        suumm=0;
        for i=1:n
            ttemp=exp(2*X(i,10)*(theta(19)*X(i,1)+theta(20)*X(i,30)));
            suumm=suumm+(2*X(i,10)*X(i,k))*(2*X(i,10)*X(i,l))*ttemp/((ttemp+1)^2*n);
        end;
        AAA(k,l)=suumm;
    end;
end;
for k=10:29
    for l=10:29
        suumm=0;
        for i=1:n
            ttemp=exp(2*X(i,10)*(theta(19)*X(i,1)+theta(20)*X(i,30)));
            suumm=suumm+(2*X(i,19)*X(i,k+1))*(2*X(i,19)*X(i,l+1))*ttemp/((ttemp+1)^2*n);
        end;
        AAA(k,l)=suumm;
    end;
end;
BBB=AAA([1,29],[1,29]);
CCC=AA(2:28,[1,29])*(BB^(-1));
ddd=eig(BBB);
answww1=ddd(1);
answww2=norm(CCC,Inf);
      

        