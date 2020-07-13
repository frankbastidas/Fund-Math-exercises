%% 1.2 Implementation Recurrence relation
clc;clear;
x = 5; % start point
Px=1.2;
p=1/Px;
iter=40;
xplt=1:(iter*2-1);
yplt=1:(iter*2-1);
for R = 1:2:(iter*2-1)
    if R > 1
        xplt(R-1)=x;
        yplt(R-1)=x;
    end
    xplt(R)=x;
    yplt(R)=1+(1-p)*(x); % equacao de recurenca 1+(1-p)*(x)
    x=yplt(R);
    disp(1+(1-p)*(x));   % print recurrence relation value
end

xRecurenca=-4:0.1:5;
yRecurenca=1+(1-p)*(xsen);
xfij=-3:0.4:5;
yfij=-3:0.4:5;
plot(xplt,yplt,'--rs',xRecurenca,yRecurencao,'-b',xfij,yfij,':g','LineWidth',2,...
                       'MarkerEdgeColor','k',...
                       'MarkerFaceColor','g',...
                       'MarkerSize',10)

%% 2. Gaussian elimination Ax=B
clc;clear;
mat=[2 3 -1;4 4 -3;2 -3 1]; %SELA(A)
B=[5 3 -1]';
[m,n] = size(mat);
numdet=1;
fil=0;
respdet=0;
if m~=n
   disp("Isn't squard matrix");
else
    Nmat=[mat;mat(1:n-1,:)];
    for i = 1:m
        for j = 1:n
            Nmat(i+j-1,j);
            numdet=Nmat(j+i-1,j)*numdet;
        end
        respdet=respdet+numdet;
        numdet=1;
    end
    for i = 1:m
        for j = 1:n
            Nmat(i+j-1,n+1-j);
            numdet=Nmat(j+i-1,n+1-j)*numdet;
        end
        respdet=respdet-numdet;
        numdet=1;
    end
    if numdet==0
        disp("Matrix with determinant zero");
    else
        disp("calculating Matrix..");
        Nmat=[mat,B];
        for j=1:n-1
            for i=j+1:n
                if Nmat(i,j)==0
                    continue;
                else
                    NV=Nmat(i,j)/Nmat(j,j);
                    Nmat(i,:)=Nmat(i,:)-NV*Nmat(j,:);
                end               
            end
        end
        Nmat
        Xn=ones(1,m);
        disp("solving ...");
        A=0;
        for i=m:-1:1
            if i~=m
                for j=i+1:m                    
                   A=A-Nmat(i,j)*Xn(j);
                end
            end
            Xn(i)=(A+Nmat(i,m+1))/Nmat(i,i);
            A=0;
        end
    end
end
% Solution
Xn


% Compare
% det(mat);
% respdet;


%%






