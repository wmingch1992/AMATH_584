%Midterm AMATH584 Q5 - Mingcheng Wang

clear all; close all; clc

%test matrix
A = rand(8)  
[L,U,P] = lu_pivot(A) 

%compare the result with the lu function in matlab
[L_2,U_2,P_2] = lu(A)

%difference between my code and the lu function in matlab 
L_2-L
U_2-U
P_2-P

function [L,U,P]=lu_pivot(A)

%get size of matrix
[n,n]=size(A);

%initialize L,U and P
L=eye(n); 
U=A;
P=eye(n); 

for k=1:n
    %get the index and value of the coloumn below the diagnoal 
    [pivot m]=max(abs(U(k:n,k)));
    m=m+k-1;
    if m~=k
        % do the exchange in U
        temp=U(k,:);
        U(k,:)=U(m,:);
        U(m,:)=temp;
        % do the excahnge in P
        temp=P(k,:);
        P(k,:)=P(m,:);
        P(m,:)=temp;
        if k >= 2
            temp=L(k,1:k-1);
            L(k,1:k-1)=L(m,1:k-1);
            L(m,1:k-1)=temp;
        end
    end
    for j=k+1:n
        %get the ratio 
        L(j,k)=U(j,k)/U(k,k);
        %zero out the value 
        U(j,:)=U(j,:)-L(j,k)*U(k,:);
    end
end

end 