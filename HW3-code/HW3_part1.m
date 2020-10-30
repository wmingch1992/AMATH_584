%HW3, Q1
clear all; close all; clc


A = magic(7);
[q1,r1]=modified_gs(A)
[q2,r2]=qrfactor(A)
I = eye(size(A));
qrerror_modified_gs = norm(q1*r1-A,inf)/norm(A,inf);
qrerror_qrfactor = norm(q2*r2-A,inf)/norm(A,inf);
ortherror_modified_gs = norm(q1'*q1-I,inf);
ortherror_qrfactor = norm(q2'*q2-I,inf);

fprintf('Well-condition matrix, ||QR-A||\n')
fprintf('modified gs method error:  %10.3e\n',qrerror_modified_gs)
fprintf('qr factor method error:  %10.3e\n',qrerror_qrfactor)
fprintf('Well-condition matrix, ||QQ-I||\n')
fprintf('modified gs method error:  %10.3e\n',ortherror_modified_gs)
fprintf('qr factor method error:  %10.3e\n',ortherror_qrfactor)
%choose a ill-condition matrix
A2 = magic(6);
cond(A2)
[q3,r3]=modified_gs(A2)
[q4,r4]=qrfactor(A2)
I = eye(size(A2));
qrerror_modified_gs = norm(q3*r3-A2,inf)/norm(A2,inf);
qrerror_qrfactor = norm(q4*r4-A2,inf)/norm(A2,inf);
ortherror_modified_gs = norm(q3'*q3-I,inf);
ortherror_qrfactor = norm(q4'*q4-I,inf);

fprintf('ill-condition matrix, ||QR-A||\n')
fprintf('modified gs method error:  %10.3e\n',qrerror_modified_gs)
fprintf('qr factor method error:  %10.3e\n',qrerror_qrfactor)
fprintf('ill-condition matrix,||QQ-I||\n')
fprintf('modified gs method error:  %10.3e\n',ortherror_modified_gs)
fprintf('qr factor method error:  %10.3e\n',ortherror_qrfactor)


function [Q,R] = qrfactor(A)

[m,n] = size(A);
Q=eye(m);
for k = 1:n
    % Find the HH reflector
    z = A(k:m,k);
    v = [ -sign(z(1))*norm(z) - z(1); -z(2:end) ];
    v = v / sqrt(v'*v);   % remoce v'*v in den
    
    % Apply the HH reflection to each column of A and Q
    for j = 1:n
        A(k:m,j) = A(k:m,j) - v*( 2*(v'*A(k:m,j)) );
    end
    for j = 1:m
        Q(k:m,j) = Q(k:m,j) - v*( 2*(v'*Q(k:m,j)) );
    end
        
end
Q = Q';
R = triu(A);  % exact triangularity
end 

function [Q,R] =  modified_gs(A)
    % Modified Gram-Schmidt.  [Q,R] = modified_gs(A);
    [n,p] = size(A);
    Q = zeros(n,p);
    R = zeros(p,p);
    for k = 1:p
        Q(:,k) = A(:,k);
        for i = 1:k-1
            R(i,k) = Q(:,i)'*Q(:,k);
            Q(:,k) = Q(:,k) - R(i,k)*Q(:,i);
        end
        R(k,k) = norm(Q(:,k))';
        Q(:,k) = Q(:,k)/R(k,k);
    end
end



