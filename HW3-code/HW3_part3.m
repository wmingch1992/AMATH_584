%Q3 part a)
clear all; close all; 

m= 100;
res = zeros(m,m);
for i=1:m
    for j=1:m
        if i>j
            A = randn(i,j);
            res(i,j) = cond(A);
        else
            res(i,j) = NaN;
        end 
    end
end

max(res)
contourf(res,[0,1,2,3,4,5,6,10,25,50,100,200])
colorbar
xlabel('n','FontSize', 20)
ylabel('m','FontSize', 20)
title('condition number for matrix(m,n), m>n','FontSize', 20)
savefig('HW3_part3_a.fig')



%%
%Q3 part b)
clear all; close all; clc
A1 = randn(50,49);
A2 = [A1 A1(:,1)];
%cond(A1)
cond(A2)
det(A2)

%%
%Q3 part c)
clear all; close all; clc
A1 = randn(50,49);
A2 = [A1 A1(:,1)];

cond(A2)

eplison = [1e-16 1e-14 1e-12 1e-10 1e-8 1e-6 1e-4];
res_cond = zeros(size(eplison,2),1);
for i=1:size(eplison,2)
    A3 = A2;
    A3(:,50) = A2(:,50) + eplison(i)*rand(50,1);
    res_cond(i) = cond(A3);
end 

res_cond

plot(eplison,res_cond)
set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')
xlabel('eplison','FontSize', 20)
ylabel('condition number','FontSize', 20)
title('condition number change with eplison','FontSize', 20)
