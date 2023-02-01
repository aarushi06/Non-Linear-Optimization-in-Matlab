%Interior Point Method
%Aarushi Mehra
%Student no 82519695
clear;
x_i = ones(16,1)*9;
mu = 1000;
[delta_x, cost, L] = ipmethod(x_i,mu)
i = 1;
alpha = 1;
continue_process = true;
while continue_process 
    if any(abs(delta_x) > 0.0001) 
        fprintf('Iteration: %d', i)
        xnew = x_i + alpha*delta_x  
        mu = mu/10000
        [change_in_x, cost, L] = ipmethod(xnew, mu)
        x_i = xnew;
        delta_x = change_in_x;
        clear change_in_x;
        clear objective;
        i = i+1;
    else
        continue_process = false
    end 
end  


function [d_x, Cost, L] = ipmethod(x_values, mu)
% total num of variables 
n = 16;
x = sym('x',[16,1]);

%setting price functions
u1 = 5 + 0.5*x(1)+ 0.2*x(1)^2;
u2 = 15 + 1.5*x(2) + 0.2*x(2)^2;
u3 = 6 + 0.6*x(3) + 0.2*x(3)^2;
u4 = 7 + 0.7*x(4) + 0.2*x(4)^2;
u5 = 6 + 0.6*x(5) + 0.2*x(5)^2;

%incremental cost function to be minimized
Cost = int(u1,x(1)) + int(u2,x(2)) + int(u3,x(3)) + int(u4,x(4)) + int(u5,x(5));
% Costx = [gradient(Cost, x)];
% Costxx = hessian(Cost, x);
%constraints
quantity = (150 - x(1)-x(2)-x(3)-x(4)-x(5))*x(16);
c1min = (x(1));%*x(16);;
c2min = (x(2));%;*x(17);;
c3min = (x(3));%*x(18);;
c4min = (x(4));%*x(19);
c5min = (x(5));%*x(20);

c1max = (x(1)-20+x(6)^2)*x(11);
c2max = (x(2)-40+x(7)^2)*x(12);
c3max = (x(3)-30+x(8)^2)*x(13);
c4max = (x(4)-80+x(9)^2)*x(14);
c5max = (x(5)-60+x(10)^2)*x(15);

lamda = 0.5;
z = mu ./ x';
X = diag(x);
Z = diag(z);
e = ones(n,1);

L = Cost + mu*(log(c1min*c2min*c3min*c4min*c5min)) + c1max+ c2max+ c3max+ c4max+ c5max+ quantity;
f = [diff(L, x(1)) ; diff(L, x(2)) ; diff(L, x(3)) ; diff(L, x(4)) ; diff(L, x(5)) ; diff(L, x(6)) ; diff(L, x(7)) ; diff(L, x(8)) ; diff(L, x(9)) ; diff(L, x(10)) ; diff(L, x(11)) ; diff(L, x(12)) ; diff(L, x(13)) ; diff(L, x(14)) ; diff(L, x(15)) ; diff(L, x(16))];
f_diff = hessian(L, x);
f_diff = vpa(subs(f_diff, x, x_values));
f = vpa(subs(f, x, x_values));
d_x = -inv(f_diff)*f;

Cost = vpa(subs(Cost, x, x_values));
L = vpa(subs(L, x, x_values));
end

