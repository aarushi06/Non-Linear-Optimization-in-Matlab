%Kuhn-Tucker Method
%Aarushi Mehra
%Student no 82519695
clear;
x_i = ones(26,1)*5;
[d_x, objctive] = kktsystem(x_i)
i = 1;
continue_process = true;
while continue_process 
    if any(abs(d_x) > 0.01) 
        fprintf('Iteration: %d', i)
        xnew = x_i + d_x;
        [change_in_x, objective] = kktsystem(xnew)
        x_i = xnew;
        d_x = change_in_x;
        clear change_in_x;
        clear objective;
        i = i+1;
    else
        continue_process = false
    end 
end  

%double the quantityyy!!!!

function [delta_x, f] = kktsystem(x_values)
clear delta_x;
clear f;
clear x;
x = sym('x',[26,1]);
%setting price functions
u1 = 5 + 0.5*x(1)+ 0.2*x(1)^2;
u2 = 15 + 1.5*x(2) + 0.2*x(2)^2;
u3 = 6 + 0.6*x(3) + 0.2*x(3)^2;
u4 = 7 + 0.7*x(4) + 0.2*x(4)^2;
u5 = 6 + 0.6*x(5) + 0.2*x(5)^2;

%incremental cost function
Cost = int(u1,x(1)) + int(u2,x(2)) + int(u3,x(3)) + int(u4,x(4)) + int(u5,x(5)) ;

%constraints
quantity = x(26)*(150 - x(1)-x(2)-x(3)-x(4)-x(5));
c1min = (-x(1)+x(6)^2)*x(16); 
c2min = (-x(2)+x(7)^2)*x(17);
c3min = (-x(3)+x(8)^2)*x(18);
c4min = (-x(4)+x(9)^2)*x(19);
c5min = (-x(5)+x(10)^2)*x(20);

c1max = (x(1)-20+x(11)^2)*x(21);
c2max = (x(2)-40+x(12)^2)*x(22);
c3max = (x(3)-30+x(13)^2)*x(23);
c4max = (x(4)-80+x(14)^2)*x(24);
c5max = (x(5)-60+x(15)^2)*x(25);

L = Cost + quantity + c1min + c2min + c3min + c4min + c5min + c1max + c2max+ c3max + c4max + c5max;
f = [diff(L, x(1)) ; diff(L, x(2)) ; diff(L, x(3)) ; diff(L, x(4)) ; diff(L, x(5)) ; diff(L, x(6)) ; diff(L, x(7)) ; diff(L, x(8)) ; diff(L, x(9)) ; diff(L, x(10)) ; diff(L, x(11)) ; diff(L, x(12)) ; diff(L, x(13)) ; diff(L, x(14)) ; diff(L, x(15)) ; diff(L, x(16)) ; diff(L, x(17)); diff(L,x(18)); diff(L,x(19)); diff(L,x(20)); diff(L,x(21)); diff(L,x(22)); diff(L,x(23)); diff(L,x(24)); diff(L,x(25)); diff(L,x(26))];
H = hessian(L, x);

inverse_H = inv(H);
delta_x = -inverse_H*f ;
delta_x = vpa(subs(delta_x, x, x_values));
f = vpa(subs(f, x, x_values));
end




