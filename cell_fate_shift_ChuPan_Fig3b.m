clc, clear, close all;
%deterministic model for S1&S2(NOT reduced lambda)
x=0:1:400; y=0:1:50; 
[x,y]=meshgrid(x,y);
global tau_R tau_G K_DG K_DR alpha_R alpha_G lambda n_R n_G
tau_R=0.13;  tau_G= 0.015;  
n_R=2;  K_DR=10;
% n_G=4;  K_DG=16.5;
n_G=5.5;  K_DG=25;
lambda=1;
alpha_R= (26.836+320.215/(1+(lambda/0.661)^4.09))*lambda*1.1;
alpha_G= (25.609+627.747/(1+(lambda/0.865)^4.635))*lambda*1.1;

K_DG_title=(K_DG/(alpha_G/lambda)); K_DR_title=(K_DR/(alpha_R/lambda));
u1=alpha_R*(tau_R+(1-tau_R)./(1+(y./K_DG).^n_G))-lambda*x; 
v1=alpha_G*(tau_G+(1-tau_G)./(1+(x./K_DR).^n_R))-lambda*y;
quiver(x,y,u1,v1,5,'filled'), hold on

x1=0:0.1:400; y1=0:0.1:400;
R=[]; G=[];
for i=1:4001
    R(i)=alpha_R/lambda*H_R(y1(i));
end

for j=1:4001
    G(j)=alpha_G/lambda*H_G(x1(j));
end
plot(R,y1,'g',x1,G,'r','linewidth',1.5);
xlabel('$[R]$','Interpreter','latex')
ylabel('$[G]$','Interpreter','latex','Rotation',0)


% %找交点（即saddle point(if exist))
% A=[R',y1']; B=[x1',G'];
% k=1; C=[];
% for p=1:50001
%     for q=1:50001
%         if abs(sqrt((A(p,1)-B(q,1))^2+(A(p,2)-B(q,2))^2))<0.9*10^(-3)
%             C(k,1)=A(p,1); C(k,2)=A(p,2);
%             D(k,1)=B(q,1); D(k,2)=B(q,2);
%             k=k+1;
%         end
%     end
% end


T=100; dt=0.1;
[t,y]=odesol([0,T],[114.047114528272,67.8377612852593],dt);

L=100; %求解时间长度
m=L/dt;
xt=y(1:m,1); yt=y(1:m,2);
figure(2);
plot(xt,yt); 
xlabel('$[R]$','Interpreter','latex')
ylabel('$[G]$','Interpreter','latex','Rotation',0)
axis([0,400,0,400]);

figure(3);
plot(t,y(:,1),'r',t,y(:,2),'g','linewidth',2);


DR=(y(1001,1)*0.197)^2;
DG=(y(1001,2)*0.216)^2


function [t,y]=odesol(inter,y0,h)
t(1)=inter(1); y(1,:)=y0;
n=(inter(2)-inter(1))/h;
for i=1:n
t(i+1)=t(i)+h;
y(i+1,:)=eulerstep(t(i),y(i,:),h);
end
end

function y=eulerstep(t,y,h)
%one step of the Euler Method
y=y+h*ydot(t,y);
end

function z=ydot(t,y)
% y(1)is Red,y(2)is Green.
global alpha_R alpha_G lambda
z(1)=alpha_R*H_R(y(2))-lambda*y(1);
z(2)=alpha_G*H_G(y(1))-lambda*y(2);
end

function H1=H_R(y)
global tau_R K_DG n_G
H1=tau_R+(1-tau_R)/(1+(y/K_DG)^n_G);
end

function H2=H_G(x)
global tau_G K_DR n_R
H2=tau_G+(1-tau_G)/(1+(x/K_DR)^n_R);
end


