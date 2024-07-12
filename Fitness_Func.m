function BsJ=pid_fm_def(kx,BsJ)
global yd y timef
ts=0.001;%仿真步长（采样时间）
sys=tf(1,[0.1015,1.6215,2.52,1],'inputdelay',0.25);%构建传递函数
dsys=c2d(sys,ts,'z');%连续系统离散化，Z变换
[num,den]=tfdata(dsys,'v');%获取传递函数数据

%初始化
u_1=0;u_2=0;u_3=0;%前三次输入
y_1=0;y_2=0;y_3=0;%前三次输出
e_1=0;e_2=0;%前两次误差
B=0;%适应度值

G=2000;%仿真步数
for k=1:1:G
timef(k)=k*ts;
yd(k)=1.0;%输入

y(k)=-den(2)*y_1-den(3)*y_2-den(4)*y_3+num(2)*u_1+num(3)*u_2+num(4)*u_3;%输出，Z变换
e(k)=yd(k)-y(k);
ey(k)=y(k)-y_1;
kp=kx(1);ki=kx(2);kd=kx(3);
%PID control 增量式
delta_u(k)=kp*(e(k)-e_1)+ki*e(k)*ts+kd*(e(k)-2*e_1+e_2)/ts;
u(k)=u_1+delta_u(k);

%参数更新
u_3=u_2;u_2=u_1;u_1=u(k);
y_3=y_2;y_2=y_1;y_1=y(k);
e_2=e_1;e_1=e(k);
end
for i=1:1:G
        Ji(i)=0.999*i*abs(e(i))+0.001*u(i)^2; %适应度函数
        B=B+Ji(i);
      %%%平方误差积分（ISE）
%       Ji(i)=e(i)^2; 
%       B=B+Ji(i);
      %%%绝对误差积分（IAE）
%       Ji(i)=abs(e(i)); 
%       B=B+Ji(i);
      %%%时间乘平方误差积分（ITSE）
%       Ji(i)=i*e(i)^2; 
%       B=B+Ji(i);
      %%%时间乘绝对误差积分（ITAE）
%       Ji(i)=i*abs(e(i)); 
%       B=B+Ji(i);     
       if e(i)<0    %惩罚项
             B=B+100*abs(e(i))+10*ey(i);
       end
end
BsJ=B;