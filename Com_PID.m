close all;clc;clear
global yd y timef
ts=0.001;%仿真步长（采样时间）
sys=tf(1,[0.1015,1.6215,2.52,1],'inputdelay',0.25);%构建传递函数
dsys=c2d(sys,ts,'z');%连续系统离散化，Z变换
[num,den]=tfdata(dsys,'v');%获取传递函数数据

u_1=0;u_2=0;u_3=0;%前两次输入
y_1=0;y_2=0;y_3=0;%前两次输出
e_1=0;e_2=0;%前两次误差

kp=7.5652; 
ki=3.3555;
kd=4.2322;

G=2000;%仿真步数
for k=1:1:G
timef(k)=k*ts;
yd(k)=1.0;%单位阶跃输入

y(k)=-den(2)*y_1-den(3)*y_2-den(4)*y_3+num(2)*u_1+num(3)*u_2+num(4)*u_3;%输出，Z变换
e(k)=yd(k)-y(k);
%PID control 增量式
delta_u(k)=kp*(e(k)-e_1)+ki*e(k)*ts+kd*(e(k)-2*e_1+e_2)/ts;
u(k)=u_1+delta_u(k);
%PID control 位置式
%e_s=e_1+e_2+e(k);
%u(k)=kp*e(k)+ki*e_s*ts+kd*(e(k)-e_1)/ts;

%参数更新
u_3=u_2;u_2=u_1;u_1=u(k);
y_3=y_2;y_2=y_1;y_1=y(k);
e_2=e_1;e_1=e(k);
end

%绘图
figure(1);%创建绘图窗口
plot(timef,yd,'r',timef,y,'k:','linewidth',2);%绘图参数
xlabel('Time(s)');ylabel('System Response');%坐标轴
legend('STEP','PID');%图例