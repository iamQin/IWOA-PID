close all;clc;clear
global yd y timef
ts=0.001;%���沽��������ʱ�䣩
sys=tf(1,[0.1015,1.6215,2.52,1],'inputdelay',0.25);%�������ݺ���
dsys=c2d(sys,ts,'z');%����ϵͳ��ɢ����Z�任
[num,den]=tfdata(dsys,'v');%��ȡ���ݺ�������

u_1=0;u_2=0;u_3=0;%ǰ��������
y_1=0;y_2=0;y_3=0;%ǰ�������
e_1=0;e_2=0;%ǰ�������

kp=7.5652; 
ki=3.3555;
kd=4.2322;

G=2000;%���沽��
for k=1:1:G
timef(k)=k*ts;
yd(k)=1.0;%��λ��Ծ����

y(k)=-den(2)*y_1-den(3)*y_2-den(4)*y_3+num(2)*u_1+num(3)*u_2+num(4)*u_3;%�����Z�任
e(k)=yd(k)-y(k);
%PID control ����ʽ
delta_u(k)=kp*(e(k)-e_1)+ki*e(k)*ts+kd*(e(k)-2*e_1+e_2)/ts;
u(k)=u_1+delta_u(k);
%PID control λ��ʽ
%e_s=e_1+e_2+e(k);
%u(k)=kp*e(k)+ki*e_s*ts+kd*(e(k)-e_1)/ts;

%��������
u_3=u_2;u_2=u_1;u_1=u(k);
y_3=y_2;y_2=y_1;y_1=y(k);
e_2=e_1;e_1=e(k);
end

%��ͼ
figure(1);%������ͼ����
plot(timef,yd,'r',timef,y,'k:','linewidth',2);%��ͼ����
xlabel('Time(s)');ylabel('System Response');%������
legend('STEP','PID');%ͼ��