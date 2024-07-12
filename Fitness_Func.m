function BsJ=pid_fm_def(kx,BsJ)
global yd y timef
ts=0.001;%���沽��������ʱ�䣩
sys=tf(1,[0.1015,1.6215,2.52,1],'inputdelay',0.25);%�������ݺ���
dsys=c2d(sys,ts,'z');%����ϵͳ��ɢ����Z�任
[num,den]=tfdata(dsys,'v');%��ȡ���ݺ�������

%��ʼ��
u_1=0;u_2=0;u_3=0;%ǰ��������
y_1=0;y_2=0;y_3=0;%ǰ�������
e_1=0;e_2=0;%ǰ�������
B=0;%��Ӧ��ֵ

G=2000;%���沽��
for k=1:1:G
timef(k)=k*ts;
yd(k)=1.0;%����

y(k)=-den(2)*y_1-den(3)*y_2-den(4)*y_3+num(2)*u_1+num(3)*u_2+num(4)*u_3;%�����Z�任
e(k)=yd(k)-y(k);
ey(k)=y(k)-y_1;
kp=kx(1);ki=kx(2);kd=kx(3);
%PID control ����ʽ
delta_u(k)=kp*(e(k)-e_1)+ki*e(k)*ts+kd*(e(k)-2*e_1+e_2)/ts;
u(k)=u_1+delta_u(k);

%��������
u_3=u_2;u_2=u_1;u_1=u(k);
y_3=y_2;y_2=y_1;y_1=y(k);
e_2=e_1;e_1=e(k);
end
for i=1:1:G
        Ji(i)=0.999*i*abs(e(i))+0.001*u(i)^2; %��Ӧ�Ⱥ���
        B=B+Ji(i);
      %%%ƽ�������֣�ISE��
%       Ji(i)=e(i)^2; 
%       B=B+Ji(i);
      %%%���������֣�IAE��
%       Ji(i)=abs(e(i)); 
%       B=B+Ji(i);
      %%%ʱ���ƽ�������֣�ITSE��
%       Ji(i)=i*e(i)^2; 
%       B=B+Ji(i);
      %%%ʱ��˾��������֣�ITAE��
%       Ji(i)=i*abs(e(i)); 
%       B=B+Ji(i);     
       if e(i)<0    %�ͷ���
             B=B+100*abs(e(i))+10*ey(i);
       end
end
BsJ=B;