clear
clc

x1 = xlsread('width vs signal.xlsx','B2:B101');
x2 = xlsread('width vs signal.xlsx','C2:C101');
y = xlsread('width vs signal.xlsx','D2:D101');

modelfun = @(beta, t)beta(1)./(1+(beta(1)-beta(2))/beta(2)*exp(-beta(3)*t))+beta(4)
beta0=[47,0,4,0];
mdl = NonLinearModel.fit(x1,y,modelfun,beta0)

xnew =[2.1:0.01:3.8]';
[yp1,ypci]=mdl.predict(xnew,'Prediction','observation');

hold on 
plot(x1,y,'k.','markersize',8)
plot(xnew, yp1,'k','linewidth', 1)
xlim([2.1 3.8])
m=xlabel('The width of peaks (log_1_0)','fontsize',12)
n=ylabel('Signal Value','fontsize',12)
h2 = legend('Observed value','Fited curve','location','northwest');
h3 = title('The relationship of signal value and peaks width','fontsize',13)
set(m,'FontName','Times New Roman')
set(n,'FontName','Times New Roman')
set(h2,'FontName','Times New Roman')
set(h3,'FontName','Times New Roman')
