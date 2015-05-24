clear
clc

fid = fopen('xiangtu.csv');
x = textscan(fid,'%f %s %f');
[x1,index] = xlsread('all_distribution.xls','A2:K128');
fclose(fid);
y = x{:,3};
x = x{:,2};


t = zeros(1,128);
t(1) = 1;
t(128) = length(x);

shuju = zeros(9000,127);

n = 2;
for i=2:length(x)
    if ~ strcmp(x{i},x{i-1})
        t(n) = i-1;
        n = n+1;
    end
end
n = 1;
for i = 1:127
    k = randperm(9000);
    m = y(t(n):t(n+1));
    m2 = m(k);
    shuju(:,n) = m2;
    n = n+1;
end
[a b] = sort(sum(shuju)/9000);
shuju = shuju(:,b);
index2 = index(b);

set(gca,'FontSize',4)
boxplot(shuju,'labels',index2,'outliersize',3,'labelorientation','inline')
% boxplot(shuju,'labels',index2,'orientation','horizontal','outliersize',3)
% title('The Distribution of Width','FontName','Times New Roman','fontsize',12)
% ylabel('EID of 127 Reference Human Epigenomes','FontName','Times New Roman','fontsize',12)
% xlabel('Peaks Width','FontName','Times New Roman','fontsize',12)
    
