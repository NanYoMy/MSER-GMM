%
% image=smallstary_one(prepath,'single',6.5,7.5,100,15);
clc;
clear;
n=20;
errorGmm=[];
errorGaussian=[];
errorSum=[];
errorPowerSum=[];
totalError=zeros(n,2,2);
for i=1:n
    error=ScriptForDouble(6.5,5.5,10,10.5,i+20,10)
    totalError(i,:,:)=error;
end

n=20;
hold off;
figure(1);
fontsize=20;
tmp=sqrt(totalError(:,1,1).^2)+(totalError(:,1,2).^2)
plot([20:n+19],tmp,'-.rd','LineWidth',3,'MarkerSize',10);
xlabel('mean of Guassian noise','fontsize',fontsize);
ylabel('error(pixel)','fontsize',fontsize);
title('error of each estimated center','fontsize',fontsize)
set(gca,'FontName','Times New Roman','FontSize',fontsize) 
xlim([20 39])

figure(2);
tmp=sqrt(totalError(:,2,1).^2)+(totalError(:,2,2).^2)
plot([1:n],tmp,'-.rd','LineWidth',3,'MarkerSize',10);
xlabel('mean of Guassian noise','fontsize',fontsize);
ylabel('error(pixel)','fontsize',fontsize);
title('error of each estimated center','fontsize',fontsize)
set(gca,'FontName','Times New Roman','FontSize',fontsize) 
xlim([20 39])