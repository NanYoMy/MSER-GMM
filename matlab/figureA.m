%
% image=smallstary_one(prepath,'single',6.5,7.5,100,15);
clc;
clear;
n=20;
errorGmm=[];
errorGaussian=[];
errorSum=[];
errorPowerSum=[];
totalError=zeros(n,4,2);

for i=1:n
    error=ScriptForSingle(6.5,7.5,i+20,10)
    totalError(i,:,:)=error;
end
totalError(:,1,1);
totalError(:,2,1);
totalError(:,3,1);
totalError(:,4,1);

hold off;
figure(1);
n=20;
linewidth=5;
markersize=10;
fontsize=30;
tmp=sqrt(totalError(:,1,1).^2+totalError(:,1,2).^2);
plot([20:19+n],tmp,'-.k<','LineWidth',linewidth,'MarkerSize',markersize);
hold on
tmp=sqrt(totalError(:,2,1).^2+totalError(:,2,2).^2);

plot([20:19+n],totalError(:,2,1),'-g*','LineWidth',linewidth,'MarkerSize',markersize);
tmp=sqrt(totalError(:,3,1).^2+totalError(:,3,2).^2);
plot([20:19+n],tmp,':rd','LineWidth',linewidth,'MarkerSize',markersize);
% plot([1:50],totalError(:,4,1),'-bo');

% legend('our method ','Gaussian fitting','COG','power COG')
h=legend('our method ','Gaussian fitting','COG','fontsize',fontsize)
set(h,'Fontsize',fontsize);
xlabel('mean of Guassian noise','fontsize',fontsize);
ylabel('error(pixel)','fontsize',fontsize);
title('The error of each algorithms','fontsize',fontsize)
set(gca,'FontName','Times New Roman','FontSize',fontsize) 
% a=[21:40];
% set(gca,'xtick',a)
xlim([20 39])
% hold off;
% figure(1);
% plot([1:50],totalError(:,1,2),'-.k<');
% hold on
% plot([1:50],totalError(:,2,2),'-g*');
% plot([1:50],totalError(:,3,2),':rd');
% % plot([1:50],totalError(:,4,1),'-bo');
% % legend('our method ','Gaussian fitting','COG','power COG')
% legend('our method ','Gaussian fitting','COG')
% xlabel('mean of Guassian noise');
% ylabel('error');