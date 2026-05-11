clear
close all
rng(1)

n=[100,500,1000,5000,10000];

c=0.5;
eff=1;

times=20;

level=0.95;
bootstrap_replications=1000;

T=100000;
time=linspace(1,times,T);


EY0s_1=0.4;
EY0s_2=0.2;

t_star=ceil(T/2);

tau=[zeros(t_star,1);ones(T-t_star,1)]*eff;

sims=1000;


eta=1+sqrt(linspace(1,times,T)/times)'-0.5*((linspace(1,times,T)/times-(1/2)).^2)';
h_1=(eta+c)/T;
h_2=(eta)/T;
h_1_treat=h_1+tau/T;

H_1=cumsum(h_1(2:end));
H_2=cumsum(h_2(2:end));

H_1_treat=cumsum(h_1_treat(2:end));


EY0_1=1-exp(-H_1)*(1-EY0s_1);
EY0_2=1-exp(-H_2)*(1-EY0s_2);

EY0_1_treat=1-exp(-H_1_treat)*(1-EY0s_1);

EY0_1=[EY0s_1;EY0_1];
EY0_2=[EY0s_2;EY0_2];
EY0_1_treat=[EY0s_1;EY0_1_treat];

%EY0_treat=1-exp(-H_1_treat)*(1-EY0s_1);

figure()
hold on

plot(linspace(1,times,T)',EY0_1_treat,'b','LineWidth',1.5)
plot(linspace(1,times,T)',EY0_2,'r','LineWidth',1.5)
plot(linspace(1,times,T)',EY0_1,'--b','LineWidth',1.5)
xl = xline((t_star)*((times-1)/(T-1))+1,'-.','$t^*$','fontsize',12,'Interpreter','latex');
hold off

%h=title(strcat('$n=$',num2str(n(l))),'Interpreter','latex');
xlabel('$t$','fontsize',12,'Interpreter','latex');
ylabel('Mean Outcome','fontsize',12,'Interpreter','latex');

lgd=legend({'$E[Y_{t,i}|G_i=1]$','$E[Y_{t,i}|G_i=2]$','$E[Y_{t,i}^{(0)}|G_i=1]$'},'Interpreter','latex','fontsize',10);
lgd.Location='southeast';
set(gcf,'position',[10,10,480,360]);
set(gca,'FontSize',12);
print(strcat('True Mean Outcomes','.jpg'),'-djpeg'); 

figure()
hold on
plot(linspace(1,times,T-1)',((T-1)/(times))*H_1_treat./(1:T-1)','b','LineWidth',1.5)
plot(linspace(1,times,T-1)',((T-1)/(times))*H_2./(1:T-1)','r','LineWidth',1.5)
plot(linspace(1,times,T-1)',((T-1)/(times))*H_1./(1:T-1)','--b','LineWidth',1.5)
xl = xline((t_star)*((times-1)/(T-1))+1,'-.','$t^*$','fontsize',12,'Interpreter','latex');
hold off


%h=title(strcat('$n=$',num2str(n)),'Interpreter','latex');
xlabel('$t$','fontsize',12,'Interpreter','latex');
ylabel('Time-Average Hazard','fontsize',12,'Interpreter','latex');

lgd=legend({'$\Delta_{t-1}R_{1,t}/(t-1)$','$\Delta_{t-1}R_{2,t}/(t-1)$','$\Delta_{t-1}R^{(0)}_{1,t}/(t-1)$'},'Interpreter','latex','fontsize',10);
lgd.Location='southeast';
set(gcf,'position',[10,10,480,360]);
set(gca,'FontSize',12);
print(strcat('True Time-Average Hazards','.jpg'),'-djpeg'); 

inc=floor(T/times);
t_star=floor(t_star/inc);

for l=1:size(n,2)

n1=n(l);
n2=n(l);

l
for m=1:sims

%Now let's draw individual-level data

y1t=zeros(times,1);
tau_true=zeros(times,1);

Y1obs=zeros(n1,times);
Y2obs=zeros(n2,times);

Y1obs(:,1)=binornd(1,EY0_1(1),n1,1);
Y2obs(:,1)=binornd(1,EY0_2(1),n2,1);

y1t(1)=EY0_1(1);
tau_true(1)=0;

for t=2:times

prob1=( EY0_1_treat((t-1)*inc) - EY0_1_treat((t-2)*inc+1) )/(1-EY0_1_treat((t-2)*inc+1));
prob2=( EY0_2((t-1)*inc) - EY0_2((t-2)*inc+1) )/(1-EY0_2((t-2)*inc+1));

Y1obs(:,t)=binornd(1,prob1,n1,1);
Y2obs(:,t)=binornd(1,prob2,n2,1);

y1t(t)=EY0_1((t-1)*inc);

tau_true(t)=EY0_1_treat((t-1)*inc)-EY0_1((t-1)*inc);
end




Y1obs=(cumsum(Y1obs,2)>0);
Y2obs=(cumsum(Y2obs,2)>0);

D=[ones(n1,1);zeros(n2,1)];
absorbed_time=[sum(Y1obs==0,2)+1;sum(Y2obs==0,2)+1];

out=durationDiD(absorbed_time,D,[],t_star,1,20,'common dynamics',bootstrap_replications,level,level);



y1=mean(Y1obs,1)';
y2=mean(Y2obs,1)';
chat_classic=mean(y1(1:t_star)-y2(1:t_star));
y1pred_classic=y2+chat_classic;
tau_classic=y1-y1pred_classic;


%bootstrap classid DiD


splits=ceil((bootstrap_replications*times*(n1+n2))/100000000); %split up the bootstrap so we don't run into memory problems

boots=ceil(bootstrap_replications/splits);

for s=1:splits

inds1=randi(n1,n1,boots);
inds2=randi(n2,n2,boots);

Y1b=Y1obs(inds1,:);
Y2b=Y2obs(inds2,:);

for b=1:boots
y1b=mean(Y1b((b-1)*n1+1:b*n1,:),1)';
y2b=mean(Y2b((b-1)*n2+1:b*n2,:),1)';
chat_classic=mean(y1b(1:t_star)-y2b(1:t_star));
y1pred_classicb(:,(s-1)*boots+b)=y2b+chat_classic;
tau_classicb(:,(s-1)*boots+b)=y1b-y1pred_classicb(:,(s-1)*boots+b);
parclassb(:,(s-1)*boots+b)=(y1b(1:t_star-1)-y2b(1:t_star-1))-(y1b(t_star)-y2b(t_star));
end
end


std_classic=sqrt(var(y1pred_classicb,[],2));
bootdiffs_classic=abs(y1pred_classicb-y1pred_classic)./std_classic;
std_classic_tau=sqrt(var(tau_classicb,[],2));
taudiffs_classic=abs(tau_classicb-tau_classic)./std_classic_tau;
std_parrclass=sqrt(var(parclassb,[],2));
parrclass=(y1(1:t_star-1)-y2(1:t_star-1))-(y1(t_star-1)-y2(t_star-1));
bootparrclass=abs(parclassb(1:t_star-2,:)-parrclass(1:t_star-2))./std_parrclass(1:t_star-2);
c_parrclass=quantile(max(bootparrclass,[],1),level,2);
CI_parrclass=[parrclass-c_parrclass*std_parrclass,parrclass+c_parrclass*std_parrclass];
c_unif_classic=quantile(max(bootdiffs_classic(t_star+1:end,:),[],1),level,2);
CI_unif_classic=[y1pred_classic-c_unif_classic*std_classic,y1pred_classic+c_unif_classic*std_classic];
c_pointwise_classic=quantile(bootdiffs_classic,level,2);
CI_pointwise_classic=[y1pred_classic-c_pointwise_classic.*std_classic,y1pred_classic+c_pointwise_classic.*std_classic];
c_unif_classic_tau=quantile(max(taudiffs_classic(t_star+1:end,:),[],1),level,2);
CI_unif_classic_tau=[tau_classic-c_unif_classic_tau*std_classic_tau,tau_classic+c_unif_classic_tau*std_classic_tau];
c_pointwise_classic_tau=quantile(taudiffs_classic,level,2);
CI_pointwise_classic_tau=[tau_classic-c_pointwise_classic_tau.*std_classic_tau,tau_classic+c_pointwise_classic_tau.*std_classic_tau];


p_values(m)=out.p_value;
CI_parrcorr(m)=all(out.CI_delta(1,2:end)<=0)*all(out.CI_delta(2,2:end)>=0);
CI_unifcorr(m)=all((out.CI_EY1_uniform(1,t_star+1:end)<=y1t(t_star+1:end)').*(out.CI_EY1_uniform(2,t_star+1:end)>=y1t(t_star+1:end)'));
CI_pointwisecorr(m)=mean((out.CI_EY1_pointwise(1,t_star+1:end)<=y1t(t_star+1:end)').*(out.CI_EY1_pointwise(2,t_star+1:end)>=y1t(t_star+1:end)'));
CI_unifcorr_tau(m)=all((out.CI_tau_uniform(1,t_star+1:end)<=tau_true(t_star+1:end)').*(out.CI_tau_uniform(2,t_star+1:end)>=tau_true(t_star+1:end)'));
CI_pointwisecorr_tau(m)=mean((out.CI_tau_pointwise(1,t_star+1:end)<=tau_true(t_star+1:end)').*(out.CI_tau_pointwise(2,t_star+1:end)>=tau_true(t_star+1:end)'));


CI_parrcorrclass(m)=all((abs(parrclass(1:t_star-1))./std_parrclass(1:t_star-1))<=c_parrclass);
CI_unifcorr_classic(m)=all((abs(y1pred_classic(t_star+1:end)-y1t(t_star+1:end))./std_classic(t_star+1:end))<=c_unif_classic);
CI_pointwisecorr_classic(m)=mean((abs(y1pred_classic(t_star+1:end)-y1t(t_star+1:end))./std_classic(t_star+1:end))<=c_pointwise_classic(t_star+1:end),1);
CI_unifcorr_classic_tau(m)=all((abs(tau_classic(t_star+1:end)-tau_true(t_star+1:end))./std_classic_tau(t_star+1:end))<=c_unif_classic_tau);
CI_pointwisecorr_classic_tau(m)=mean((abs(tau_classic(t_star+1:end)-tau_true(t_star+1:end))./std_classic_tau(t_star+1:end))<=c_pointwise_classic_tau(t_star+1:end),1);


y1preded(:,m)=out.EY1_imputed(t_star+1:end)';
taued(:,m)=out.tau(t_star+1:end)';
y1pred_classiced(:,m)=y1pred_classic(t_star+1:end);
tau_classiced(:,m)=tau_classic(t_star+1:end);
y1true(:,m)=y1t(t_star+1:end);
true_tau(:,m)=tau_true(t_star+1:end);
end

bias(l)=mean(abs(mean(y1preded-y1true,2)),1);
MSE(l)=mean(mean((y1preded-y1true).^2,2),1);

parr_rej(l)=1-mean(CI_parrcorr);
parr_rejclass(l)=1-mean(CI_parrcorrclass);

coverage_unif(l)=mean(CI_unifcorr);
coverage_point(l)=mean(CI_pointwisecorr);

bias_classic(l)=mean(abs(mean(y1pred_classiced-y1true,2)),1);
MSE_classic(l)=mean(mean((y1pred_classiced-y1true).^2,2),1);

coverage_unif_classic(l)=mean(CI_unifcorr_classic);
coverage_point_classic(l)=mean(CI_pointwisecorr_classic);

bias_tau(l)=mean(abs(mean(taued-true_tau,2)),1);
MSE_tau(l)=mean(mean((taued-true_tau).^2,2),1);

coverage_unif_tau(l)=mean(CI_unifcorr_tau);
coverage_point_tau(l)=mean(CI_pointwisecorr_tau);

bias_classic_tau(l)=mean(abs(mean(tau_classiced-true_tau,2)),1);
MSE_classic_tau(l)=mean(mean((tau_classiced-true_tau).^2,2),1);

coverage_unif_classic_tau(l)=mean(CI_unifcorr_classic_tau);
coverage_point_classic_tau(l)=mean(CI_pointwisecorr_classic_tau);


%Figure with prallel trends test

x_fill =[(2:t_star), flipud((2:t_star)')'];      
band_fill=[out.CI_delta(1,2:end), flipud(out.CI_delta(2,2:end)')']; 

figure()
hold on
plot((2:t_star)',out.delta(2:end),'g','LineWidth',1.5)
fill(x_fill, band_fill, 1,'facecolor', [0.5,0.5,0.5], 'edgecolor', 'none', 'facealpha', 0.4);

xl = xline(t_star,'-.','$t^*$','fontsize',12,'Interpreter','latex');

yl = yline(0,'--');
hold off


h=title(strcat('$n=$',num2str(n(l)),'$, T=$',num2str(times)),'Interpreter','latex');
xlabel('$t$','fontsize',12,'Interpreter','latex');

lgd=legend({'$\hat{\delta}_t$','95% Uniform CIs'},'Interpreter','latex','fontsize',10);
lgd.Location='southeast';

set(gcf,'position',[10,10,480,360]);
set(gca,'FontSize',12);


print(strcat('Parrallel Trends Test',num2str(n(l)),num2str(times),'.jpg'),'-djpeg'); 

%figure with uniform bands

figure()
hold on
plot((1:times),out.H1,'b','LineWidth',1.5)
plot((1:times),out.H2,'r','LineWidth',1.5)
plot((1:times),out.H1_imputed,'--b','LineWidth',1.5)

xl = xline(t_star,'-.','$t^*$','fontsize',12,'Interpreter','latex');

yl = yline(0,'--');
hold off


h=title(strcat('$n=$',num2str(n(l)),'$, T=$',num2str(times)),'Interpreter','latex');
xlabel('$t$','fontsize',12,'Interpreter','latex');
ylabel('$\Delta_{t-1}\hat{R}_{k,t}/(t-1)$','fontsize',12,'Interpreter','latex');

lgd=legend({'Treated','Untreated','Imputed'},'Interpreter','latex','fontsize',10);
lgd.Location='southeast';

set(gcf,'position',[10,10,480,360]);
set(gca,'FontSize',12);


print(strcat('Estimated Time-Average Hazards',num2str(n(l)),num2str(times),'.jpg'),'-djpeg'); 

%figure with uniform bands
x_fill =[(t_star+1:times), flipud((t_star+1:times)')'];      
band_fill=[out.CI_tau_uniform(1,t_star+1:end), flipud(out.CI_tau_uniform(2,t_star+1:end)')']; 


figure()
hold on
plot((1:times),out.tau,'k','LineWidth',1.5)
plot((1:times)',tau_true,'--k','LineWidth',1.5)
fill(x_fill, band_fill, 1,'facecolor', [0.5,0.5,0.5], 'edgecolor', 'none', 'facealpha', 0.4);

xl = xline(t_star,'-.','$t^*$','fontsize',12,'Interpreter','latex');
hold off

h=title(strcat('$n=$',num2str(n(l)),'$, T=$',num2str(times)),'Interpreter','latex');
xlabel('$t$','fontsize',12,'Interpreter','latex');
ylabel('$\tau_t$','fontsize',12,'Interpreter','latex');

lgd=legend({'Estimated Effect','True Effect','95% Uniform CIs'},'fontsize',10);
lgd.Location='northwest';

set(gcf,'position',[10,10,480,360]);
set(gca,'FontSize',12);


print(strcat('DiDsimunif_tau',num2str(n(l)),num2str(times),'.jpg'),'-djpeg'); 


%figure with uniform bands
x_fill =[(t_star+1:times), flipud((t_star+1:times)')'];      
band_fill=[CI_unif_classic_tau(t_star+1:end,1)', flipud(CI_unif_classic_tau(t_star+1:end,2))']; 

figure()
hold on
plot((1:times)',tau_classic,'k','LineWidth',1.5)
plot((1:times)',tau_true,'--k','LineWidth',1.5)
fill(x_fill, band_fill, 1,'facecolor', [0.5,0.5,0.5], 'edgecolor', 'none', 'facealpha', 0.4);

xl = xline(t_star,'-.','$t^*$','fontsize',12,'Interpreter','latex');
hold off

h=title(strcat('$n=$',num2str(n(l)),'$, T=$',num2str(times)),'Interpreter','latex');
xlabel('$t$','fontsize',12,'Interpreter','latex');
ylabel('$\tau_t$','fontsize',12,'Interpreter','latex');

lgd=legend({'Estimated Effect','True Effect','95% Uniform CIs'},'fontsize',10);
lgd.Location='southwest';

set(gcf,'position',[10,10,480,360]);
set(gca,'FontSize',12);

print(strcat('DiDsimunif_classic_tau',num2str(n(l)),num2str(times),'.jpg'),'-djpeg'); 

end
Varnames={'n','Mean Absolute Bias', 'MSE','Uniform Band Coverage','Pointwise Band Coverage','Mean Absolute Bias: Classic DiD','MSE: classic Did', 'Uniform Band Coverage: Classic DiD','Pointwise Band Coverage: Classic DiD'};
results_table=table(n',bias',MSE',coverage_unif',coverage_point',bias_classic',MSE_classic',coverage_unif_classic',coverage_point_classic','VariableNames',Varnames);
writetable(results_table,'MCresults_pred.xlsx');

Varnames={'n','Mean Absolute Bias', 'MSE','Uniform Band Coverage','Pointwise Band Coverage','Mean Absolute Bias: Classic DiD','MSE: classic Did', 'Uniform Band Coverage: Classic DiD','Pointwise Band Coverage: Classic DiD','Duration Parallel Trend Rejection','Classic Parallel Trend Rejection'};
results_table=table(n',bias_tau',MSE_tau',coverage_unif_tau',coverage_point_tau',bias_classic_tau',MSE_classic_tau',coverage_unif_classic_tau',coverage_point_classic_tau',parr_rej',parr_rejclass','VariableNames',Varnames);
writetable(results_table,'MCresults.xlsx');