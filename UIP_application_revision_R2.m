clear
close all
rng(1)

data=readtable("restud_data.csv");

burn_in=152; %Initial periods for which alpha=0 (i.e., periods not used in estimation of c)
t_star=202; %Final pre-treatment period
cohort_end=273; %When benefits run out, used to define cohorts
extrapolation_end=365; %We extrapolate up to one year from initial unemployment

reform_day=min(data.begin(logical(data.after))); %The date of the reform
durations=data.dur*7; %The unemployment durations, in days
periods=linspace(1,extrapolation_end,extrapolation_end);

%% A cross-cohort study

bootstrap_replications=1000; %Number of bootstrap replications
level=0.95; %Level for the confidence bands
parallel_level=0.6; %Level of parallel trends test


%First, let's select our sample. We include individuals who are
%eligible for the PBD extension (data.t39=1) but not for the increased
%replacement rate (data.tr=0). Thus we form Eligibility a vector of dummies
%for whether an individual satisfies these criteria.
Eligibility=(1-data.tr).*data.t39;

%Next we define our treated and untreated cohorts.
 
%Untreated are those who became unemployed between 365 + t_star and cohort_end   
%days prior to the reform:
Untreated_cohort=(data.begin>reform_day-t_star-365).*(data.begin<=reform_day-cohort_end);
%Treated are those who became unemployed between t_star days prior to
%the reform and 365-cohort_end days after the reform.
Treated_cohort=(data.begin<=reform_day+365-cohort_end).*(data.begin>reform_day-t_star);

%We only wish to keep individuals who satisfy the eligibility criteria and
%who belong to one of the two cohorts above. Let's define a vector Keep
%that contains dummies for whether we keep an individual in the sample.
keep=Eligibility.*(Untreated_cohort+Treated_cohort);

%Now we define our data matrix
absorbed_time=data.dur*7; %Vector containing the earliest time each individual is observed in the absorbing state. Here just the duration measured in days.
D=Treated_cohort;
%X=D.*floor(data.begin/7)+(1-D).*floor((data.begin+365)/7); %Matrix of controls. In this setting we control for the start day of the year of unemployment spells using the covariate banalcing method.
X=D.*(data.begin)+(1-D).*(data.begin+365); %Matrix of controls. In this setting we control for the start day of the year of unemployment spells using the covariate banalcing method.


%Let's implement our sample selection (drop rows with 0 values for Keep)
absorbed_time=absorbed_time(logical(keep));
X=X(logical(keep));
D=D(logical(keep));

%Print numbers of treated and untreated individuals
sum(D)
sum(1-D)

%Let's plot the raw means:

Y=periods>absorbed_time;
y1=mean(Y(D==1,:),1);
y2=mean(Y(D==0,:),1);

figure()
hold on
 plot(periods,y1,'b','LineWidth',1.5)
 plot(periods,y2,'r','LineWidth',1.5)
 xline(periods(1,210),'-.');
  xline(periods(1,cohort_end));
 legend('Treated Cohort','Untreated','Location','southeast')
 xlabel('$t$','fontsize',12,'Interpreter','latex');
ylabel('$\bar{Y}_{k,t}$','fontsize',12,'Interpreter','latex');
set(gcf,'position',[10,10,480,360]);
set(gca,'FontSize',12);
hold off
print(strcat('UIP_fig_levels','.jpg'),'-djpeg');


%We will loop through the twp specifications.

specifications={'common dynamics','proportional hazards'};

for j=1:2

[E1_Y1,E1_Y2,tau,EY1_imputed,delta,CI_tau_uniform,CI_EY1_uniform,CI_delta,CI_tau_pointwise,CI_EY1_pointwise,H1,H2,H1_imputed,p_value]=durationDiD(absorbed_time,D,X,t_star,burn_in,extrapolation_end,specifications{j},bootstrap_replications,level,parallel_level);

%p_value from parallel trends test:
p_value

%Figure with full Hs in all periods (up to a year)

figure()
hold on
 plot(periods,H1,'b','LineWidth',1.5)
 plot(periods,H2,'r','LineWidth',1.5)
 xline(periods(1,210),'-.');
 xline(273)
 legend('Treated Cohort','Untreated','Location','southeast')
  xlabel('$t$','fontsize',12,'Interpreter','latex');
ylabel('$\Delta_{t-1}\hat{R}_{k,t}/(t-1)$','fontsize',12,'Interpreter','latex');
set(gcf,'position',[10,10,480,360]);
set(gca,'FontSize',12);
hold off
print(strcat(num2str(j),'UIP_fig_Hs','.jpg'),'-djpeg');


%Figure with the imputed Hs

figure()
hold on
plot(periods(burn_in:end),H1(burn_in:end),'b','LineWidth',1.5)
plot(periods(1,burn_in:end),H1_imputed(1,burn_in:end),'--b','LineWidth',1.5)
plot(periods(burn_in:end),H2(burn_in:end),'r','LineWidth',1.5)

 xline(periods(1,t_star),'-.','$t^*$','fontsize',12,'Interpreter','latex');
 xline(273)
legend({'Treated','Imputed','Untreated'},'Location','southwest');

xlabel('$t$','fontsize',12,'Interpreter','latex');
ylabel('$\Delta_{t-1}\hat{R}_{k,t}/(t-1)$','fontsize',12,'Interpreter','latex');
set(gcf,'position',[10,10,480,360]);
set(gca,'FontSize',12);
hold off
print(strcat(num2str(j),'UIP_fig_H0','.jpg'),'-djpeg');


%Figure with TEs

x_fill =[periods(1,t_star+1:end), flipud(periods(1,t_star+1:end)')'];     
band_fill=[CI_tau_uniform(1,t_star+1:end), flipud(CI_tau_uniform(2,t_star+1:end)')'];    
band_fillP=[CI_tau_pointwise(1,t_star+1:end), flipud(CI_tau_pointwise(2,t_star+1:end)')']; 

figure()
hold on
plot(periods(t_star+1:end),tau(t_star+1:end),'k','LineWidth',1.5)
fill(x_fill, band_fill, 1,'facecolor', [0.5,0.5,0.5], 'edgecolor', 'none', 'facealpha', 0.4);
fill(x_fill, band_fillP, 1,'facecolor', [0.2,0.2,0.2], 'edgecolor', 'none', 'facealpha', 0.4);
yline(0,'--k','LineWidth',1.5)
 xline(273)
xlabel('$t$','fontsize',12,'Interpreter','latex');
ylabel('$\tau_t$','fontsize',12,'Interpreter','latex');

lgd=legend({'Estimated Effect',strcat(num2str(100*level),'% Uniform CIs'),strcat(num2str(100*level),'% Pointwise CIs')},'fontsize',10);
lgd.Location='southeast';

set(gcf,'position',[10,10,480,360]);
set(gca,'FontSize',12);
print(strcat(num2str(j),'UIP_fig_taus','.jpg'),'-djpeg'); 


%Figure with imputed counterfactual mean outcomes and CIs

x_fill =[periods(1,t_star+1:end), flipud(periods(1,t_star+1:end)')'];     
band_fill=[CI_EY1_uniform(1,t_star+1:end), flipud(CI_EY1_uniform(2,t_star+1:end)')']; 
band_fillP=[CI_EY1_pointwise(1,t_star+1:end), flipud(CI_EY1_pointwise(2,t_star+1:end)')']; 

figure()
hold on
plot(periods(burn_in:end),1-E1_Y1(burn_in:end),'b','LineWidth',1.5)
plot(periods(1,burn_in+1:end),EY1_imputed(1,burn_in+1:end),'--b','LineWidth',1.5)
plot(periods(burn_in:end),1-E1_Y2(burn_in:end),'r','LineWidth',1.5)

fill(x_fill, band_fill, 1,'facecolor', [0.5,0.5,0.5], 'edgecolor', 'none', 'facealpha', 0.4);
fill(x_fill, band_fillP, 1,'facecolor', [0.2,0.2,0.2], 'edgecolor', 'none', 'facealpha', 0.4);

xline(periods(1,t_star),'-.','$t^*$','fontsize',12,'Interpreter','latex');
xline(273)
legend({'Treated','Imputed','Untreated',strcat(num2str(100*level),'% Uniform CIs'),strcat(num2str(100*level),'% Pointwise CIs')},'Location','southeast');
xlabel('$t$','fontsize',12,'Interpreter','latex');
ylabel('Proportion Re-employed','fontsize',12,'Interpreter','latex');
set(gcf,'position',[10,10,480,360]);
set(gca,'FontSize',12);
print(strcat(num2str(j),'UIP_fig_point','.jpg'),'-djpeg');

hold off


%Parallel trends test figure

x_fill =[periods(1,burn_in:t_star), flipud(periods(1,burn_in:t_star)')'];     
band_fill=[CI_delta(1,:), flipud(CI_delta(2,:)')']; 

figure()
hold on

plot(periods(1,burn_in:t_star),delta,'g','LineWidth',1.5)
fill(x_fill, band_fill, 1,'facecolor', [0.5,0.5,0.5], 'edgecolor', 'none', 'facealpha', 0.4);
 xline(periods(1,t_star),'-.','$t^*$','fontsize',12,'Interpreter','latex');
yl = yline(0,'--');
hold off

xlabel('$t$','fontsize',12,'Interpreter','latex');

lgd=legend({'$\hat{\delta}_t$',strcat('$',num2str(100*parallel_level),'\%$ Uniform CIs')},'Interpreter','latex','fontsize',10);
lgd.Location='northeast';

set(gcf,'position',[10,10,480,360]);
set(gca,'FontSize',12);
print(strcat(num2str(j),'UIP_fig_test','.jpg'),'-djpeg'); 

end

