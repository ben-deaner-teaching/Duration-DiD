function out = durationDiD(absorbed_time,D,X,t_star,burn_in,extrapolation_end,spec,bootstrap_replications,level,parallel_level,sampling_weights)

    arguments
        absorbed_time %Vector of times at which the individual is first in the absorbing state
        D %Vector of treatment-group indicators (D=1 if individual is in treated group)
        X %Vector of covariate values
        t_star %Final pre-treatment time
        burn_in = 1 %Only use observations at time burn_in onwards when estimating c
        extrapolation_end = floor(max(absorbed_time))  %extrapolate only up to this period
        spec (1,:) char {mustBeMember(spec,{'common dynamics','proportional hazards'})} = 'common dynamics' %Specification
        bootstrap_replications = 1000 %Number of bootstrap replications
        level = 0.95 %Level of confidence bands
        parallel_level = 0.6 %Level of parallel/proportional trends test
        sampling_weights=ones(size(D)); %additional sampling weights if needed.
    end


    %Compute weights
    if ~isempty(X)
        weights=covariate_balancing(absorbed_time,X,D,0).*(sampling_weights/mean(sampling_weights));
    else
        weights=sampling_weights;
    end

    %Generate vector of time-periods
    periods=linspace(1,floor(extrapolation_end),floor(extrapolation_end));
    

    %generate a panel of outcomes.
    Y=periods>=absorbed_time;

    %Calculated weighted mean survivals for groups 1 and 2
    E1_Y1=mean(weights(logical(D)).*(1-Y(logical(D),:)),1);
    E1_Y2=mean(weights(logical(1-D)).*(1-Y(logical(1-D),:)),1);


    %Calculate corresponding values of \Delta_{t-1}R_{k,t}/(t-1)
    H1=-(log(E1_Y1)-log(E1_Y1(1)))./(periods-1);
    H2=-(log(E1_Y2)-log(E1_Y2(1)))./(periods-1);
    

    %Compute c, imputed \Delta_{t-1}R^{(0)}_{1,t}/(t-1), and \delta_t.
    if strcmp(spec,'common dynamics') %common dynamics case
        c=nanmean(H1(burn_in:t_star)-H2(burn_in:t_star));
        H1_imputed=c+H2;
        delta=(H1(burn_in:t_star)-H2(burn_in:t_star))-(H1(t_star)-H2(t_star));
    else %proportional hazards case
        c=nanmean(H1(burn_in:t_star)./H2(burn_in:t_star));
        H1_imputed=c*H2;
        delta=(H1(burn_in:t_star)./H2(burn_in:t_star))-(H1(t_star)./H2(t_star));
    end

    %Compute imputed counterfactual R_{1,t}, countefactual group-specific
    %mean outcomes, and ATTs
    R1_imputed=H1_imputed.*(periods-1)-log(E1_Y1(1));
    EY1_imputed=1-E1_Y1;
    EY1_imputed(t_star+1:end)=1-exp(-R1_imputed(t_star+1:end));
    tau=(1-E1_Y1)-EY1_imputed;

    %Draw Bootstrap samples:
    boot_sample=randi(size(D,1),size(D,1),bootstrap_replications);

    %Now let's bootstrap!
    for b=1:bootstrap_replications
        
        %Bootstrap samples of each variable
        absorbed_timeb=absorbed_time(boot_sample(:,b),:);
        Db=D(boot_sample(:,b),:);

        if ~isempty(X)
            Xb=X(boot_sample(:,b),:);
            weightsb=covariate_balancing(absorbed_timeb,Xb,Db,1);
        else
            weightsb=ones(size(D));
        end
        %Now let's repeat the steps above but on the bootstrap sample.
        Yb=periods>=absorbed_timeb;

        %Calculated weighted mean survivals for groups 1 and 2
        E1_Y1b=mean(weightsb(logical(Db)).*(1-Yb(logical(Db),:)),1);
        E1_Y2b=mean(weightsb(logical(1-Db)).*(1-Yb(logical(1-Db),:)),1);

        %Calculate corresponding values of \Delta_{t-1}R_{k,t}/(t-1)
        H1b=-(log(E1_Y1b)-log(E1_Y1b(1)))./(periods-1);
        H2b=-(log(E1_Y2b)-log(E1_Y2b(1)))./(periods-1);

        %Compute c, imputed \Delta_{t-1}R^{(0)}_{1,t}/(t-1), and \delta_t.
        if strcmp(spec,'common dynamics')
            cb=nanmean(H1b(burn_in:t_star)-H2b(burn_in:t_star));
            H1_imputedb=cb+H2b;
            deltab(b,:)=(H1b(burn_in:t_star)-H2b(burn_in:t_star))-(H1b(t_star)-H2b(t_star));
        else
            cb=nanmean(H1b(burn_in:t_star)./H2b(burn_in:t_star));
            H1_imputedb=cb*H2b;
            deltab(b,:)=(H1b(burn_in:t_star)./H2b(burn_in:t_star))-(H1b(t_star)./H2b(t_star));
        end

        %Compute imputed counterfactual R_{1,t}, countefactual group-specific
        %mean outcomes, and ATTs
        R1_imputedb(b,:)=H1_imputedb.*(periods-1)-log(E1_Y1b(1));
        EY1_imputedb(b,:)=1-E1_Y1b;
        EY1_imputedb(b,t_star+1:end)=1-exp(-R1_imputedb(b,t_star+1:end));
        taub(b,:)=(1-E1_Y1b)-EY1_imputedb(b,:);

    end


    %Duration DiD Standard Errors
    std_EY1=sqrt(var(EY1_imputedb,[],1));
    std_tau=sqrt(var(taub,[],1));
    std_delta=sqrt(var(deltab,[],1));

    %Form absolute deviations for constructing bootstrap confidence
    %sets
    boot_EY1=abs(EY1_imputedb-EY1_imputed)./std_EY1;
    boot_tau=abs(taub-tau)./std_tau;
    boot_delta=abs(deltab-delta)./std_delta;

    %Calculate uniform critical values
    c_EY1=quantile(max(boot_EY1(:,t_star+1:end),[],2),level,1);
    c_tau=quantile(max(boot_tau(:,t_star+1:end),[],2),level,1);
    c_delta=quantile(max(boot_delta,[],2),parallel_level,1);

    %Construct uniform confidence bands
    out.CI_EY1_uniform=[zeros(2,t_star),[EY1_imputed(:,t_star+1:end)-c_EY1*std_EY1(:,t_star+1:end);EY1_imputed(:,t_star+1:end)+c_EY1*std_EY1(:,t_star+1:end)]];
    out.CI_tau_uniform=[zeros(2,t_star),[tau(:,t_star+1:end)-c_tau*std_tau(:,t_star+1:end);tau(:,t_star+1:end)+c_tau*std_tau(:,t_star+1:end)]];
    out.CI_delta=[delta-c_delta*std_delta;delta+c_delta*std_delta];

    %Calculate p-value for parallel trends test
    out.p_value=mean(max(boot_delta,[],2)>abs(max(delta./std_delta)));

    %Calculate pointwise critical values
    c_EY1=quantile(boot_EY1,level,1);
    c_tau=quantile(boot_tau,level,1);

    %Construct pointwise confidence bands
    out.CI_EY1_pointwise=[EY1_imputed-c_EY1.*std_EY1;EY1_imputed+c_EY1.*std_EY1];
    out.CI_tau_pointwise=[tau-c_tau.*std_tau;tau+c_tau.*std_tau];

    %Output results
    out.delta=delta;
    out.H1_imputed=H1_imputed;
    out.tau=tau;
    out.EY1_imputed=EY1_imputed;
    out.E1_Y1=E1_Y1;
    out.E1_Y2=E1_Y2;
    out.H1=H1;
    out.H2=H2;
end


function [weights]=covariate_balancing(absorbed_time,X,D,silent)

    % Form covariate balancing weights
    [~, ~, inds] = unique(X, 'rows');

    treated = D == 1 & absorbed_time >= 1;
    untreated = D == 0 & absorbed_time >= 1;

    counts_treated = accumarray(inds(treated), 1, [max(inds), 1]);
    counts_untreated = accumarray(inds(untreated), 1, [max(inds), 1]);

    %Compute totals
    n_treated = sum(D == 1 & absorbed_time >= 1);
    n_untreated = sum(D == 0 & absorbed_time >= 1);

    % Weights: density ratio for untreated, 1 for treated
    weights = (counts_treated(inds) / n_treated) ./ (counts_untreated(inds) / n_untreated);
    weights(D == 1) = 1;

    % Put zero weight on cells with fewer than two untreated and upweight the
    % remaining cells accordingly
    censor=counts_untreated(inds) < 2;
    weights(censor)=0;
    weights(D==1)=weights(D==1)/mean(1-censor(D==1));
    weights(D==0)=weights(D==0)/mean(1-censor(D==0));

    %Report the number of dropped cells:
    if any((counts_untreated<2)==1)&&~silent
        warning(strcat('Dropped  ',num2str(sum(counts_untreated<2)),' cells out of  ',num2str(size(counts_untreated)),'left with',num2str(sum(1-censor(D==1))),' treated observations and',num2str(sum(1-censor(D==0))),' utreated.'))
    end

end