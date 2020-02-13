% Plots the posterior for all estimated parameters after SMC is run
figure;
title1 = ["a"; "Th"; "lambda"];

modeltype= ["Beta Binomial type 2 functional response model";
    "Beta Binomial type 3 functional response model";
    "Beta Binomial Zero handling time type 2 function response model";
    "Beta Binomial Zero handling time type 3 function response model";
    "Binomial type 2 functional response model";
    "Binomial type 3 functional response model";
    "Binomial Zero handling time type 2 function response model";
    "Binomial Zero handling time type 3 function response model"];

index = find(all(theta(:,:,M)~=-Inf)==1);
for ind=index
    subplot(2,2,ind)
    ksdensity(exp(theta(:,ind,M))); %% SMC posterior
    title(char(title1(ind)),'FontSize',12);
end
ax=axes('Units','Normal','Position',[.085 .085 .85 .85],'Visible','off');
set(get(ax,'Title'),'Visible','on');
title(char(modeltype(Models(M))),'FontSize',8);
h=get(ax,'Title');

