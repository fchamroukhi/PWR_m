function show_PWR_results(t, y, PWR, yaxislim)

set(0,'defaultaxesfontsize',14);
%colors = {'b','g','r','c','m','k','y'};
colors = {[0.8 0 0],[0 0 0.8],[0 0.8 0],'m','c','k','y',[0.8 0 0],[0 0 0.8],[0 0.8 0],'m','c','k','y',[0.8 0 0],[0 0 0.8],[0 0.8 0]};
style =  {'r.','b.','g.','m.','c.','k.','y.','r.','b.','g.','m.','c.','k.','y.','r.','b.','g.'};

if (nargin<4)||isempty(yaxislim)
    yaxislim = [mean(y)-2*std(y), mean(y)+2*std(y)];
end

if size(t,1)~=1, t=t';end
if size(y,2)~=1, y=y';end


n= length(y);
K = size(PWR.param.betak, 2);

gammak = PWR.param.gammak; % index of the transition time points (ou t(gammak): the time points)

set(0,'defaultaxesfontsize',14);


%% time series, regressors, and segmentation
scrsz = get(0,'ScreenSize');
figure('Position',[50 scrsz(4) 560 scrsz(4)/2]);
plot(t,y,'Color',[0.5 0.5 0.5]);%black')%
for k=1:K
    model_k = PWR.stats.regressors(:,k);
    
    active_model_k = model_k(PWR.stats.klas==k); 
    active_period_model_k = t(PWR.stats.klas==k);
    
    inactive_model_k = model_k(PWR.stats.klas ~= k); 
    inactive_period_model_k = t(PWR.stats.klas ~= k); 
    if (~isempty(active_model_k))
        hold on,
        plot(inactive_period_model_k,inactive_model_k,style{k},'markersize',0.001);
        hold on,
        plot(active_period_model_k, active_model_k,'Color', colors{k},'linewidth',3.5);
    end
end
title('Time series, PWR regimes, and segmentation')
ylabel('y');
xlabel('t');
ylim(yaxislim);

%% time series, estimated regression function, and optimal segmentation
figure('Position',[scrsz(4)/1.2 scrsz(4) 560 scrsz(4)/2]);
plot(t,y,'Color',[0.5 0.5 0.5]);%black')%
gammak=[0 gammak n];
for k=1:K
    Ik = gammak(k)+1:gammak(k+1);
    segmentk= PWR.stats.mean_function(Ik);
    hold on, plot(t(Ik),segmentk,'r','linewidth',2);
end
hold on,plot([t(PWR.param.gammak);t(PWR.param.gammak)],[min(y)-3*std(y) max(y)+3*std(y)],'k--','linewidth',2);
title('Time series, PWR function, and segmentation')
ylim([240 600]);
xlabel('t');ylabel('y');
ylim(yaxislim);



