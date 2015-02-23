
function [ ElevRavgMedian ] = ElevRunningAvg_LRAUV( elevAngle, elevUnits, outPlot )

% ElevExpandingAvg_LRAUV.m
% Last modified Jan 12, 2015
% Ben Raanan



% Method 1: Determain offset by computing expanding avg.
%--------------------------------------------------------------------------
eleR = elevAngle(~isnan(elevAngle));

if strcmp(elevUnits,'rad')
 eleR = eleR.*(180/pi);
end
 
eleRsum = zeros(size(eleR));
eleRavg = NaN(size(eleR));
eleRmed = NaN(size(eleR));
dragVar = NaN(size(eleR));
for k=2:length(eleR)
    
    eleRsum(k) = eleRsum(k-1)+eleR(k);
    eleRavg(k) = eleRsum(k)/k;              % compute running mean
    eleRmed(k) = nanmedian(eleRavg(1:k));   % compute running median (of means)
    if k>=201
    dragVar(k) = nanvar(eleRavg(k-200:k));  % compute sliding variance (dragging)
    end
end
eleRavg(isnan(eleRavg)) = [];


% Calculate required sample size for population mean
%--------------------------------------------------------------------------
% http://www.r-tutor.com/elementary-statistics/interval-estimation/sampling-size-population-mean
% Necessary Sample Size = (Z-score)^² * sigma^² / (margin of error)^²
% 90% Z Score = 1.645, 95% Z Score = 1.96, 99 ZScore = 2.326
sigma = nanstd(eleRavg);    % standered deviation of means    
ZScore = 2.326;             % 99%
E = 0.02;                   % margin of error (in deg)     

n = ( ZScore*sigma/E )^2;
%--------------------------------------------------------------------------

% Median od means
ElevRavgMedian = median(eleRavg);


if outPlot==1
    figure;
    set(gcf,'Units','normalized','Position',[0.1 0.3 0.8 0.5],...
    'PaperPositionMode','auto')
    plot(eleRavg); hold on;
    plot(ones(size(eleRavg))*ElevRavgMedian,'--k','linewidth',2)
    plot(eleRmed,'linewidth',2)
    plot([n n],[min([min(dragVar) min(eleRavg)]) max(eleRavg)],'--r',...
        'linewidth',2)
    plot(dragVar,'linewidth',2)
    
    ylabel('Elevator angle (deg)','fontWeight','bold','fontsize',16);
    xlabel('Sample size','fontWeight','bold','fontsize',16);
    legend('Running mean',...
        ['Running mean meadian (' num2str(ElevRavgMedian,3) '^o)'],...
        'Running median of means',...
        ['Sample size needed for 99% CI (E=' num2str(E) '; n=' num2str(fix(n),6) ')'],...
        'Variance (200 dp window)','location','best');
    set(gca,'layer','top','fontWeight','bold','fontsize',14);
    axis tight; box on; grid on;   
end


if strcmp(elevUnits,'rad')
 ElevRavgMedian = ElevRavgMedian.*(pi/180);
end

% Method 2: Get elev angle vals for apex and nadir
%--------------------------------------------------------------------------
% 1) vertical mode/commanded speed?
% 2) 3 consecituve < abs(p)
% 3) p<0 nadir / p>0 apex
end
