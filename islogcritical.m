function [ccomp, critical_timestamp] = islogcritical( vehicle, year, matpath )

% Locate critical faults associated with log of intrest
% Last modified Jan 15, 2015
% Ben Raanan


load syslog_AllVehicles_comp

assert(numel(matpath)==1,'Function can process only one log at a time!')
matpath = char(matpath);

% get log name from path


ftmp = fileparts(matpath);
slashi=strfind(ftmp,filesep);
log=ftmp(slashi(end)+1:end) ;


% index log in CRITICAL message list
list = syslogs.(vehicle).(['y' num2str(year)]).Fault.data.compFilt;
ind = strcmp(list(:,3),['D' log]);


if any(ind)
    % collect CRITICAL error timestamps
    critical_timestamp = zeros(sum(ind),3);
    critical_timestamp(:,1) = cell2mat(list(ind,10));
    for j=1:sum(ind)
        critical_timestamp(j,2) = critical_timestamp(j,1)-5/(24*60*60);
        critical_timestamp(j,3) = critical_timestamp(j,1)+5/(24*60*60);
    end
    
    % get CRITICAL error metadata
    if verLessThan('matlab','8.4')
        ccomp = list(ind,:);
    else
        ccomp = syslogs.(vehicle).(['y' num2str(year)]).Fault.compFiltTable(ind,:);
    end
else
    critical_timestamp=[];
    ccomp=[];
end
end
