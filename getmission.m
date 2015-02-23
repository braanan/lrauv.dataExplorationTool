function [ mission ] = getmission( syslog, time, z )


if ~isempty(syslog)
    stf = strfind(syslog,'Started mission');
    sti = ~cellfun('isempty',stf);
    
    start = syslog(sti);
    
    mname = cell(size(start));
    stime = NaN(size(start));
    for k=1:numel(start)
        
        dataS = start{k};
        
        mi = strfind(dataS,'Started mission ')+length('Started mission ');
        mname{k} = dataS(mi:end);
        
        % extract date/time
        indZ = strfind(dataS,'Z'); % flag 'Z'
        
        year = str2double(dataS(1:4));
        month = str2double(dataS(6:7));
        day = str2double(dataS(9:10));
        hh = str2double(dataS(12:13));
        mm = str2double(dataS(15:16));
        ss = str2double(dataS(18:indZ-1));
        
        
        % log time data
        stime(k) = datenum(year,month,day,hh,mm,ss);
    end
    
    
    u = unique(mname);
    Xname = cell(size(z));
    X = NaN(length(u),length(z));
    m = mname(1);
    for q = 2:length(stime)
        
        ti = (time>=stime(q-1) & time<=stime(q));
        Xname(ti) = m;
        X(strcmp(u,m),ti) = z(ti);
        m = mname(q);
        
    end;
    
    
    mission.time  = time;
    mission.z     = X;
    mission.name  = Xname;
    mission.namelist = u;
    
else
    
    warning('Could not retrieve mission names: syslog is empty') 
    mission.time  = time;
    mission.z     = z;
    mission.name  = [];
    mission.namelist = [];
    
end





