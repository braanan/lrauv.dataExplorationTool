function [ outfname ] = processLargeMAT_LRAUV( vehicle, year, searchterm, outmatfolder )

% This function extracts specific data fields (deignated in cell array
% varfields.mat) from large mat files located on smb://atlas.shore.mbari.org/LRAUV/.
% It interpolates variables to a fixed time grid and saves parameters of
% interest to separate mat file.
%
% INPUTS:
% vehicle: 'Tethys', 'Daphne', 'Makai'
% year:     2010 - 2014 (e.g., yr=2013;)
% searchterm:   string containing mission, log or file name (e.g., '20130920_20131007')
%               input searchterm=[] if no search term is desired (will process
%               all mat files for vehicle and year).
%
% outmatfolder: path for saving processed datasets
%               (e.g., [ '/Volumes/Passport/MBARI/' num2str(yr) '/mat/workver/' ])
%
% EXAMPLES:
% year        = 2013;
% vehicle     = 'Tethys';
% searchterm  = '20130920_20131007';     % <- this is a mission folder name
% searchterm  = '201309121813_201309140344'; % <- this is a log folder name
% searchterm  = []; % <- no search term
%
% OUTPUT:
% outfname: cell arrey of new datasets paths
%
% Calls findmat_LRAUV.m and fixTimeseries.m
% Last modified Jan 02, 2015
% Ben Raanan
%--------------------------------------------------------------------------



% Load .mat file names and paths
%--------------------------------------------------------------------------
[list,listf]  = findmat_LRAUV( vehicle, year , searchterm );

% c=load(list{:});

% load field list
%--------------------------------------------------------------------------
load varfields.mat

% to change var list edit varfields and save
% save('~/Documents/MATLAB/MBARI/LoadAndFix/ServerMatFiles/mat/varfields.mat', 'varfields')



% Extract specified varfields from .mat files and interp to fixed time-grid
%--------------------------------------------------------------------------
outfname = cell(size(list));
for k = 1:length(list)
    % index file
    fileIndex = k;
    
    % check if already exists in outmatfolder
    %----------------------------------------------------------------
    d = dir(outmatfolder);
    isub = [d(:).isdir];
    if ~any(~cellfun('isempty',regexp({d(~isub).name}, listf{fileIndex})))
        
        % prompt user if file size is over 250MB
        %------------------------------------------------------------
        clear d
        d = dir(list{fileIndex});
        fsize = d.bytes/10^6;
        if fsize>=250
            str='';
            prompt = [listf{fileIndex} ' file size is ' num2str(fsize,3) 'MB. Proceed? YES/NO [YES]: '];
            while ~strcmpi(str,'NO') && ~strcmpi(str,'YES')    
                str = input(prompt,'s');
                if isempty(str)
                    str = 'YES';
                end
            end
            if strcmpi(str,'NO')
                warning([ listf{fileIndex} ' processing aborted...'])
                continue
            end
        end
        
        % gather info of mat file contents
        %------------------------------------------------------------
        varlist = whos('-file',list{fileIndex});
        varlistf = struct2cell(varlist); varlistf = varlistf(1,:)';
        
        
        % construct a variable index list for vars defiened in varfields
        %------------------------------------------------------------
        varInt = zeros(size(varfields));
        for c=1:length(varfields)
            if sum(strcmp((varfields{c}),varlistf))~=0
                
                varInt(c) = find(strcmp((varfields{c}),varlistf)==1);
            else
                continue
            end
        end; clear c
        varInt(varInt==0)=[];
        
        
        % load to struct
        if ~isempty(varInt)
            vars=load(char(list(fileIndex)),varlist(varInt).name);
        else
            continue
        end
        
        % extract cycle time-stamps and elevator servo offset angle
        %------------------------------------------------------------
        if isfield(vars,'SpeedControl')
            vars.VerticalControl.speedCmd = vars.SpeedControl.speedCmd;
        end
        if isfield(vars.ElevatorServo,'offsetAngle')
            vars.ElevatorOffsetAngle = vars.ElevatorServo.offsetAngle;
        end
        if isfield(vars,'CycleStarter')
            vars.CycleStarter = vars.CycleStarter;
        end
        vars = rmfield(vars,{'ElevatorServo','CycleStarter','SpeedControl'});
        
        
        % get parameters on same time grid:
        %------------------------------------------------------------
        tRes = 1/2.5; % sec, new time grid sampling rate (for interp)
        % tRes = 1/10; % sec, new time grid sampling rate (for interp)
        % create new time grid
        t1 = round2(min(vars.depth.time),tRes/(3600*24)); % first time
        t2 = max(vars.depth.time); % last time
        time=t1:tRes/(3600*24):t2; % new 1 sec grid
        
        % prompt user if log time span is under 30 minutes
        tspan=24*60*(time(end)-time(1));
        if tspan<30
        str='';
            prompt = [listf{fileIndex} ': log time span is only ' num2str(tspan,3) ' min long. Proceed? YES/NO [YES]: '];
            while ~strcmpi(str,'NO') && ~strcmpi(str,'YES')    
                str = input(prompt,'s');
                if isempty(str)
                    str = 'YES';
                end
            end
            if strcmpi(str,'NO')
                error([ listf{fileIndex} ' processing aborted...'])
                continue
            end
        end
        
        % list new variable names
        fixVars = vars; % make copy
        fields = fieldnames(fixVars);
        
        % in some datasets time grids vectors are not strictly monotonic increasing
        % this loop fixes time vectoer and intrepolates data to a new time grid
        % with sampling intervals specifies by tRes
        for c = 1:numel(fields)
            
            % case VerticalControl struct (to get pitchCmd info)
            if strcmp('VerticalControl', fields{c})
                
                interpVars.Cmd = ineterpCmd_LRAUV(fixVars.(fields{c}), time);
                
            else
                temp = fixVars.(fields{c}); % extract
                % fix time and interpolate to new time grid
                if numel(temp.value)>1
                    [temp.time, temp.value] = fixTimeseries(temp.time, temp.value);
                    interpVars.(fields{c}) = interp1(temp.time, temp.value, time);
                    
                else
                    warning(['Log ' char(listf(fileIndex)) ': ' fields{c} ' record missing!'])
                end
            end
        end; clear c temp
        
        % save
        %------------------------------------------------------------
        %
        outfname{k} = [ outmatfolder 'LRAUV_SIM_' char(listf(fileIndex)) ];
        save(outfname{k},'vars','fixVars','interpVars','time');
        % uisave({'vars','fixVars','interpVars','time'},outfname{k});
        clear vars fixVars interpVars time
        %}
    else
        outfname{k} = [ outmatfolder 'LRAUV_SIM_' char(listf(fileIndex)) ];
        warning(['mat file alredy exists in output folder: ' outfname{k}])
    end
end
end

