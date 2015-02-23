function [ matpath, matname ]  = findmat_LRAUV( vehicle, year, varargin )

% This function returns .mat file paths found on smb://atlas.shore.mbari.org/LRAUV/
% 
% 
% INPUTS:
% vehicle: 'Tethys', 'Daphne', 'Makai'
% year: 2010 - 2014 (e.g., yr=2013;) 
% search term: string containing mission, log or file name 
%              (e.g., '20130920_20131007')
% 
% EXAMPLE: 
% year = 2013;
% vehicle  = 'Tethys';
% searchterm  = '20130920_20131007';     % <- this is a mission folder name
% searchterm  = '201309121813_201309140344'; % <- this is a log folder name
% searchterm  = []; % <- no search term (or leave blank) 
% 
% [ matpath ]  = findmat_LRAUV( vehicle, year, searchterm )
%
% Last modified Jan 02, 2015
% Ben Raanan
%--------------------------------------------------------------------------

% check inputs
%--------------------------------------------------------------------------
assert(ischar(vehicle),'Input vehicle must be a string.')
assert(isnumeric(year),'Input year must be numeric.')
if nargin>2 
    insearch = varargin{1};
    if ~isempty(insearch)
        assert(ischar(insearch),'Search term must be a string.')
    end
end

% get .mat path and name lists 
%--------------------------------------------------------------------------
% workd = '~/Documents/MATLAB/MBARI/LoadAndFix/ServerMatFiles/';
% load([workd 'mat/LRAUVmatFiles.mat']); 

yr = ['y' num2str(year)];


% check if LRAUVmatFiles.mat exists
%---------------------------------------------------------------------
try
    load LRAUVmatFiles.mat  % generated by: matFilePaths_LRAUV.m
catch
    str='';
    prompt = 'Could not find file LRAUVmatFiles.mat. Crawl server for mat files? YES/NO [YES]:';
    while ~strcmpi(str,'NO') && ~strcmpi(str,'YES')
        str = input(prompt,'s');
        if isempty(str)
            str = 'YES';
        end
    end
    if strcmpi(str,'YES')
        
        display('matFilePaths_LRAUV is crawling server! May take a few minutes to run...')
        serverpath = uigetdir(pwd,'Select server folder: \\atlas\LRAUV' );
        matFilePaths_LRAUV( serverpath );
        load LRAUVmatFiles.mat
    elseif strcmpi(str,'NO')
        error('Breaking out of findmat_LRAUV');
    end
end



% check if LRAUVmatFiles.mat is outdated
%---------------------------------------------------------------------
d=dir(which('LRAUVmatFiles.mat'));
if (datenum(clock)-d.datenum)>30
    str='';
    prompt = ['Mat file record is ' num2str((datenum(clock)-d.datenum),2)...
        ' days old. Crawl server for mat files? YES/NO [YES]:'];
    while ~strcmpi(str,'NO') && ~strcmpi(str,'YES')
        str = input(prompt,'s');
        if isempty(str)
            str = 'YES';
        end
    end
    if strcmpi(str,'YES')
        display('matFilePaths_LRAUV is crawling server! May take a few minutes to run...')
        serverpath = uigetdir(pwd,'Select server folder: atlas/LRAUV' );
        matFilePaths_LRAUV( serverpath );
        load LRAUVmatFiles.mat
    elseif strcmpi(str,'NO')
        load LRAUVmatFiles.mat
    end
end




list = LRAUVmatFiles.(vehicle).matFilePaths.(yr);
listf= LRAUVmatFiles.(vehicle).matFileNames.(yr);

% find .mat files and log paths
%--------------------------------------------------------------------------
if nargin>2 && ~isempty(insearch)
    f=~cellfun('isempty',(strfind(list,insearch)));
    [f,~] = find(f);
    if any(f)
        matpath = list(f,2);
        matname = listf(f);
    else
        error(['Your search - ' insearch ' - did not match any .mat files or paths'])
    end
else
    matpath = list(:,2);
    matname = listf;
end
end