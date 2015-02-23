function [ syslog ] = getsyslog( matpath )

% Read-in syslog associated with .mat file
% Last modified Jan 15, 2015
% Ben Raanan

matpath = char(matpath);
ftmp = fileparts(matpath);
if ispc 
    sl='\';
else
    sl='/';
end

% readin syslog to struct
fileID = fopen([ftmp sl 'syslog']);

if fileID~=-1
    
    formatSpec = '%s';
    txt = textscan(fileID,formatSpec,'Delimiter', '\n');
    syslog = txt{1,1};
    
    if ~isempty(txt)
        fclose(fileID);
    end
    
else
    syslog=[];
    warning([ftmp 'syslog : Invalid file identifier - check server connection or path'])
end
