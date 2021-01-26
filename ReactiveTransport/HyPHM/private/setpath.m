%> @file setpath.m Adds related directories to the Matlab path.

mypath = mfilename('fullpath');
mypath = mypath(1:end-length('/private/')-length(mfilename));

% intT and other tools where also the modified genpath is stored
addpath([mypath, filesep, 'tools']);

% recursively add standard directories
addpath(genpath('scripts'))
addpath(genpath('classes'))
addpath(genpath('domains'))
addpath(genpath('opt'))
addpath(genpath('symbolic'))
