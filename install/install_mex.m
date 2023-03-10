% INSTALL_MEX Compile all MEX files in package
%
% Usage
%    install_mex();
%
% Description
%    Compiles all MEX files that are part of the ASPIRE package. This function
%    should only be called once when ASPIRE is first installed.

function install_mex()
    current_path = pwd;

    traversedirrectorytree(aspire_root(), @run_makemex);

    cd(current_path);
end

% RUN_MAKEMEX Run the `makemex` script pointed to by `fullFileName`
%
% Usage
%    run_makemex(fullFileName);
%
% Input
%    fullFileName: A file name of a potential `makemex` script
%
% Description
%    If `fullFileName` is the path of a script with the name `makemex`, it will
%    be executed in its directory.

function run_makemex(fullFileName)
    [mexPath, fname] = fileparts(fullFileName);

    if strcmp(fname, 'makemex')
        cd(mexPath);
        fprintf('Running %s\n', fullFileName);
        run(fname);
    end
end

% TRAVERSEDIRRECTORYTREE Apply "fileOperation" to a directory tree.
%
%   TRAVERSEDIRRECTORYTREE(path,operation) traverse a directory tree whose
%   root is path and apply "fileOperation" to each node. The function handle
%   operation has the signature operation(fullFileName).
%
% Example:
%   dispfile = @(str)(disp(str));
%   traversedirrectorytree(path, dispfile);
%
% Yoel Shkolnisky, September 2013

function traversedirrectorytree(path, fileOperation)
    if ~strcmp(path(end), filesep)
        path(end+1) = filesep;
    end

    dirInfo = dir(path);
    files =~ [dirInfo.isdir];
    fileNames = {dirInfo(files).name};

    if ~isempty(fileNames)
        for i=1:length(fileNames)
            % Do whathever (show file names?)
            fileOperation(fullfile(path, fileNames{i}));
        end
    end

    % For each subdir, call itself again
    isDir = [dirInfo.isdir];
    dirNames = {dirInfo(isDir).name};
    dirNames(strcmp(dirNames, '.') | strcmp(dirNames, '..')) = [];

    for i=1:length(dirNames)
        traversedirrectorytree([path dirNames{i} filesep], fileOperation);
    end
end
