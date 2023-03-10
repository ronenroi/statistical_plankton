function test_cryo_normalize_background_outofcore
%
% Test the function cryo_normalized_backgorund_outofcore by comapring its
% output to cryo_normalize_background. Both functions shouls have the same
% output.
%
% Yoel Shkolnisky, May 2016.

K=10;   % Number of projection images 
n=89;   % Each projection image is of size nxn.
max_shift_2d=0;     % Not used.
shift_step_2d=1;    % and max_shift_2d with steps of  shift_step_2d pixels.
snr=1;  % Not used. Only clean proejctions are used.

% Generate proejctions. Use only the generated clean projections "projs".
initstate;
projs=cryo_gen_projections(n,K,snr,max_shift_2d,shift_step_2d);

% Write and read the stack to make sure we are using precisely the same
% projections.
instackname=tempname;
WriteMRC(single(projs),1,instackname); % Save the projections to an MRC file.
projs=ReadMRC(instackname);
projs=double(projs); % Although we read single precision numbers, cast to 
    % double precision since cryo_normalize_background_outofcore uses
    % double precision internally, and we want both
    % cryo_normalize_background and cryo_normalize_background_outofcore to
    % have exactly the same roundoff error.

% Downsample the images in-memory
projs_normalized1=cryo_normalize_background(projs);
projs_normalized1=single(projs_normalized1); % The numbers in 
    % cryo_normalize_background_outofcore are read below as single
    % precision.

% Downsample the images from an MRC file
outstackname=tempname;
cryo_normalize_background_outofcore(instackname,outstackname)
projs_normalized2=ReadMRC(outstackname);

err=norm(projs_normalized1(:)-projs_normalized2(:));  % Should be zero.
fprintf('diff between two methods = %5.3e\n',err);
if err<1.0e-10
    fprintf('test ok\n');
else
    fprintf('test FAILED\n');
end
