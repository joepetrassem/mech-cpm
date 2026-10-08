function setup_paths()
%SETUP_PATHS Add the MechoCPM source folders to the MATLAB path.
%   Call once per session (the examples and tests do this for you).
root = fileparts(mfilename('fullpath'));
addpath(genpath(fullfile(root, 'src')));
end
