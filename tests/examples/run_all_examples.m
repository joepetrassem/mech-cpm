function run_all_examples(include_slow)
%RUN_ALL_EXAMPLES Run every example script in this folder.
%   RUN_ALL_EXAMPLES() runs examples 1-6 (about a minute in total).
%   RUN_ALL_EXAMPLES(true) also runs the slow example 7 (paper Figs. 4-7,
%   several minutes).
if nargin < 1, include_slow = false; end
slow  = {'example_07_external_stress_comparison.m'};
here  = fileparts(mfilename('fullpath'));
files = dir(fullfile(here, 'example_*.m'));
for i = 1:numel(files)
    if ~include_slow && ismember(files(i).name, slow)
        fprintf('\n##### %s skipped (run_all_examples(true) to include) #####\n', files(i).name);
        continue
    end
    fprintf('\n##### %s #####\n', files(i).name);
    run_one(fullfile(here, files(i).name));
end
end

function run_one(path)
% Separate workspace, so the examples' "clear" does not affect the loop.
run(path);
end
