function add_iFluid_path(iFluid_root_path)
% Add all iFluid subfolders to Matlab path

    dirs = {'Fluid', 'iFluidZero', 'models', 'modules', 'solvers', 'utils'};

    for i = 1:numel(dirs)
        addpath(genpath([iFluid_root_path filesep dirs{i}]));
    end

end
   