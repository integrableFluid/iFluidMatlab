function use_GPU = GPU_mode_on()

    global enable_iFluid_GPU

    if exist('enable_iFluid_GPU', 'var')
        if enable_iFluid_GPU
            use_GPU = true;
        else 
            use_GPU = false;
        end
    else
        use_GPU = false;
    end

end