function env_bin = init_python_env_matlab()
    %-----------------------------------------------------------------
    % Set the Python environment (Python 3.11)
    % Windows format
    env_bin = 'C:\Users\selim\miniconda3\envs\scanpy_env_311\python.exe';
    %env_bin = 'C:\Local_install\miniconda3\envs\scanpy_env_311\python.exe';
    if ispc
        env_bin = strrep(env_bin,"\","\\");
    end
    % Linux format
    %env_bin = "/home/ssromerogon/packages/scanpy_env/bin/python3";
    %-----------------------------------------------------------------
    % Load python environment
    % Clear any existing Python environment to force reinitialization
    pe = pyenv('Version', env_bin);
    
    % Check if the environment is loaded
    if pe.Status ~= "Loaded"
        fprintf("Reinitializing Python environment...\n");
        pe = pyenv('Version', env_bin);
        %pause(20);  % Optional: Wait for 1 second?
        % Load the environment by executing a simple Python command
        py.exec('import sys');
    end
    % Display the environment details
    disp(pyenv);

end