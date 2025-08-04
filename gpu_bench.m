% GPU vs CPU Matrix Multiplication Benchmark

% --- Configuration ---
matrixSizes = [100, 200, 500, 1000, 1500, 2000, 2500, 3000, 4000, 5000, 10000]; % Square matrix dimensions to test
numRepetitions = 5; % Number of times to repeat each multiplication for averaging

% --- Initialization ---
cpuTimes = zeros(length(matrixSizes), 1);
gpuTimesWithTransfer = zeros(length(matrixSizes), 1);
gpuTimesComputeOnly = zeros(length(matrixSizes), 1);

fprintf('Starting Matrix Multiplication Benchmark...\n');
fprintf('-------------------------------------------\n');
fprintf('%12s | %15s | %20s | %20s\n', 'Matrix Size', 'CPU Time (s)', 'GPU Time (Total, s)', 'GPU Time (Compute, s)');
fprintf('-------------------------------------------\n');

% --- Benchmarking Loop ---
for i = 1:length(matrixSizes)
    N = matrixSizes(i);
    fprintf('%10dx%-1d | ', N, N);

    % Generate random matrices on the CPU
    A_cpu = rand(N, N);
    B_cpu = rand(N, N);

    % --- CPU Benchmark ---
    tempCpuTime = 0;
    for k = 1:numRepetitions
        tic;
        C_cpu = A_cpu * B_cpu;
        tempCpuTime = tempCpuTime + toc;
    end
    cpuTimes(i) = tempCpuTime / numRepetitions;
    fprintf('%15.4f | ', cpuTimes(i));

    % --- GPU Benchmark (with data transfer) ---
    tempGpuTimeWithTransfer = 0;
    try
        % Warm-up GPU for the first iteration of a new size (optional but good practice)
        if i == 1 || (i > 1 && matrixSizes(i-1) < 1000 && N >=1000) % Heuristic for warm-up
            A_gpu_warmup = gpuArray(rand(100,100));
            B_gpu_warmup = gpuArray(rand(100,100));
            C_gpu_warmup = A_gpu_warmup * B_gpu_warmup;
            wait(gpuDevice); % Wait for GPU to finish
            clear A_gpu_warmup B_gpu_warmup C_gpu_warmup;
        end

        for k = 1:numRepetitions
            tic;
            A_gpu = gpuArray(A_cpu); % Move A to GPU
            B_gpu = gpuArray(B_cpu); % Move B to GPU
            C_gpu = A_gpu * B_gpu;   % Perform multiplication on GPU
            C_result_from_gpu = gather(C_gpu); % Move result back to CPU
            wait(gpuDevice); % Ensure all GPU operations are complete before stopping timer
            tempGpuTimeWithTransfer = tempGpuTimeWithTransfer + toc;
        end
        gpuTimesWithTransfer(i) = tempGpuTimeWithTransfer / numRepetitions;
        fprintf('%20.4f | ', gpuTimesWithTransfer(i));

    catch ME
        fprintf('GPU Error: %s. Skipping GPU for this size.\n', ME.message);
        gpuTimesWithTransfer(i) = NaN; % Mark as Not a Number if GPU fails
        gpuTimesComputeOnly(i) = NaN;
        fprintf('%20s | %20s\n', 'N/A', 'N/A');
        continue; % Skip to next matrix size if GPU fails (e.g., out of memory)
    end


    % --- GPU Benchmark (compute only, minimizing transfer overhead in timing) ---
    tempGpuTimeComputeOnly = 0;
    try
        A_gpu = gpuArray(A_cpu); % Data transfer outside the timed loop for this part
        B_gpu = gpuArray(B_cpu);
        wait(gpuDevice); % Ensure data is on GPU

        for k = 1:numRepetitions
            tic;
            C_gpu = A_gpu * B_gpu;   % Perform multiplication on GPU
            wait(gpuDevice); % Ensure GPU computation is complete
            tempGpuTimeComputeOnly = tempGpuTimeComputeOnly + toc;
        end
        gpuTimesComputeOnly(i) = tempGpuTimeComputeOnly / numRepetitions;
        % Don't need to gather C_gpu here for timing computation, but you would for using the result
        fprintf('%20.4f\n', gpuTimesComputeOnly(i));
        clear A_gpu B_gpu C_gpu C_result_from_gpu; % Clear GPU memory

    catch ME
        fprintf('GPU Compute-Only Error: %s. Skipping for this size.\n', ME.message);
        gpuTimesComputeOnly(i) = NaN;
        fprintf('%20s\n', 'N/A');
        % No continue here as the withTransfer part might have worked
    end
end

fprintf('-------------------------------------------\n');
fprintf('Benchmark Complete.\n');

% --- Plotting Results ---
figure;
semilogy(matrixSizes, cpuTimes, 'b-o', 'LineWidth', 1.5, 'MarkerSize', 6);
hold on;
semilogy(matrixSizes, gpuTimesWithTransfer, 'r-s', 'LineWidth', 1.5, 'MarkerSize', 6);
semilogy(matrixSizes, gpuTimesComputeOnly, 'g--x', 'LineWidth', 1.5, 'MarkerSize', 6);
hold off;

title('CPU vs GPU Matrix Multiplication Benchmark');
xlabel('Matrix Size (N x N)');
ylabel('Average Execution Time (seconds, log scale)');
legend('CPU', 'GPU (with Data Transfer)', 'GPU (Compute Only)', 'Location', 'NorthWest');
grid on;
set(gca, 'XTick', matrixSizes); % Ensure all matrix sizes are shown as ticks
xtickangle(45); % Angle ticks if they overlap

% --- Display GPU Information ---
try
    gpuInfo = gpuDevice;
    fprintf('\n--- GPU Information ---\n');
    fprintf('GPU Name: %s\n', gpuInfo.Name);
    fprintf('GPU Compute Capability: %s\n', gpuInfo.ComputeCapability);
    fprintf('GPU Total Memory: %.2f GB\n', gpuInfo.TotalMemory / (1024^3));
    fprintf('-----------------------\n');
catch
    fprintf('\nNo GPU detected or Parallel Computing Toolbox not available.\n');
end