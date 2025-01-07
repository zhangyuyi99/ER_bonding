function visualize_1d()
    % VISUALIZE  Shows a time-lapse of the 1D particle distribution as a histogram
    %   This script:
    %   1) Loads 'test_result.mat' which contains xout and param
    %   2) Creates a figure for each time step (or a subset),
    %   3) Displays the histogram of particle positions,
    %   4) Optionally saves the frames into a video or GIF.
    %
    % Adjust bin size, frequency of frames, etc. as needed.

    % 1) Load results
    data = load('result\test_result.mat', 'xout', 'param');
    if ~isfield(data, 'xout') || ~isfield(data, 'param')
        error('Missing xout or param in test_result.mat');
    end
    xout = data.xout;       % size (N, ndims, steps)
    param = data.param;

    % Make sure we are in 1D
    if param.ndims ~= 1
        error('visualize: This example is for 1D only.');
    end

    % Setup for binning
    L = param.L_ER;
    nBins = 100;          % number of histogram bins (feel free to change)
    edges = linspace(0, L, nBins+1);   % bin edges from 0..L
    midpoints = 0.5 * (edges(1:end-1) + edges(2:end));

    % We'll produce a figure showing the distribution
    fig = figure('Name','ER Simulation','Color','white');

    % Optional: Create a VideoWriter or GIF
    doVideo = true;   % or false if you just want to view
    if doVideo
        % Option 1: save to an MP4
        v = VideoWriter('result\ER_sim.mp4','MPEG-4');
        open(v);
        % (Alternatively, you can do a GIF approach using imwrite in a loop.)
    end

    % We'll skip frames if too many steps
    totalSteps = size(xout,3);
    skipFrame = 1;      % e.g. visualize every 1 step
    % If many steps, consider skipFrame = 10 or so

    for st = 1:skipFrame:totalSteps
        % Extract positions at this step (N x 1)
        x = xout(:,:,st);

        % Optionally wrap positions in [0,L) for consistency
        % x = mod(x, L);

        % Build histogram (count how many fall in each bin)
        counts = histcounts(x, edges);

        % Convert counts -> density if you like
        % e.g. density = counts / (sum(counts) * binWidth) or something
        % Here we just plot counts
        plot(midpoints, counts, 'b-', 'LineWidth', 2);
        title(sprintf('Step %d / %d', st, totalSteps));
        xlabel('Position along ER');
        ylabel('Particle count');
        ylim([0, 20]);  % a little headroom

        drawnow;
        pause(0.01);   % short pause for interactive viewing

        if doVideo
            frame = getframe(fig);
            writeVideo(v, frame);
        end
    end

    if doVideo
        close(v);
        fprintf('Video saved to ER_sim.mp4\n');
    end
end
