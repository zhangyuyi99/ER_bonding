function visualize_1d_bind()
    % VISUALIZE_1D_BIND  Shows a time-lapse of the 1D particle distribution, 
    % separating bound vs. unbound particles into two lines.
    %
    %   1) Loads 'test_result.mat' which should contain xout, param, and isBoundout.
    %   2) Creates a figure for each time step (or a subset),
    %   3) Displays two histograms on the same plot: one for bound, one for unbound,
    %   4) Optionally saves the frames into a video.
    %
    % Adjust bin size, frequency of frames, etc. as needed.

    % 1) Load results
    data = load('result\test_result.mat', 'xout', 'param', 'isBoundout');
    if ~isfield(data, 'xout') || ~isfield(data, 'param') || ~isfield(data,'isBoundout')
        error('Missing xout or param or isBoundout in test_result.mat');
    end
    xout       = data.xout;       % size (N, ndims, steps)
    param      = data.param;
    isBoundout = data.isBoundout; % size (N, steps)

    % Make sure we are in 1D
    if param.ndims ~= 1
        error('visualize: This example is for 1D only.');
    end

    % Setup for binning
    L     = param.L_ER;
    nBins = 50;                     % number of histogram bins 
    edges = linspace(0, L, nBins+1);% bin edges from 0..L
    midpoints = 0.5 * (edges(1:end-1) + edges(2:end));

    % We'll produce a figure showing two lines: bound / unbound
    fig = figure('Name','ER Simulation','Color','white');

    % Optional: Create a VideoWriter or skip
    doVideo = true;  % or false 
    if doVideo
        v = VideoWriter('result\ER_sim_bind.mp4','MPEG-4');
        open(v);
    end

    totalSteps = size(xout, 3);
    skipFrame  = 1;  % e.g. visualize every step

    for st = 1:skipFrame:totalSteps
        % Extract positions at this step (N x 1)
        x   = xout(:,:,st);
        % Extract which particles are bound or unbound
        % isBoundout is (N, steps), so get column st
        isBoundNow = isBoundout(:,:,st);
        isBoundNow = logical(isBoundNow);

        % Separate positions
        xBound    = x(isBoundNow);
        xUnbound  = x(~isBoundNow);

        % Build histograms
        countsBound   = histcounts(xBound,   edges);
        countsUnbound = histcounts(xUnbound, edges);

        % Plot two lines on the same axes
        plot(midpoints, countsBound,   'r-', 'LineWidth', 2); 
        hold on
        plot(midpoints, countsUnbound, 'b-', 'LineWidth', 2);
        hold off

        xlabel('Position along ER');
        ylabel('Particle count');
        legend({'Bound','Unbound'}, 'Location','best');
        title(sprintf('Step %d / %d', st, totalSteps));

        % Adjust y-limits if needed 
        % e.g. if you know the max count might be ~N
        maxCount = max([countsBound, countsUnbound]);
        ylim([0, 400]);

        drawnow;
        pause(0.01);  % short pause for interactive viewing

        if doVideo
            frame = getframe(fig);
            writeVideo(v, frame);
        end
    end

    if doVideo
        close(v);
        fprintf('Video saved to ER_sim_bind.mp4\n');
    end
end
