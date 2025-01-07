function [t,str] = timekeeper(t,i,N,unit)
  %a helper function to deal with keeping track of time and predict ETA for a for loop
  %Before the for loop begins, run t = timekeeper; (without argument)
  %In the for loop, use the function at the end of the for loop, provide i, the current step (or more broadly the current index, physical time etc), N, the final step number (or more broadly, the final index, physical time etc)
  %t is a struct that keeps the information from one step to the next.
  %t.start is the tic returned at the initialization, starttime is the datetime returned at the initialization
  %t.elapsed is the elapsed time in seconds since the beginning
  %t.nrun is the number of timekeeper invoked after initialization)
  %t.step is the time taken since the last invoke
  %t.avg is the average time per invoke so far
  %t.istart is the index of the first step. if saved intermittently, user needs to modify manually
  %t.avgperstep is the average time per step so far
  %the three above are duration.
  %t.ETA is the datetime of the predicted ETA (only computed if N is provided)
  %use unit to set the unit to display: ('seconds', 'minutes', 'hours')
  %if not provided, unit is determined automatically by t.avg, 'seconds' if t.avg < 300 s, 'hours' if t.avg > 300 min, 'minutes' if 6min < t.avg < 300min
  %str will be displayed if not output
  if nargin < 2
    %initialize
    t.start = tic;
    t.starttime = datetime;
    t.elapsed = seconds(0);
    t.nrun = 0;
    t.istart = [];
  else
    if isempty(t.istart)
      t.istart = i;
    end
    tnew = seconds(toc(t.start));
    tstep = tnew-t.elapsed;
    %note that the user doesn't have to invoke timekeeper at each step, especially if one step is very fast. Instead, user can do this intermittently.
    %tos tavg is actually the average time per invoke.
    t.nrun = t.nrun + 1;
    tavg = tnew / t.nrun;
    t.elapsed = tnew;
    t.step = tstep;
    t.avg = tavg;
    avgperstep = tnew / (i-t.istart+1);
    t.avgperstep = avgperstep;

    if nargin < 4 || isempty(unit)
      if tavg < seconds(300)
        unit = 'seconds';
      elseif tavg > minutes(300)
        unit = 'hours';
      else
        unit = 'minutes';
      end
    end
    switch unit
    case 'seconds'
      unitalias = 's';
    case 'minutes'
      unitalias = 'min';
    case 'hours'
      unitalias = 'hr';
    otherwise
      error('unsupported unit');
    end

    str = ['Step ',num2str(i),' took ',sprintf('%4.1f',feval(unit,tstep)),' ',unitalias,...
           ', avg ',sprintf('%4.1f',feval(unit,tavg)),' ',unitalias];
    if nargin > 2
      t.ETA = avgperstep*(N-i)+tnew+t.starttime;
      str = [str,', ETA ',datestr(t.ETA,'mm/dd HH:MM AM')];
    end
    if nargout < 2
      disp(str);
    end
  end
