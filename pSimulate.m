function [xhat, yhat, that] = pSimulate(this,tspan,x0,NameValueArgs)
% 
% edit Toren Arginteanu, January 2025 
% Previous version resulted in errors when simulating a neural state space
% autonomous system in MATLAB 2023b.
% 
% idNeuralStateSpace hidden function.

%   Author: R. Chen, R. Singh
%   Copyright 2022 The MathWorks, Inc.

% pSimulate simulates state and output trajectories using an ODE45 solver based
% on initial state and input sequence at the specific time points.
%
%   Use one of the following syntax depending on system configuration:
%
%       [x, y, t] = sim(nss, tspan, x0)
%       [x, y, t] = sim(nss, tspan, x0, Input=u)
%       [x, y, t] = sim(nss, tspan, x0, Input=u, InputInterSample=method)
%
%   where inputs are:
%
%           nss: an idNeuralStateSpace object.
%         tspan: a two-element vector [T0 TFINAL] or a vector with
%                several time points [T0 T1 ... TFINAL]. For a
%                continuous-time system, if you specify more than two time
%                points, ode45 returns interpolated solutions at the
%                requested times.  For a discrete time system, time points
%                must be multiples of the system sample time "sys.Ts".
%            x0: a vector for initial state
%             u: sequence of plant input signal described by a timetable, a
%                timeseries or a structure with Time and Value fields.
%                Only required when system has inputs.
%        method: "zoh","foh","cubic","makima","pchip", or "spline",
%                indicating how to interpolate input signal during
%                integration (default = "zoh").
%
%   Outputs "x" and "y" are state and output trajectories and "t" is the
%   corresponding time vector.  When there is no output argument, response
%   plots are generated.

arguments
    this
    tspan double {mustBeVector, mustBeFinite, mustBeGreaterThanOrEqual(tspan,0)}
    x0 double {mustBeVector}
    NameValueArgs.Input = []
    NameValueArgs.InputInterSample (1,1) string {mustBeMember(NameValueArgs.InputInterSample,["previous","linear","cubic","makima","pchip","spline","zoh","foh","bl"])} = "previous"
end
%% state network
% validate inputs
if this.HasInputs
    [Vu,Tu] = checkDataSetForValidation(NameValueArgs.Input);
else
    Tu = [];
end
% get options
options = getOptions(this);
switch NameValueArgs.InputInterSample
    case "zoh"
        options.InputInterSample = 'previous';
    case "foh"
        options.InputInterSample = 'linear';
    case "bl"
       error('BL not supported (requires ideal interp)')
    otherwise
        options.InputInterSample = NameValueArgs.InputInterSample;
end
if numel(Tu)>1
   options.Interpolant = griddedInterpolant(Tu,Vu',options.InputInterSample);
end
if this.HasInputs
    % predict state trajectory
    [xhat, that] = simState(this.StateMLPParameters,options,x0,tspan);
    % predict output trajectory
    if this.HasIdentityOutputFcn
        yhat = xhat;
    else
        yhat = [xhat; predict_Y_XU(this.OutputMLPParameters,options,xhat,that)];
    end
else
    % predict state trajectory
    [xhat, that] = simState(this.StateMLPParameters,options,x0,tspan);
    % predict output trajectory
    if this.HasIdentityOutputFcn
        yhat = xhat;
    else
        yhat = [xhat; predict_Y_X(this.OutputMLPParameters,options,xhat,that)];
    end
end
if nargout==0
    % plot state trajectory
    figure
    for ct=1:this.NumStates
        subplot(this.NumStates,1,ct)
        plot(that,xhat(ct,:),'b-');
        ylabel("x("+ct+")");
        if ct==1
            title('state trajectory')
        end
        if ct==this.NumStates
            xlabel("time ("+this.TimeUnit+")");
        end
    end
    % plot output trajectory
    figure
    for ct=1:this.NumOutputs
        subplot(this.NumOutputs,1,ct)
        plot(that,yhat(ct,:),'b-');
        ylabel("y("+ct+")");
        if ct==1
            title('output trajectory')
        end
        if ct==this.NumOutputs
            xlabel("time ("+this.TimeUnit+")");
        end
    end
else
    xhat = xhat';
    yhat = yhat';
    that = reshape(that,length(that),1);
end
