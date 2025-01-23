function [xhat, that] = simState(Learnables,Options,x0,tspan)
% 
% edit Toren Arginteanu, January 2025 
% Previous version resulted in errors when simulating a neural state space
% autonomous system in MATLAB 2023b.
% 
% idNeuralStateSpace private function.

%   Author: R. Chen
%   Copyright 2022-2023 The MathWorks, Inc.

%%
% This function uses "ode45" solver (continuous-time) on the trained state
% function to generate state trajectory.  Only one experiment is allowed.
% This function is used by analysis commands such as "sim" and "predict".
%%
% x0: initial state (nx-by-1)
% tspan (length>=2):
%     when length(tspan) =2, they are start and end time points 
%        In CT, ode solver will return xhat and that at points in between
%        In DT, it returns trajectory at "start:Ts:end" time points
%     when length(tspan) >2, they are the intended time points
% Tu: input trajectory time points (1-by-#in_timepoints) when nu>0
% Vu: input trajectory values (nu-by-#input_timepoints) when nu>0
% xhat: predicted state trajectory (nx-by-#state_timepoints)
% that: column vector of time points corresponding to xhat

if length(tspan)==1
    if tspan==0
        that = tspan; xhat = x0;
    else
        keyboard
    end
else

if Options.IsContinuousTime
   % Use ode45 to compute the solution
   if Options.HasInputs
      fcn = @(t,x) forwardCTXU(t,x,Learnables,Options);
   else
      fcn = @(t,x) forwardCTX(t,x,Learnables,Options);
   end
   [that,xhat] = ode45(fcn, tspan, x0);
   % when tspan has two elements, xhat and that are determined by the ode
   % solver.  that is the reason we need to return both xhat and that
   % because tspan can be the first and last time points.
   xhat = xhat';
else
   % here tspan is desired time vector but we need to perform discrete
   % time integration to tspan(end) and then pick up the corresponding
   % values to return
   xhat(:,1) = x0;
   if length(tspan)==2
      Tx = (tspan(1):Options.Ts:tspan(end))';
   else
      Tx = tspan;
   end
   PP = Options.Interpolant;
   for ct=1:(length(Tx)-1)
      if Options.HasInputs
         UatTx = PP(Tx(:))';
         if Options.IsTimeInvariant
            xhat(:,ct+1) = forwardDT([xhat(:,ct);UatTx(:,ct)],Learnables,Options);
         else
            xhat(:,ct+1) = forwardDT([(ct-1)*Options.Ts;xhat(:,ct);UatTx(:,ct)],Learnables,Options);
         end
      else
         if Options.IsTimeInvariant
            xhat(:,ct+1) = forwardDT(xhat(:,ct),Learnables,Options);
         else
            xhat(:,ct+1) = forwardDT([(ct-1)*Options.Ts;xhat(:,ct)],Learnables,Options);
         end
      end
   end
   that = Tx;
end

end
