function [lgraph, dlnet, numLearnables] = DLN_BasalGanglia_v3(numChan, N, showplot)
% 
% v1 was not working with idNeuralStateSpace. This version removes the
% indirect pathway to make the graph a single branch. 
% 
% Create Deep Learning Network Architecture
% Network was designed based on the structure of the basal ganglia. 
% Script for creating the layers for a deep learning network with the following 
% properties:
%
% 
%  Number of layers: 13
%  Number of connections: 12
%
% 
% Run the script to create the layers in the output variable |lgraph|.
% 
% To learn more, see <matlab:helpview('deeplearning','generate_matlab_code') 
% Generate MATLAB Code From Deep Network Designer>.
% 
% Auto-generated by MATLAB on 15-Jan-2025 02:00:45
% 
% Inputs: 
%   numChan: number of ECoG channels. This determines the dimension of
%            input and output. Default = 63, i.e. 21x3 grid with no
%            missing channels.
%   N: scaling factor that can proportionally increase the dimension of all
%      layers except input and output. Default = 1
%   showplot: if true [default], lgraph will be plotted in a new figure
% 
% Outputs: 
%   lgraph: LayerGraph object 
%   dlnet: dlnetwork object without output (regression) layer 
%   numLearnables: # of learnables 
% 
%% handle incomplete inputs 
if nargin < 3
    showplot = true;
    if nargin < 2
        N = [];
        if nargin < 1
            numChan = [];
        end
    end
end
if isempty(numChan)
    numChan = 63;
end
if isempty(N)
    N = 1;
end

%% Create Array of Layers

initval = @(sgn, h, w) sgn/(w+h) + sqrt(2/(w+h))*randn(h,w);

layers = [
    featureInputLayer(numChan,"Name","ECoG")
    fullyConnectedLayer(N*70,"Name","fcCortexInput","Weights",initval(0,N*70,numChan))
    tanhLayer("Name","actCortexInput")
    fullyConnectedLayer(N*70,"Name","fcStriatum","Weights",initval(1,N*70,N*70))
    tanhLayer("Name","actStriatum")
    fullyConnectedLayer(N*30,"Name","fcGPinput","Weights",initval(-1,N*30,N*70))
    tanhLayer("Name","actGPinput")
    fullyConnectedLayer(N*10,"Name","fcGPout","Weights",initval(0,N*10,N*30))
    tanhLayer("Name","actGPout")
    fullyConnectedLayer(N*250,"Name","fcThalamusVA","Weights",initval(-1,N*250,N*10))
    tanhLayer("Name","actThalamusVA")
    fullyConnectedLayer(numChan,"Name","fcCortexOutput","Weights",initval(1,numChan,N*250))
    regressionLayer("Name","regressionoutput")];

%% Create Layer Graph
% Create the layer graph variable to contain the network layers.

lgraph = layerGraph(layers);

%% dlnet and # of learnables 

% for some reason, matlab doesn't want to output the number of learnables
% of the lgraph, but it will do so for the dlnet 
% and for some reason, conversion does not work unless the regressionoutput
% layer is removed 
dlnet = dlnetwork(removeLayers(lgraph, 'regressionoutput'));
numLearnables = [dlnet.Learnables.Value]; % cell array for each layer 
numLearnables = arrayfun(@(l) numel(l{:}), numLearnables); % total # at each layer
numLearnables = sum(numLearnables); % grand total # 

%% Plot Layers

if showplot
    figure;
    plot(lgraph);
end

end