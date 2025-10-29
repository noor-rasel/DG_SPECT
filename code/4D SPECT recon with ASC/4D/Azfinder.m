% erwan Gravier
% research project
% 06/09/04
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this function is the implementation of the 
% Channelized Observer. Will calculate the SNRcho
% and then deduce Az from it.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% used for MOTION COMPENSATED EM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% H1: no lesion
% H2: lesion present
% Rule: choose H2 if Likelyhood (Liky) >= lambda
%              H1 otherwise
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Az=Azfinder(stocknoles,stockles,noiselevel)
% Constant and variables
Sstocknoles=size(stocknoles);
numchannels=Sstocknoles(2);    % number of channels

% calculation of internal noise
internnoise=eye(numchannels).*noiselevel;

% calculation of the templates.
templnoles=mean(stocknoles)';
templles=mean(stockles)';
deltas=templles-templnoles;


% calculation of the covariance matrix
covnoles=cov(stocknoles);
covles=cov(stockles);
covuse=0.5*(covnoles+covles)+internnoise;
covuse=pinv(covuse);

% now calculating SNRcho
SNRcho=sqrt((deltas'*covuse*deltas));

% calculating Az, the area under the ROC curve
Az=0.5*(1+erf(SNRcho/2));



