function [staircase, stairParams] = usePalamedesStaircase(stairParams,response)

% Purpose:  Initialize and update adaptive procedures.
%           This function simply calls the Palamedes toolbox, but makes
%           it easier to set and change parameters.
%
%           You can pass in a stairParams structure.
%           e.g., [staircase stairParams] = usePalamedesStaircase(stairParams);
%
%           If run with two inputs, the first input must be 'staircase'
%           (i.e., the output of initialization) and the second input must
%           be either 1 or 0 for correct and incorrect responses, respectively.
%           This will update the adaptive method.
%           e.g., staircase = usePalamedesStaircase(staircase,1);
%
% ----------------------------------------------------------------------
% Input(s)
% stairParams           : structure containing parameters (listed below) to control adaptive method
%     whichStair        : 1=best PEST; 2=QUEST     (default=1)
%     alphaRange        : vector of possible threshold stimulus values  (default=0.01:0.01:1)
%     fitBeta           : slope of underlying psychometric function     (default=2)
%     fitLambda         : lapse rate (i.e., 1-upper asymptote) of psychometric function (default=0.01)
%     fitGamma          : guess rate (i.e., lower asymptote) of psychometric function   (default=0.5)
%     threshPerformance : target threhsold performance (must be specified if using arbWeibull, see PF)
%     lastPosterior     : posterior distribution from earlier run. when inputted, adaptive method will continue where previous run stopped (default=[])
%     PF                : shape of underlying psychometric function. can be: PAL_Weibull, PAL_Quick, PAL_Gumbel, PAL_HyperbolicSecant, PAL_Logistic, or arbWeibull(default)
%     updateAfterTrial  : if >0, the adaptive method will not use posterior-based estimates until the trial number matches the input for this variable (default=0)
%     preUpdateLevels   : these are the stimulus levels that will be tested before the adaptive method is updated (only works if updateAfterTrial>1)
% response              : 1=correct; 0=incorrect
%
% ----------------------------------------------------------------------
% Output(s)
% staircase             : structure controlling Palamedes adaptive method
% stairParams           : structure of parameters used to initialize Palamedes
% ----------------------------------------------------------------------
% Function created by Michael Jigo
% Last update : Feb. 17 2021
% ----------------------------------------------------------------------

if nargin == 0 % Use default staircase parameters
    paramNames = {'whichStair' 'alphaRange' 'fitBeta' 'fitLambda' 'fitGamma' 'threshPerformance' 'lastPosterior' 'PF' 'updateAfterTrial' 'preUpdateLevels'};
    for p = 1:numel(paramNames), stairParams.(paramNames{p}) = []; end
    stairParams = validateParams(stairParams);
    staircase = initStaircase(stairParams);
    
elseif nargin==1 % if a structre of parameters are passed in
    stairParams = validateParams(stairParams);   
    staircase = initStaircase(stairParams);
    
elseif nargin==2 % Update the staircase
%     the first input is the initialized staircase structure (NOT the stairParams structure) 
%     and the second is the response (0=incorrect or 1=correct)
    
    if ismember(func2str(stairParams.PF),{'PAL_Weibull', 'PAL_Quick', 'PAL_Gumbel', 'PAL_HyperbolicSecant', 'PAL_Logistic'})
        staircase = PAL_AMRF_updateRF(stairParams,stairParams.xCurrent,response);
    else
        staircase = update_arbitrary_pf(stairParams, stairParams.xCurrent, response, stairParams.threshPerformance);
    end
end

