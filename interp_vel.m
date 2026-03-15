function [uwork, vwork] = interp_vel(XLag, YLag, stringweek, stringyear) 

% This function permits us to calculate the velocity of the sea ice with the Ice Interpolent data. 
% We need to switch XLag and YLag to not lose Active tracers and sea ice. 
% 
% It takes from LITS (python) the positions and the time to calculate 
% the velocity fields. 
% 

    SIVPath = [pwd '/SIV_V4_scatteredInterpolant_wSIC_V4_wSIV_V4/'];
    addpath(SIVPath);

    load([SIVPath 'IceMotionInterp_scatteredInterpolant_Week_' stringweek '.' stringyear '.mat']) ;
    % load([SICPath 'IceMotionInterpolant_Week_' stringweek '.' stringyear '.mat'])

    % For SIV_V4_scatteredInterpolant :
    uwork = interpu(XLag, YLag);  % Interpolate sea-ice veolocities
    vwork = interpv(XLag, YLag); % from the EASE grid onto our current exact

end