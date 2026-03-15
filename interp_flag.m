function [ActiveFlag_out] = interp_flag(sicEASE, XLag, YLag) 

% This function is used to determine the ActiveFlag after each run.
% It takes from LITS (in python) ActiveFlace, sicEASE (concentration), XLag, YLag
% and n to determine the ActiveFlag at time m. 
% 
% The vectors ActiveFlag, sicEASE, XLag, YLag are already formated to the right n.
% 
% Here, XLag and YLag are switched to have the right orientation, if we put them the "right" way, 
% we are losing too much ice and Active tracers. 


    eps = 0.15 ;

    ActiveFlag_out = (interp2(1:361, 1:361, sicEASE, YLag, XLag) > eps);
