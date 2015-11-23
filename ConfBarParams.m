%-------------------------------------------------------------------------%
% ConfBarParams
% Confidence Bar Coordinates
%-------------------------------------------------------------------------%
function [ConfidenceBarCoords] = ConfBarParams(xCenter, yCenter,screenXpixels,screenYpixels)

% Here we set the size of our confidence bar
BarLength = 0.5*screenXpixels; % 50% of the width screen
HeightOfBarWhisks = 0.05*screenYpixels; % 5% of the height of screen

LeftExtent  = xCenter-BarLength/2;
MidLeftExtent = xCenter-BarLength/4;
RightExtent = xCenter+BarLength/2 ;
MidRightExtent = xCenter+BarLength/4;
BottomExtent = yCenter+HeightOfBarWhisks/2 ;
TopExtent   =  yCenter- HeightOfBarWhisks/2 ;

HorizontalBarCoords   = [LeftExtent RightExtent; yCenter yCenter];

LeftWhisk             = [LeftExtent LeftExtent ; BottomExtent TopExtent];
RightWhisk            = [RightExtent RightExtent ; BottomExtent TopExtent];
MidLeftWhisk          = [MidLeftExtent MidLeftExtent ; BottomExtent TopExtent];
MidRightWhisk         = [MidRightExtent MidRightExtent ; BottomExtent TopExtent];
CenterWhisk           = [xCenter xCenter ; BottomExtent TopExtent];

ConfidenceBarCoords   = [HorizontalBarCoords LeftWhisk RightWhisk ...
    MidLeftWhisk MidRightWhisk CenterWhisk];
end