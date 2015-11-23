%-------------------------------------------------------------------------%
% sendMarker
% sends LabStreamLayer markers to NIC on stimulation computer
%-------------------------------------------------------------------------%
function sendMarker(marker,LSLOutlet)
    ret=MatNICMarkerSendLSL(marker,LSLOutlet);
    if ret<0
        warning('error sending marker')
    end
end