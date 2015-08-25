function activeMouseID = getMouseID
%
% script to mouse response pad

d = PsychHID('Devices');

x = strfind({d.usageName},'Mouse');
activeMouseID = find(~cellfun(@isempty,x));
