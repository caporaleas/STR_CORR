function [perdiff] = calc_perdiff(inputArg_50,inputArg_100)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here
perdiff = ((inputArg_50-inputArg_100)./inputArg_100).*100;
end