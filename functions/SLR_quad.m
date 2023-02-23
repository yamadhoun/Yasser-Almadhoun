function [SLRquad] = SLR_quad(SLRise)
%
%
% SLR_quad
%
% SLR quadratic curve from start time point to end time point
% close all;
% clear all;
% clc;

%% INPUT:
%      timestart       : the start time of the SLR curve  in [eg year]
%      timeend         : the end time of the SLR curve  in [eg year]
%      SLRstart        : the SLR at the start time of the SLR curve  in [m]
%      SLRend          : the SLR at the end time of the SLR curve  in [m]

%% OUTPUT:
%      times           : the row of times of the SLR curvce in [eg year]
%      SLR             : the row of SLRs of the SLR curvce in [m]

%% CALCULATIONS:
% y = a*(x-h)^2 + k
% P1(h, k) = P1(timestart, SLRstart)
% P2(x, y) = P1(timeend, SLRend)

h = SLRise.timestart;
k = SLRise.SLRstart;
x = SLRise.timeend;
y = SLRise.SLRend;

a = (y - k) / (x - h)^2;

x = floor(linspace(h,x,SLRise.SWLength));
y = a*(x-h).^2 + k;

times = x';
SLRvalue = y';

SLRquad = [times,SLRvalue];

end
