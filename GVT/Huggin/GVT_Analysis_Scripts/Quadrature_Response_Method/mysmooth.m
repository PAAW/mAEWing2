function Gs = mysmooth(G,Nwin)
% function Gs = mysmooth(G,Nwin)
%
% This function smooths frequency response data.  The smoothed
% value Gs at frequency i is obtained by taking the average value of G
% over a window i-Nwin:i+Nwin.
%
% Inputs:
%    G := Frequency response data from GVT_Load.
%    Nwin := Window for smoothing. Larger values of Nwin increase
%      the smoothing. [Default = 3]
%
% Outputs:
%    Gs := Smoothed frequency response. Gs has same dimensions as G.
%


% Input parsing
if nargin==2
    Nwin = 3;
end

% Smoothing
[Nomeg,Nresp] = size(G(:,:));
Gs = zeros(Nomeg,Nresp);
for i=1:Nomeg
    sidx = max(1,i-Nwin):min(i+Nwin,Nomeg);
    Gs(i,:) = mean( G(sidx,:) );
    %Gs(i,:) = median( G(sidx,:) );
end
Gs = reshape(Gs,size(G));