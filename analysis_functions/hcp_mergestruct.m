function s = hcp_mergestruct(varargin)

% HCP_MERGESTRUCT combines multiple structures that can have different fields
% into a single structure array.
%
% Use as
%   s = hcp_mergestruct(a, b, ...)
% 
% The resulting struct "s" contains all fields of all input structures combined.

% Copyright (C) 2004, Robert Oostenveld
%
% $ Log$

s = struct;
for i=1:length(varargin)
  fn = fieldnames(varargin{i});
  for j=1:length(fn)
    s = setfield(s, fn{j}, getfield(varargin{i}, fn{j}));
  end
end

