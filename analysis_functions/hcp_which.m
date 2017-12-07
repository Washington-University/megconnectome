function filename = hcp_which(original)

% HCP_WHICH is a replacement function for the builtin MATLAB which
% function that works in deployed mode.

if ~isempty(dir(original))
  % the file exists as it is, no need for any changes
  filename = original;
  return
end

% try to find the file on the path
p = tokenize(path, ':');

for i=1:length(p)
  filename = fullfile(p{i}, original);
  if ~isempty(dir(filename))
    break
  end
end

% still not found, try it with a *.m extension
for i=1:length(p)
  filename = fullfile(p{i}, [original '.m']);
  if ~isempty(dir(filename))
    break
  end
end

% still not found, try it with a *.mat extension
for i=1:length(p)
  filename = fullfile(p{i}, [original '.mat']);
  if ~isempty(dir(filename))
    break
  end
end

% still not found, return the original name
filename = original;
warning('could not find the file "%s" on the path', filename);

