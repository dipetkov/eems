

function dlmcell(fileName,cellArray,varargin)


format = '%f';
delimiter = ' ';

i = 0;
while (i<length(varargin))
  i = i + 1;
  if (strcmp(varargin{i},'precision'))
    precision = varargin{i+1};
    if (isnumeric(precision))
      format = strcat('%.',num2str(precision),'f');
    else
      format = precision;
    end
  elseif (strcmpi(varargin{i},'delimiter'))
    delimiter = varargin{i+1};
  end
  i = i + 1;
end


nr = size(cellArray,1);
nc = size(cellArray,2);


if (nc==1)
  %% Open file to write
  datei = fopen(fileName, 'w');
  for r=1:nr
    v = cellArray{r,1};
    if ~isvector(v)
      disp('Error: Input one-dimensional cell array of vectors.');
      return;
    end
    printrowvec(datei,v,format,delimiter);
  end
  % Close file
  fclose(datei);
elseif (nr==1)
  %% Open file to write
  datei = fopen(fileName, 'w');
  for c=1:nc
    v = cellArray{1,c};
    if ~isvector(v)
      disp('Error: Input one-dimensional cell array of vectors.');
      return;
    end
    printrowvec(datei,v,format,delimiter);
  end
  % Close file
  fclose(datei);
else
  disp('Error: Input one-dimensional cell array of vectors.');
  return;
end


function printrowvec(filept,vector,format,delimiter)
n = length(vector);
for i=1:n
  fprintf(filept,format,vector(i));
  if (i<n)
    fprintf(filept,delimiter);
  else
    fprintf(filept,'\n');
  end
end
