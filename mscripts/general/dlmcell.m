

function dlmcell(fileName,cellArray)


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
    printrowvec(datei,v);
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
    printrowvec(datei,v);
  end
  % Close file
  fclose(datei);
else
  disp('Error: Input one-dimensional cell array of vectors.');
  return;
end


function printrowvec(filept,vector)
n = length(vector);
for i=1:n
  fprintf(filept,'%f',vector(i));
  if (i<n)
    fprintf(filept,' ');
  else
    fprintf(filept,'\n');
  end
end

