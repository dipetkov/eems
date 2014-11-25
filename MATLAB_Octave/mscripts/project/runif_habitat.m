

function coord = runif_habitat(n,habitat)


xmin = habitat.xmin;
xmax = habitat.xmax;
ymin = habitat.ymin;
ymax = habitat.ymax;

coord = zeros(n,2);

for k = 1:n
  in = 0;
  while (in==0)
    x = xmin+(xmax-xmin)*rand(1);
    y = ymin+(ymax-ymin)*rand(1);
    in = is_in_habitat(habitat,[x,y]);
  end
  coord(k,1) = x;
  coord(k,2) = y;
end
