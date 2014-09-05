

function ans = is_in_habitat(habitat,coord)
%% An irregularly shaped habitat is represented as a black-and-white image %%
%% So to check whether a tile center is inside the habitat,                %%
%% first assign the point to a pixel in the image and then check the color %%
%% of the pixel (black - inside, white - outside)                          %%


xmin = habitat.xmin;
xmax = habitat.xmax;
ymin = habitat.ymin;
ymax = habitat.ymax;
xcoord = coord(:,1);
ycoord = coord(:,2);
n = nrow(coord);
ans = zeros(n,1);

if (habitat.isregular)
  ans = (xcoord>=xmin) & (xcoord<=xmax) & ...
	(ycoord>=ymin) & (ycoord<=ymax);
else
  %% First assign the point to a pixel in the image %%
  vi = round(habitat.nx*(xcoord-xmin)/(xmax-xmin));
  vj = round(habitat.ny*(ycoord-ymin)/(ymax-ymin));
  for k = 1:n
    i = vi(k);
    j = vj(k);
    if ( (i>=1) && (i<=habitat.nx) && ...
         (j>=1) && (j<=habitat.ny) )
      ans(k) = habitat.image(j,i);
    end
  end
end
