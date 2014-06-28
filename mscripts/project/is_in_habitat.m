

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

ans = (xcoord>=xmin) & (xcoord<=xmax) & ...
      (ycoord>=ymin) & (ycoord<=ymax);
