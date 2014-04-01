

function [habitat,inPops] = read_rectangular_habitat(datapath,Vcoord,xPop,yPop)


dimns = dlmread(strcat(datapath,'.dimns'));

xmin = dimns(1,1);
xmax = dimns(1,2);
ymin = dimns(2,1);
ymax = dimns(2,2);

habitat = struct('xmin',{xmin},'xmax',{xmax},...
                 'ymin',{ymin},'ymax',{ymax},...
		 'xPop',{xPop},'yPop',{yPop});
inPops = ones(nrow(Vcoord),1);
