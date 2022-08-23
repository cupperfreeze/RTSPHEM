% Compute liquid phase distribution in 3D image bw with saturation desSat
% bw is expected to encode pore-space as 0 and matrix structure as 1
function bw = MorphModel(bw, desSat)


caseSwitch = 1;			%Case1 = maximal Spheres, Case0 = wetting layer of uniform thickness
bw=1-bw;

if caseSwitch
Imsize = size(bw);			
[X,Y,Z]=ndgrid(1:Imsize(1),1:Imsize(2),1:Imsize(3));
initialPore = sum(bw(:)==1);
currentSat = 1;
breakVar = false;
distMap = bwdist(1-bw,'euclidean'); 				%Calculate Euclidean distance amp
[ordered,indizes] = sort(distMap(:),'descend');
InverseIndizes(indizes) = 1:length(indizes);
idx=[];
while currentSat>desSat  
   
    %distMap(idx) = 0;
    %[radius, center] = max(distMap(:));
    ordered(InverseIndizes(idx)) = 0;
    [row,col,radius]= find(ordered,1);
    center = indizes(row);
    
    radius =radius-1;
    if radius<1
       disp('Warning: Desired saturation Level could not be reached');
       return
    end
    [centerIdx1,centerIdx2,centerIdx3] = ind2sub(Imsize,center);
    center = [centerIdx1,centerIdx2,centerIdx3];
    while true
     %  idx = ((center(1)-X(:)).^2 + (center(2)-Y(:)).^2 + (center(3)-Z(:)).^2 - radius^2)<0;
        % Start working locally on Cube around current Sphere
        xIdx = max((center(1)-ceil(radius)),1):min((center(1)+ceil(radius)),Imsize(1)); 
        yIdx = max((center(2)-ceil(radius)),1):min((center(2)+ceil(radius)),Imsize(2)); 
        zIdx = max((center(3)-ceil(radius)),1):min((center(3)+ceil(radius)),Imsize(3)); 
        [X,Y,Z]=ndgrid(xIdx,yIdx,zIdx);
        idxloc = ((center(1)-X(:)).^2 + (center(2)-Y(:)).^2 + (center(3)-Z(:)).^2 - radius^2)<0;
        idx = X(idxloc)+(Y(idxloc)-1)*Imsize(2) + (Z(idxloc)-1)*Imsize(2)*Imsize(3);         
       
        
        bwOld = bw;
        currentSatOld = currentSat;
        overlap = sum(bw(idx)==0);
        bw(idx)=0;
        currentSat = currentSat - (length(idx)-overlap)/initialPore;
        if currentSat>desSat
            if breakVar 
                return 
            end
            break
        else										%Case when current sphere would undershoot desSat
            bw=bwOld;
            currentSat = currentSatOld;
            radius=radius-0.5;
            breakVar = true;							
        end
    end
end

else
	Imsize = size(bw);
	[X,Y,Z]=meshgrid(1:Imsize(1),1:Imsize(2),1:Imsize(3));
	initialPore = sum(bw(:)==1);
	currentSat = 1;
	distMap = bwdist(1-bw,'euclidean');    
	[sorted,indizes] = sort(distMap(:));
	breakIndex = round(prod(Imsize)-(1-desSat)*initialPore);
	bw(indizes(breakIndex:end))=0;
end
end
