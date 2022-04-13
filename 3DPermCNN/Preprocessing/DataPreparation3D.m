function [counter,locations] = DataPreparation3D(image, counter, TotalSize, grid, MiniCubeSize, SlideShare, relPath)

%Extract subimages of dimension MiniCubeSize from larger scan and account for disconnected porespace.
%Save Subsample and 2 rotations in 1bit enconded raw file

% auxiliary variables
numV = uint32(grid.nodes);
coords = grid.coordinates;
neighborIndices = uint32([grid.pull{1, 1, 1}, grid.pull{2, 1, 1}, grid.pull{1, 1, 2}, grid.pull{2, 1, 2},grid.pull{1, 1, 3}, grid.pull{2, 1, 3}]);
mexEnabled = true;

% Test wheater mex is available
try 
     subimage=image((1:MiniCubeSize), (1:MiniCubeSize),(1:MiniCubeSize));   
     [im,~,isConnected] = fillholes3D_mex(numV,coords,neighborIndices, subimage(:),uint8(1));
catch
     mexEnabled = false;
end

shift = round(SlideShare * MiniCubeSize);
numberCutsPerDim = floor(1+(TotalSize-MiniCubeSize)./shift);
locations = nan(prod(numberCutsPerDim),3);                                  % save coordinates of subsample 'counter' within original cube


% iterate over subcubes
 for i=1:numberCutsPerDim(1)
     for j=1:numberCutsPerDim(2)
         for k=1:numberCutsPerDim(3)
            
                %extract subimage and the 2 rotated variants
                %close disconnected porespace
                %if no connected porespace wrt main direction --> disregard
                %if ~(mod(i,2) ==1 & mod(i,2) ==1 & mod(i,2) ==1)

                subimage=image(((i-1)*shift(1))+(1:MiniCubeSize(1)), ((j-1)*shift(2))+(1:MiniCubeSize(2)),((k-1)*shift(3))+(1:MiniCubeSize(3)));   
                if mexEnabled
                    [im,~,isConnected] = fillholes3D_mex(numV,coords,neighborIndices, subimage(:),uint8(1));
                else
                    [im,~,isConnected] = fillholes3D(numV,coords,neighborIndices, subimage(:),uint8(1));
                end            
                if isConnected
                    locations((i-1)*numberCutsPerDim(2)*numberCutsPerDim(3)+(j-1)*numberCutsPerDim(3)+k,1) = counter;
                    nameConnected = strcat('out_connected' ,num2str(counter));
                    nameOriginal = strcat('out_original' ,num2str(counter));
                    mat2rawrwd(reshape(1-im,MiniCubeSize),nameConnected, 1,[0,0,0], TotalSize./TotalSize(1), relPath);
                    mat2rawrwd(reshape(1-subimage(:),MiniCubeSize),nameOriginal, 1,[0,0,0], TotalSize./TotalSize(1), relPath);
                    counter = counter + 1;
                end

                if mexEnabled
                    [im,~,isConnected] = fillholes3D_mex(numV,coords,neighborIndices, rotXtoY(subimage),uint8(1));
                else
                    [im,~,isConnected] = fillholes3D(numV,coords,neighborIndices, rotXtoY(subimage),uint8(1));
                end
                if isConnected
                    locations((i-1)*numberCutsPerDim(2)*numberCutsPerDim(3)+(j-1)*numberCutsPerDim(3)+k,2) = counter;
                    nameConnected = strcat('out_connected' ,num2str(counter));
                    nameOriginal = strcat('out_original' ,num2str(counter));
                    mat2rawrwd(reshape(1-im,MiniCubeSize),nameConnected, 1,[0,0,0], TotalSize./TotalSize(1), relPath);
                    mat2rawrwd(reshape(1 - rotXtoY(subimage),MiniCubeSize),nameOriginal, 1,[0,0,0], TotalSize./TotalSize(1), relPath);
                    counter = counter + 1;
                end

                if mexEnabled
                    [im,~,isConnected] = fillholes3D_mex(numV,coords,neighborIndices, rotXtoZ(subimage),uint8(1));
                else
                    [im,~,isConnected] = fillholes3D(numV,coords,neighborIndices, rotXtoZ(subimage),uint8(1));
                end
                if isConnected
                    locations((i-1)*numberCutsPerDim(2)*numberCutsPerDim(3)+(j-1)*numberCutsPerDim(3)+k,3) = counter;
                    nameConnected = strcat('out_connected' ,num2str(counter));
                    nameOriginal = strcat('out_original' ,num2str(counter));
                    mat2rawrwd(reshape(1-im, MiniCubeSize),nameConnected, 1,[0,0,0], TotalSize./TotalSize(1), relPath); 
                    mat2rawrwd(reshape(1 -rotXtoZ(subimage), MiniCubeSize),nameOriginal, 1,[0,0,0], TotalSize./TotalSize(1), relPath);
                    counter = counter + 1;
                end
            end
         end
     %end 
 end

end
 
  
 function out = rotXtoZ(in)
   out=permute(in,[3 2 1]);
    out=out(:);
 end
 
  function out = rotXtoY(in)
    out = permute(in,[2 1 3]);
        out=out(:);
  end
