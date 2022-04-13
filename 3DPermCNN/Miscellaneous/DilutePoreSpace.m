% Inflate or deflate pore-space using level-set algorithms
% Requires RTSPHEM/LevelSetTriplePoint to be on file path

p=genpath('./RawConversion');
addpath(p);
p=genpath('./VoxelData/SampleBENTHEIMER');
addpath(p);
p=genpath('./Preprocessing');
addpath(p);


% setup underlying grid

MiniCubeSize = [100,100,100];
grid=CartesianGrid3D(3,[0 1 0 1 0 1],MiniCubeSize-1);
numV = uint32(grid.nodes);
coords = grid.coordinates;
neighborIndices = uint32([grid.pull{1, 1, 1}, grid.pull{2, 1, 1}, grid.pull{1, 1, 2}, grid.pull{2, 1, 2},grid.pull{1, 1, 3}, grid.pull{2, 1, 3}]);
relPath = './VoxelData/VaryPoro/';

counter = 0;

    % read in file
	
    nameRead = strcat('out_connected' ,num2str(0),'_100x100x100.rwd');
    A=rwd2mat(nameRead);
 
    % turn binary image into level-set representation

    levelSetInflate = reinitializeLevelSet( grid, double(A(:))-0.5, true, inf );
    levelSetDeflate = levelSetInflate;
    
    % Compute sequence of successively deflating the pore space until sample becoms impermeable    

     for i=1:10000
        levelSetDeflate = levelSetEquationTimeStep( 0.005, 0, levelSetDeflate, ...
        grid, 1 );												%one time step level set equation --> move boundary by 0.005
        nameWrite = strcat('VaryPoro' ,num2str(counter));
        [im,~,isConnected] = fillholes3D(numV,coords,neighborIndices, 1-uint8(levelSetDeflate>0),uint8(1));   	%eliminate dead pores
        if ~isConnected
            break
        end
        mat2rawrwd(reshape(1-im,MiniCubeSize),nameWrite,1,[0,0,0],[1,1,1],relPath);					%write file
        sum(1-im(:)>0)/prod(MiniCubeSize)
        counter=counter+1;
    end

    % Compute sequence of successively inflate the pore space, same number of step as before     

    for i=counter:10000
    	levelSetInflate = levelSetEquationTimeStep( 0.005, 0, levelSetInflate, ...
        grid, -1 );												%one time step level set equation --> move boundary by 0.005
        nameWrite = strcat('VaryPoro' ,num2str(i));
        [im,~,~] = fillholes3D(numV,coords,neighborIndices, 1-uint8(levelSetInflate>0),uint8(1));  		%eliminate dead pores      
        mat2rawrwd(reshape(1-im,MiniCubeSize),nameWrite,1,[0,0,0],[1,1,1],relPath);					%write file
        sum(1-im(:)>0)/prod(MiniCubeSize)
        if 2*counter == (i+1)
            break
        end
    end
    
   
    
    
