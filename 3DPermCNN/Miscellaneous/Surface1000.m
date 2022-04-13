% Script to compute the interior surface area of large pore spaces
% Add contributions from 'subsamples' that can be computed independently

%Berea 355.4878
%Bentheimer 
%Castle 334.6657
p=genpath('./RawConversion');
addpath(p);

bigCubeSize = [1000,1000,1000];												%Size of sample of interst
imSize = [100,100,100];														%Size of subsamples
shift=1;
[X,Y,Z] = meshgrid(1:imSize(1),1:imSize(2),1:imSize(3));
input=rwd2mat('Bentheimer_2d25um_binary_1000x1000x1000.rwd');				%Sample read-in file
area=cell(bigCubeSize(1)./imSize(1));
[area{:}]= deal(zeros(prod(bigCubeSize(2:3)./imSize(2:3)),1));

parfor a=1:(bigCubeSize(1)/imSize(1))
	for b=1:(bigCubeSize(2)/imSize(2))
		for c=1:(bigCubeSize(3)/imSize(3))

			linInd = c+(b-1)*10+(a-1)*imSize(2)*imSize(3);
			A=1-input(((a-1)*shift)+(1:imSize(1)), ((b-1)*shift)+(1:imSize(2)),((c-1)*shift)+(1:imSize(3)));       	%Cut subsample 
			out = isosurface(X,Y,Z,A);										%triangulate pore-space walls
			B = zeros(size(out.faces,1), 3,3);
			for j=1:size(out.faces,1)										%compute triangle surface areas
   				B(j,:,:) = out.vertices(out.faces(j,:),:);
			end

			B(:,2,:) = B(:,2,:) - B(:,1,:);
			B(:,3,:) = B(:,3,:) - B(:,1,:);
			area{a}( c+(b-1)*10 ) = 0.5*sum(vecnorm(cross(squeeze(B(:,2,:)), squeeze(B(:,3,:)),2),2,2));
		end
	end
end

out = cell2mat(area);														%add triangle surfaces and rescale to [mm^2]
out=sum(out(:))*(2.25*10^(-3)).^2
