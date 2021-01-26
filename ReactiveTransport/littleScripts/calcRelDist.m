function out = calcRelDist(coarseGrid, coarseData, fineGrid, fineData, rangeL, rangeR)
%calculate L1 difference between P0 data on different meshes in the xrange
%[rangeL, rangeR]

grob = Variable(coarseGrid, Stepper(0:1), ...
    'Ca^(2+)', 'P0');
grob.setdata(coarseData);
fein = Variable(fineGrid, Stepper(0:1), ...
    'Ca^(2+)', 'P0');
fein.setdata(fineData);
wfun = grob.getTSI(0);
out = distanceL1(fineGrid, fineData, @(x) wfun(x(1), x(2)), rangeL, rangeR); %/distanceL1(fineGrid, fineData, @(x) 0,rangeL, rangeR);

end