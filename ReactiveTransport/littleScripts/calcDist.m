%auxiliary file
%Calculate norm between two refinement levels within certain x range

function out = calcDist(coarseGrid, coarseData, fineGrid, fineData, rangeL, rangeR)

coarse = Variable(coarseGrid, Stepper(0:1), ...
    'Ca^(2+)', 'P0');
coarse.setdata(coarseData);
fine = Variable(fineGrid, Stepper(0:1), ...
    'Ca^(2+)', 'P0');
fine.setdata(fineData);
wfun = coarse.getTSI(0);
out = distanceL1(fineGrid, fineData, @(x) wfun(x(1), x(2)), rangeL, rangeR); %/distanceL1(fineGrid, fineData, @(x) 0,rangeL, rangeR);

end