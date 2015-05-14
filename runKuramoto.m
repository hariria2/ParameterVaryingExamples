function k = runKuramoto() 
    Tend  = 15;
    Dt    = 0.01;
    N     = 50;
    W     = abs(2*pi + 1*randn(1,N));
    K     = abs(40 + 10*randn(1,N));
    IC    = 2*pi*rand(1,N)-pi;

    k = Kuramoto(Tend, Dt, W, K, IC);
    
    k.Simulate();
    k.PlotTrajectory('All', 1);
    k.PlotOnCircle(1, 2);
    %k.getEigen();
    %k.timePlotEigen('All',3,'on')
    %k.makeMovie(3,'circle');
    %SaveVid(k.CircleFrames,'./Results/Kuramoto6.avi')
end