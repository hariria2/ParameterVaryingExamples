function k = runKuramoto() 
    Tend  = 5;
    Dt    = 0.05;
    N     = 20;
    W     = abs(2*pi + pi*randn(1,N));
    K     = abs(10 + 0*randn(1,N));
    IC    = 2*pi*rand(1,N)-pi;

    k = Kuramoto(Tend, Dt, W, K, IC);
    
    k.Simulate();
    k.PlotTrajectory('All', 1);
    k.PlotOnCircle(1, 2);
    k.getEigen();
    k.timePlotEigen('All',3,'on')
    k.makeMovie(3,'circle');
    SaveVid(k.CircleFrames,'./Results/Kuramoto6.avi')
end