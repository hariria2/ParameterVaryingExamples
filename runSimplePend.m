function runSimplePend() 
    Tend = 25;
    Dt   = 0.01;
    L    = @(t) abs(0.1+sin(0.1*t)); %1+0.2*t;
    M    = @(t) ones(size(t));
    B    = @(t) zeros(size(t)); %1.5;
    IC   = [pi/6; 4.7];

    %L    = @(t) ones(size(t)); %0.5+0.5*t;
    %M    = @(t) ones(size(t));
    %B    = @(t) zeros(size(t));%exp(-3*t);
    %IC   = [pi/6; 0];

    sp = SimplePendulum(Tend, Dt, L, M, B, IC);
    %B    = @sp.getB;
    sp.b = B;
    sp.Simulate();
    sp.getEigen();

    sp.PlotTrajectory('All', 1);
    sp.PlotPhase('All',2);
    sp.PlotPhysicalSpace('All', 3);

    sp.timePlotEigen('All', 4, 'on');
    sp.argonPlotEigen();
    %sp.makeMovie(6, 'phase');
    %sp.makeMovie(7, 'physicalspace');
    sp.PlotParameters('All', 8);

    %SaveVid(sp.PhaseFrames,'./Results/PhasePlane2.avi')
    %SaveVid(sp.PhysicalFrames,'./Results/PhysicalPlane2.avi')
end