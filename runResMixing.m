function rm = runResMixing() 
    Tend  = 100;
    Dt    = 0.01;
    eps   = 0.03;
    b     = 0.0;
    omega = 0;
    IC    = [0.1,0,0];

    rm = ResonantMixing(Tend, Dt, eps, b, omega, IC);
    
    rm.Simulate();
    rm.getEigen();
    rm.getPsi();
    rm.PlotTrajectory('All', 1);
    rm.PlotSlowPlane('All', 2);
    rm.PlotPhysicalSpace('All', 3);
    rm.timePlotEigen('All', 4, 'on');
    %sp.makeMovie(6, 'phase');
    %sp.makeMovie(7, 'physicalspace');
    %SaveVid(sp.PhaseFrames,'./Results/PhasePlane2.avi')
    %SaveVid(sp.PhysicalFrames,'./Results/PhysicalPlane2.avi')
end