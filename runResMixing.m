function rm = runResMixing() 
    Tend  = 100;
    Dt    = 0.001;
    eps   = 0.03;
    b     = 0.003;
    omega = 2.5;
    IC    = [0.3,0.3,-0.3];
    
    rm = ResonantMixing(Tend, Dt, eps, b, omega, IC);
    
    rm.Simulate();
    rm.getEigen();
    rm.getPsi();
    rm.PlotTrajectory('All', 1);
    rm.PlotSlowPlane('All', 2);
    rm.PlotPhysicalSpace('All', 3);
    rm.timePlotEigen('All', 4, 'off');
    %sp.makeMovie(6, 'phase');
    %sp.makeMovie(7, 'physicalspace');
    %SaveVid(sp.PhaseFrames,'./Results/PhasePlane2.avi')
    %SaveVid(sp.PhysicalFrames,'./Results/PhysicalPlane2.avi')
end