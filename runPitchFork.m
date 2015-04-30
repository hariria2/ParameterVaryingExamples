function pf = runPitchFork() 
    Tend  = 2;
    Dt    = 0.01;
    b     = @(t) t-1;
    IC    = [0.1];

    pf = PitchForkBifurc(Tend, Dt, b, IC);
    
    pf.Simulate();
    pf.getEigen();

    pf.PlotTrajectory('All', 1);
    pf.timePlotEigen('All', 4, 'on');
    %sp.makeMovie(6, 'phase');
    %sp.makeMovie(7, 'physicalspace');
    %SaveVid(sp.PhaseFrames,'./Results/PhasePlane2.avi')
    %SaveVid(sp.PhysicalFrames,'./Results/PhysicalPlane2.avi')
end