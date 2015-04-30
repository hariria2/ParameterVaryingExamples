classdef SimplePendulum < handle

    properties
        tend   = 1;
        dt     = 0.1;
        t      = 0;
        ic     = [];
        l      = @(t) 1;
        b      = @(t) 1;
        m      = @(t) 1;
        result = [];
        g      = 9.81;
        lam;
        V;
        PhaseFrames;
        TrajectoryFrames;
        PhysicalFrames;
    end
    
    methods
        function sp = SimplePendulum(Tend, Dt, L, M, B, IC)
            sp.tend = Tend;
            sp.dt   = Dt;
            sp.t    = 0:sp.dt:sp.tend;
            sp.l    = L;
            sp.m    = M;
            sp.b    = B;
            sp.ic   = IC;
        end
        function r = getB(~, t)
            r = zeros(size(t));
            for ii = 1:length(t)
               if t(ii)<10
                   r(ii) = 0;
               else
                   r(ii) = 5;
               end
            end
        end
        function Simulate(obj)
            [~,y] = ode45(@obj.Flow, obj.t, obj.ic);
            obj.result = y;
        end
        function dy = EulerUpdate(obj, y)
            dy = y + obj.flow(obj.t,y)*obj.dt;
        end
        function dy = Flow(obj, t, y)
            dy = zeros(size(y));
            dy(1) = y(2);
            dy(2) = -obj.b(t)*y(2) - obj.g/obj.l(t) * sin(y(1));
        end
        function J = Jacobian(obj, t0)
            id = find(t0-obj.dt/2 <= obj.t & obj.t < t0+obj.dt/2); 
            x1  = obj.result(id,1);
            x2  = obj.result(id,2);
            
            L  = obj.l(t0);
            B  = obj.b(t0);
            %M  = obj.m(t0);
            
            J = [               0, 1;...
                 -obj.g/L*cos(x1), -B];
        end 
        function getEigen(obj)
            obj.lam = zeros(2,length(obj.t));
            obj.V = zeros(2,2,length(obj.t));
            for ii = 1:length(obj.t)
                J = obj.Jacobian(obj.t(ii));
                [v,d] = eig(J);
                obj.lam(:,ii) = diag(d);
                obj.V(:,:,ii) = v;
                obj.lam(abs(obj.lam)<1e-5)=0;
            end
        end
        function h = PlotTrajectory(obj, idx, fignum)
            
            if idx == 'All'
                idx = length(obj.t);
            end
            
            h = figure(fignum);
            subplot(2,1,1)
            plot(obj.t(1:idx), obj.result(1:idx,1),'k','linewidth', 3);
            grid on
            xlabel('Time','Interpreter', 'latex', 'fontsize',18)
            ylabel('$\theta$','Interpreter', 'latex', 'fontsize',18)
            ax = gca;
            %ax.YTickLabel = { '-3\pi','-2\pi','-\pi','0','\pi','2\pi', '3\pi'};
            
            subplot(2,1,2)
            plot(obj.t(1:idx), obj.result(1:idx,2),'k','linewidth', 3);
            grid on
            xlabel('Time','Interpreter', 'latex', 'fontsize',18)
            ylabel('$\dot{\theta}$','Interpreter', 'latex', 'fontsize',18)
            ax = gca;
            %ax.YTickLabel = {'-3\pi','-2\pi','-\pi','0','\pi','2\pi','3\pi', '4\pi'};
        end
        function h = PlotPhase(obj, idx, fignum)
            if idx == 'All'
                idx = length(obj.t);
            end
            
            [x,y] = meshgrid(-5:0.5:15,-10:0.5:10);
            u = y;
            v = -obj.b(obj.t(idx)) - obj.g/obj.l(obj.t(idx)) * sin(x);
            
            h = figure(fignum);
            quiver(x,y,u,v,2,'k'); hold on
            
            plot(obj.result(1:idx,1), obj.result(1:idx,2), 'k','linewidth', 3);
            grid on
            xlabel('$\theta$','Interpreter', 'latex', 'fontsize',18)
            ylabel('$\dot{\theta}$','Interpreter', 'latex', 'fontsize',18)
            %set(gcf,'units','normalized','outerposition',[0 0 1 1])
        end
        function h = PlotPhysicalSpace(obj, idx, fignum)
            
            if idx == 'All'
                idx = length(obj.t);
            end
            
            x = obj.l(obj.t(1:idx)).*sin(obj.result(1:idx,1)');
            y = obj.l(obj.t(1:idx)).*(-cos(obj.result(1:idx,1)))'; % - max(obj.l(obj.t));
            
            h = figure(fignum);
            plot(0,0,'ko','markersize', 9); hold on
            line([0,x(end)],[0, y(end)],'Color', 'k'); hold on
            plot(x(1), y(1),'ko','markersize', 9); hold on
            plot(x(end), y(end),'kx','markersize', 9); hold on
            plot(x, y, 'k','linewidth', 3);
            xlim([-max(obj.l(obj.t)),max(obj.l(obj.t))]);
            ylim([-max(obj.l(obj.t)),max(obj.l(obj.t))]);
            grid on
            xlabel('$x$','Interpreter', 'latex', 'fontsize',18)
            ylabel('$y$','Interpreter', 'latex', 'fontsize',18)
            hold off
        end
        function timePlotEigen(obj, idx, fignum, transpose)
            
            if idx == 'All'
                idx = length(obj.t);
            end
            
            bgcolor= [0.8,0.8,0.8];
                   
            lam1 = obj.lam(1,1:idx);
            lam2 = obj.lam(2,1:idx);
            
            figure(fignum)
            subplot(2,2,1)
            h1 = plot(obj.t(1:idx),real(lam1),'Color','k','LineWidth',3); 
            xlabel('Time','FontSize',18)
            ylabel('Real Part','FontSize',18)
            title('\lambda_1','FontSize',18)
            grid on
          
            ax1 = gca; % current axes
            ax1_pos  = ax1.Position;
            
            
            if strcmpi(transpose, 'on')
                ax2 = axes('Position',ax1_pos,...
                           'YAxisLocation','right',...
                           'Color','none');
                ax2.YColor = bgcolor;
                %ax2.YTickLabel = {'-5\pi', '-4\pi', '-3\pi','-2\pi','-\pi','0','\pi','2\pi','3\pi', '4\pi', '5\pi'};
                h2 = line(obj.t(1:idx), obj.result(1:idx,1),'Parent',ax2);
                set(h2,'linewidth', 3, 'Color',bgcolor);
            end
        
            hold off;
            
            subplot(2,2,2)
            
            h1 = plot(obj.t(1:idx),real(lam2),'k','LineWidth',3);
            xlabel('Time','FontSize',18)
            ylabel('Real Part','FontSize',18)
            title('\lambda_2','FontSize',18)
            grid on
            ax1 = gca; % current axes
            ax1_pos  = ax1.Position;
            if strcmpi(transpose, 'on')
                ax2 = axes('Position',ax1_pos,...
                           'YAxisLocation','right',...
                           'Color','none');
                ax2.YColor = bgcolor;
                %ax2.YTickLabel = {'-5\pi', '-4\pi', '-3\pi','-2\pi','-\pi','0','\pi','2\pi','3\pi', '4\pi', '5\pi'};
                h2 = line(obj.t(1:idx), obj.result(1:idx,2), 'Parent', ax2);
                set(h2,'linewidth', 3, 'Color',bgcolor);
            end
            
            subplot(2,2,3)
            
            h1 = plot(obj.t(1:idx),imag(lam1),'k','LineWidth',3);
            xlabel('Time','FontSize',18)
            ylabel('Imaginary Part','FontSize',18)
            grid on
            ax1 = gca; % current axes
            ax1_pos  = ax1.Position;
            
            if strcmpi(transpose, 'on')
                ax2 = axes('Position',ax1_pos,...
                           'YAxisLocation','right',...
                           'Color','none');
                ax2.YColor = bgcolor;
                %ax2.YTickLabel = {'-5\pi', '-4\pi', '-3\pi','-2\pi','-\pi','0','\pi','2\pi','3\pi', '4\pi', '5\pi'};
                h2 = line(obj.t(1:idx), obj.result(1:idx,1));
                set(h2,'linewidth', 3, 'Color',bgcolor);
            end
            
            subplot(2,2,4)
            
            h1 = plot(obj.t(1:idx),imag(lam2),'k','LineWidth',3);
            xlabel('Time','FontSize',18)
            ylabel('Imaginary Part','FontSize',18)
            grid on
            ax1 = gca; % current axes
            ax1_pos  = ax1.Position;
            if strcmpi(transpose, 'on')
                ax2 = axes('Position',ax1_pos,...
                           'YAxisLocation','right',...
                           'Color','none');
                ax2.YColor = bgcolor;
                %ax2.YTickLabel = {'-5\pi', '-4\pi', '-3\pi','-2\pi','-\pi','0','\pi','2\pi','3\pi', '4\pi', '5\pi'};
                h2 = line(obj.t(1:idx), obj.result(1:idx,2)); 
                set(h2,'linewidth', 3, 'Color',bgcolor);
            end
        end
        function argonPlotEigen(obj)
            lam1 = obj.lam(1,:);
            lam2 = obj.lam(2,:);
            
            figure
           
            plot(real(lam1),imag(lam1),'k.',...
                 real(lam2),imag(lam2),'b.')
            grid on
            title('\lambda_1 and \lambda_2','FontSize',20)
            xlabel('Real Part','FontSize',18)
            ylabel('Imaginary Part','FontSize',18)
            leg = legend('\lambda_1','\lambda_2');
            set(leg,'FontSize',18)
        end
        function makeMovie(obj,fignum,type, varargin)
            frames(length(obj.t)) = struct('cdata',[],'colormap',[]);
            for ii = 1:length(obj.t)
                if strcmpi(type,'phase')
                    h = obj.PlotPhase(ii,fignum);
                elseif strcmpi(type, 'trajectory')
                    h = obj.PlotTrajectory(ii,fignum);
                elseif strcmpi(type, 'physicalspace')
                    h = obj.PlotPhysicalSpace(ii,fignum);
                end    
                frames(ii) = getframe(h);
                close(h)
            end
            if strcmpi(type,'phase')
                obj.PhaseFrames = frames;
            elseif strcmpi(type, 'trajectory')
                obj.TrajectoryFrames = frames;
            elseif strcmpi(type, 'physicalspace')
                obj.PhysicalFrames = frames;
            end
        end
        function h = PlotParameters(obj, idx, fignum)
          if idx == 'All'
                idx = length(obj.t);
          end
          
          h = figure(fignum);
          suptitle('Parameters')
          subplot(3,1,1)
          plot(obj.t(1:idx),obj.l(obj.t(1:idx)),'k','linewidth',3)
          xlabel('Time', 'fontsize', 18)
          ylabel('Length', 'fontsize', 18)
          
          subplot(3,1,2)
          plot(obj.t(1:idx),obj.b(obj.t(1:idx)),'k','linewidth',3)
          xlabel('Time', 'fontsize', 18)
          ylabel('Damping Coefficient', 'fontsize', 18)
          
          subplot(3,1,3)
          plot(obj.t(1:idx),obj.m(obj.t(1:idx)),'k','linewidth',3)
          xlabel('Time', 'fontsize', 18)
          ylabel('Mass', 'fontsize', 18)
          
        end
            
    end
end

