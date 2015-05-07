classdef ResonantMixing < handle

    properties
        tend   = 1;
        dt     = 0.1;
        t      = 0;
        ic     = [];
        eps    = 0.03;
        b      = 0.03;
        dim    = 3;
        omega  = 1;
        psi    = 0;
        T      = 0;
        result = [];
        lam;
        V;
        PhaseFrames;
        TrajectoryFrames;
        PhysicalFrames;
    end
    
    methods
        function rm = ResonantMixing(Tend, Dt, E, B, W, IC)
            rm.tend  = Tend;
            rm.dt    = Dt;
            rm.t     = 0:rm.dt:rm.tend;
            rm.eps   = E;
            rm.b     = B;
            rm.omega = W;
            rm.ic    = IC;
            rm.dim   = length(IC);
        end
        function Simulate(obj)
            options = odeset('RelTol',1e-12,'AbsTol',1e-12);
            [~,y] = ode45(@obj.Flow, obj.t, obj.ic, options);
            obj.result = y;
        end
        function dy = EulerUpdate(obj, y)
            dy = y + obj.flow(obj.t,y)*obj.dt;
        end
        function dy = Flow(obj, t, y)
            dy = zeros(size(y));
            
            dy(1) = -cos(pi*y(1))*sin(pi*y(2)) + ...
                     pi*obj.b*sin(pi*y(1))*sin(pi*y(2))*sin(obj.omega*t) + ...
                     obj.eps*sin(2*pi*y(1))*sin(pi*y(3));
            dy(2) =  sin(pi*y(1))*cos(pi*y(2)) + ...
                     pi*obj.b*cos(pi*y(1))*cos(pi*y(2))*sin(obj.omega*t) + ...
                     obj.eps*sin(2*pi*y(2))*sin(pi*y(3));
            dy(3) = 2*obj.eps*cos(pi*y(3))*...
                    (cos(2*pi*y(1))+cos(2*pi*y(2)));
        end
        function J = Jacobian(obj, t0)
            id = find(t0-obj.dt/2 <= obj.t & obj.t < t0+obj.dt/2); 
            x  = obj.result(id,1);
            y  = obj.result(id,2);
            z  = obj.result(id,2);
            
            E  = obj.eps;
            B  = obj.b;
            W  = obj.omega;
            
            J = [pi*sin(pi*x)*sin(pi*y) + pi^2*B*cos(pi*x)*sin(pi*y)*sin(W*t0) + 2*pi*E*cos(2*pi*x)*sin(pi*z),...
                -pi*cos(pi*x)*cos(pi*y) + pi^2*B*sin(pi*x)*cos(pi*y)*sin(W*t0),...
                 pi*E*sin(2*pi*x)*cos(pi*z);... 
                 
                 pi*cos(pi*x)*cos(pi*y) - pi^2*B*sin(pi*x)*cos(pi*y)*sin(W*t0),...
                -pi*sin(pi*x)*sin(pi*y) - pi^2*B*cos(pi*x)*sin(pi*y)*sin(W*t0) + 2*pi*E*cos(2*pi*y)*sin(pi*z),...
                 pi*E*sin(2*pi*y)*cos(pi*z);...
                 
                -4*pi*E*cos(pi*z)*sin(2*pi*x),...
                -4*pi*E*cos(pi*z)*sin(2*pi*y),...
                -2*pi*E*sin(pi*z)*(cos(2*pi*x) + cos(2*pi*y));
                ];
        end
        function getPsi(obj)
            if isempty(obj.result)
                disp('error: You have to run the simulation first dummy');
                return;
            end
            obj.psi = cos(pi*obj.result(:,1)).*cos(pi*obj.result(:,2));
        end
        function getPeriod(obj)
            if isempty(obj.result)
                disp('error: You have to run the simulation first dummy');
                return;
            end
            obj.T = cos(pi*obj.result(:,1)).*cos(pi*obj.result(:,2));
        end
        function getEigen(obj)
            obj.lam = zeros(obj.dim,length(obj.t));
            obj.V = zeros(obj.dim,obj.dim,length(obj.t));
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
            label = 'xyz';
            for ii = 1:obj.dim
                subplot(obj.dim,1,ii)
                plot(obj.t(1:idx), obj.result(1:idx,ii),'k','linewidth', 3);
                grid on
                xlabel('Time','Interpreter', 'latex', 'fontsize',18)
                ylabel(label(ii),'Interpreter', 'latex', 'fontsize',18)
            
            end
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
            
            x = obj.result(1:idx,1);
            y = obj.result(1:idx,2);
            z = obj.result(1:idx,3);
            
            h = figure(fignum);
            plot3(x,y,z,'k','linewidth', 3);
            xlabel('$x$','Interpreter', 'latex', 'fontsize',18)
            ylabel('$y$','Interpreter', 'latex', 'fontsize',18)
            zlabel('$y$','Interpreter', 'latex', 'fontsize',18)
        end
        function h = PlotSlowPlane(obj,idx,fignum)
            if idx == 'All'
                    idx = length(obj.t);
            end
            
            figure(fignum)
            plot(obj.psi(1:idx), obj.result(1:idx,3),'k','linewidth',2);
            xlabel('$\psi$','interpreter','latex','fontsize',18)
            ylabel('$z$','interpreter','latex','fontsize',18)
            grid on
        end
        function timePlotEigen(obj, idx, fignum, transpose)
            
            if idx == 'All'
                idx = length(obj.t);
            end
            
            bgcolor= [0.8,0.8,0.8];
  
            
            figure(fignum)
            for ii = 1:obj.dim
                subplot(obj.dim,2,2*ii-1)
                plot(obj.t(1:idx),real(obj.lam(ii,1:idx)),'Color','k','LineWidth',3); 
                xlabel('Time','FontSize',18)
                ylabel(['$\lambda$',sprintf('%d',ii)],'interpreter','latex','FontSize',18)
                title('Real Part','FontSize',18)
                grid on
          
                ax1 = gca; % current axes
                ax1_pos  = ax1.Position;
           
                if strcmpi(transpose, 'on')
                    ax2 = axes('Position',ax1_pos,...
                               'YAxisLocation','right',...
                               'Color','none');
                    ax2.YColor = bgcolor;
                    h2 = line(obj.t(1:idx), obj.result(1:idx,ii),'Parent',ax2);
                    set(h2,'linewidth', 3, 'Color',bgcolor);
                end
        
                hold off;
            
            
                subplot(obj.dim,2,2*ii)

                plot(obj.t(1:idx),imag(obj.lam(ii,1:idx)),'k','LineWidth',3);
                title('Imaginary Part','FontSize',18)
                xlabel('Time','FontSize',18)
                ylabel(['$\lambda$',sprintf('%d',ii)],'interpreter','latex','FontSize',18)
                
                grid on
                ax1 = gca; % current axes
                ax1_pos  = ax1.Position;

                if strcmpi(transpose, 'on')
                    ax2 = axes('Position',ax1_pos,...
                               'YAxisLocation','right',...
                               'Color','none');
                    ax2.YColor = bgcolor;
                    h2 = line(obj.t(1:idx), obj.result(1:idx,ii));
                    set(h2,'linewidth', 3, 'Color',bgcolor);
                end
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
            
    end
end

