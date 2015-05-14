classdef Kuramoto < handle

    properties
        tend   = 1;
        dt     = 0.1;
        t      = 0;
        ic     = [];
        dim    = 3;
        w      = 1;
        K      = 1;
        result = [];
        lam;
        V;
        CircleFrames;
        TrajectoryFrames;
    end
    
    methods
        function k = Kuramoto(Tend, Dt, W, K, IC)
            k.tend  = Tend;
            k.dt    = Dt;
            k.t     = 0:k.dt:k.tend;
            k.w     = W;
            k.K     = K;
            k.ic    = IC;
            k.dim   = length(IC);
        end
        
        %% ODE solve
        function Simulate(obj)
            options = odeset('RelTol',1e-12,'AbsTol',1e-12);
            [~,y] = ode45(@obj.Flow, obj.t, obj.ic, options);
            obj.result = y;
        end
        function dy = EulerUpdate(obj, y)
            dy = y + obj.flow(obj.t,y)*obj.dt;
        end
        function dy = Flow(obj, ~, y)
            dy = zeros(obj.dim,1);
            
            for ii = 1:obj.dim
                dy(ii) = obj.w(ii) + obj.K(ii)/obj.dim * sum(sin(y - y(ii)));
            end
            
        end
        
        %% Obtaining Eigenvalues
        function J = Jacobian(obj, t0)
            id = find(t0-obj.dt/2 <= obj.t & obj.t < t0+obj.dt/2); 
            x  = obj.result(id,1);
            y  = obj.result(id,2);
            z  = obj.result(id,2);
            J = zeros(obj.dim);
            
            for ii = 1:obj.dim
                for jj = 1:obj.dim
                    thetai = obj.result(id,ii);
                    thetaj = obj.result(id,jj);
                    thetas = obj.result(id,:);
                    if ii == jj 
                        J(ii,jj) = 1 - obj.K(ii)/obj.dim*sum(cos(thetas - thetai));
                    else
                        j(ii,jj) = obj.K(ii)/obj.dim*cos(thetaj - thetai);
                    end
                end
            end
            
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
        
        %% Plotting Functions
        function h = PlotTrajectory(obj, idx, fignum)
            if idx == 'All'
                idx = length(obj.t);
            end
                   
            h = figure(fignum);
            plot(obj.t(1:idx), sin(obj.result(1:idx,:)),'linewidth', 3);
            grid on
            xlabel('Time','Interpreter', 'latex', 'fontsize',18)
            ylabel('Agents','Interpreter', 'latex', 'fontsize',18)
        end
        function h = PlotOnCircle(obj,idx,fignum)
           if idx == 'All'
                idx = length(obj.t);
           end
           spread = linspace(0.85,1.15,obj.dim);
           theta = linspace(0, 2*pi, 100);
           h = figure(fignum);
           plot(cos(theta),sin(theta),'k'); hold on;
           p = plot(spread.*cos(obj.result(idx,:)), spread.*sin(obj.result(idx,:)),'o');
           set(p,'Marker','o', 'MarkerFaceColor','k', 'MarkerEdgeColor','k', 'MarkerSize',8)
           hold off
           grid on
           xlabel('$X$','Interpreter', 'latex', 'fontsize',18)
           ylabel(['$Y$'],'Interpreter', 'latex', 'fontsize',18)
           axis([-1.2,1.2,-1.2,1.2])
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
                    h2 = line(obj.t(1:idx), sin(obj.result(1:idx,ii)),'Parent',ax2);
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
                if strcmpi(type,'circle')
                    h = obj.PlotOnCircle(ii,fignum);
                elseif strcmpi(type, 'trajectory')
                    h = obj.PlotTrajectory(ii,fignum);
                end    
                frames(ii) = getframe(h);
                close(h)
            end
            if strcmpi(type,'circle')
                obj.CircleFrames = frames;
            elseif strcmpi(type, 'trajectory')
                obj.TrajectoryFrames = frames;
            end
        end
            
    end
end

