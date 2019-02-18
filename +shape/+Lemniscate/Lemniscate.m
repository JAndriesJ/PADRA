classdef Lemniscate< shape.C2boundary
    % Class for regular n-polygon shape 
    properties(SetAccess = protected)
        foci = [] ;
        radius = 1 ;
        scale  = 1 ; % the scale of the 
		phi    = 0 ; % orientation of the domain
        degree = 0 ;
    end
    
    methods
      function obj = Lemniscate(foci,nbPoints, dspl)
        if nargin < 1 || isempty(foci)
            nbFoci = randi(2)+1;   
            piRange = pi*1 ; 
            foci = (0.5+0.5*rand)*[cos(piRange*rand(1,nbFoci));0.5+sin(piRange*rand(1,nbFoci))];
            centerFocus = [mean(foci(1,:));mean(foci(2,:))] ;
            foci = foci -centerFocus + 0.25*rand(2,1)+[0.25;0.25] ;
            nbPoints = 1000 ;
            dspl = 1 ;
        elseif nargin < 2 
            nbPoints = 1000 ;
            dspl = 1 ;
        elseif nargin < 3
            dspl = 1 ;
        end

        [radius,poly] =  shape.Lemniscate.Lemniscate.genpoly(foci) ;
            
        %Initialize
        ODEfun =  shape.algdom.getHamiltonian(poly); 
        tstep = 1/nbPoints    ;
        tspan = (0:tstep:9)   ;  


        function [value, isterminal, direction] = EventsFcn1(Te,Y)
        value      =   10*(max(abs(Y(1))+abs(Y(2))) - (10^(-1))) + (Te < tstep*(5*10^2))  ; % Minimizer
        isterminal =   1         ; % Halt integration 
        direction  =   0         ; % approach zero from negative side
        end  

        options1 = odeset('Events',@EventsFcn1,'RelTol',1e-12);
        %%ODE solver                              
        [T,Y] = ode113(ODEfun,tspan,[0;0],options1);     
        meanPointDist = sum(diff(vecnorm(Y(1:50,:)')))/50;
        lastbitofT =  norm(Y(end,:))/meanPointDist ; 
        meanTimeDist = T(50)/50 ;   
        tend = T(end) + meanTimeDist*lastbitofT;
        
%         OrbitalNormalizer =  2*pi/tend ;
        options = odeset('RelTol',1e-12) ;            
        tspan = 0:(tend/nbPoints):tend ;

        % %ODE solver                              
        [~,z] = ode113(ODEfun,tspan,[0;0],options); 
        %points
        points = z(1:end-1,:)';
        points = circshift(points,[0,3*(floor(nbPoints/(14/5)))]);
            %% velocity 
%             VeloFun = shape.algdom.getVelocity(poly,OrbitalNormalizer) ;
%             tvec = zeros(2,length(points));
%             for iter = 1:length(points)
%                 tvec(:,iter) =VeloFun(points(:,iter)) ;
%             end
 
            %% normal vector
%             normal = [[0 1];[-1 0]]*tvec ;
% 			tt = sqrt(normal(1,:).^2+normal(2,:).^2);
% 			normal = [normal(1,:)./ tt; normal(2,:)./ tt];       
            
            %% acceleration
%             AccFun =  shape.algdom.getAccelleration(poly,OrbitalNormalizer) ;
%             avec =  zeros(2,length(points));
%             for iter = 1:length(points)
%                 avec(:,iter)  = AccFun(points(:,iter)) ; 
%             end
  
            %% Cop-out approach          
            theta = 0:(1/length(points)):(1-(1/length(points))) ;	
            [points, tvec, avec, normal] = shape.C2boundary.rescale(points, theta, nbPoints, [], dspl);            
            %% object creation  
			obj= obj@shape.C2boundary(points, tvec, avec, normal, [0;0], 'Lemniscate');
            obj.foci = foci;
            obj.radius = radius ; 
            obj.degree = 2*size(foci,2);
      end

        
       function obj = mtimes(obj, s)
			% Overload of the operator *
			obj = mtimes@shape.C2boundary(obj, s);
			obj.scale  = obj.scale*s ; 
            obj.coefvec =  shape.algdom.Transpoly([s, 0; 0, s], obj.coefvec) ;
       end
		
	   function obj = lt(obj, phi)
			% Redefine the < operator as the rotation of the boundary
			obj = lt@shape.C2boundary(obj,phi);
			obj.phi = obj.phi+phi;
        end
     end        
        methods (Static)   
        % polynomail
        function [radius,poly] = genpoly(foci) 
        syms  x y
        poly = 1 ;
        radius = prod(foci(1,:).^2 + foci(2,:).^2) ;
        
        for iter1 = 1:size(foci,2)
            poly = poly*((x-foci(1,iter1))^2 + (y-foci(2,iter1))^2);             
        end
        poly = poly - radius ;
        end
        
        % Hamiltonian
        function Hamiltonian = getHamiltonian(poly)
        syms x y t
        GradPoly   = [diff(poly,x);diff(poly,y)];
        RotGrad   = [0,-1;1,0]*GradPoly;                                     
        fRotGrad = matlabFunction(RotGrad,'Vars',[t x y]);
        Hamiltonian = @(t,z) fRotGrad(t,z(1),z(2));
        end
        
        % Velocity
        function Velocity = getVelocity(poly,OrbitalNormalizer)
        syms x y t    
        RotGradientPoly   = [-diff(poly,y) ; diff(poly,x)]./OrbitalNormalizer;
        Velocity = matlabFunction(RotGradientPoly,'Vars',[x y]) ;
        Velocity = @(z) Velocity(z(1),z(2)) ;
        end
        
        % Acceleration
        function Accelleration = getAccelleration(poly,OrbitalNormalizer)
        syms x y            
        d = double(feval(symengine, 'degree', poly, x)) ; % obtaining degree of polynomial

        Hessian = [-diff(diff(poly,y),x),-diff(diff(poly,y),y);
                    diff(diff(poly,x),x), diff(diff(poly,x),y)];
        Velocity = shape.algdom.getVelocity(poly,OrbitalNormalizer) ;

        if double(d)/2 > 1
            Accelleration  = matlabFunction(Hessian./OrbitalNormalizer ,'Vars',[x y] );
            Accelleration = @(z) Accelleration(z(1),z(2))*Velocity(z) ;                      
        else                    
            Accelleration  = matlabFunction(Hessian./OrbitalNormalizer,'Vars',[x y] );
            Accelleration = @(z) Accelleration(z(1),z(2))*Velocity(z)  ;  
        end
  
        end
        
            
    end
end