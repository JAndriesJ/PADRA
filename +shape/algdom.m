classdef algdom < shape.C2boundary
	 
	properties(SetAccess = protected)
        scale  = 1 ; % the scale of the 
		phi    = 0 ; % orientation of the domain
        coefvec = [] ;
        isLoop = false ; 
    end
     
	methods
		function obj =algdom(Coefvec,nbPoints, dspl)
          if nargin < 3
              dspl = 1 ;
          end            
            poly = shape.algdom.genpoly(Coefvec)     ;  
            
            %Initialize
            ODEfun =  shape.algdom.getHamiltonian(poly); 
            tstep = 1/nbPoints  ;
            tspan = (0:tstep:10)   ;  
            %FirstStep = norm(ODEfun(tstep^2,[0,0]),2) ; % ,FirstStep
            
            function [value, isterminal, direction] = EventsFcn1(Te,Y)
            value      =   abs(Y(1))+abs(Y(2)) - 0.5*(10^(-1)) + (Te < tstep*(10^2))   ; % Minimizer
            isterminal =   1         ; % Halt integration 
            direction  =   0         ; % approach zero from negative side
            end  
                   
                                       
            options1 = odeset('Events',@EventsFcn1,'RelTol',1e-10 ) ;
            %%ODE solver                              
            [T,Y] = ode113(ODEfun,tspan,[0;0],options1) ;
            tend = T(end) ;
  
            %OrbitalNormalizer =  2*pi/tend ;
                                
            options = odeset('RelTol',1e-12) ;            
            tspan = 0:(tend/nbPoints):tend ;
            
            % %ODE solver                              
            [~,z] = ode113(ODEfun,tspan,[0;0],options); 
            %points
            points = z';
                        
%             % velocity 
%             VeloFun = shape.algdom.getVelocity(poly,OrbitalNormalizer) ;
%             tvec = zeros(2,length(points));
%             for iter = 1:length(points)
%                 tvec(:,iter) =VeloFun(points(:,iter)) ;
%             end
%  
%             % normal vector
%             normal = [[0 1];[-1 0]]*tvec ;
% 			tt = sqrt(normal(1,:).^2+normal(2,:).^2);
% 			normal = [normal(1,:)./ tt; normal(2,:)./ tt];       
%             
%             % acceleration
%             AccFun =  shape.algdom.getAccelleration(poly,OrbitalNormalizer) ;
%             avec =  zeros(2,length(points));
%             for iter = 1:length(points)
%                 avec(:,iter)  = AccFun(points(:,iter)) ; 
%             end
            %% Cop-out approach
            %points = circshift(points,[0,floor(nbPoints/3)])    ; %this is so that the begining and end are centered on an edge
            theta = 0:(1/length(points)):(1-(1/length(points))) ;	
            [points, tvec, avec, normal] = shape.C2boundary.rescale(points, theta, nbPoints, [], dspl);
            
            
            %% object creation  
			obj= obj@shape.C2boundary(points, tvec, avec, normal, [0;0], 'Algebraic Domain');
            obj.scale = 1;
            obj.phi = 0;
            obj.coefvec = Coefvec ;
            
            if norm(points(:,1) - points(:,end)) < 0.1
                obj.isLoop = true ;
            end
            
        end
        
        %scaling
        function obj = mtimes(obj, s)
			% Overload of the operator *
			obj = mtimes@shape.C2boundary(obj, s);
			obj.scale  = obj.scale*s ; 
            obj.coefvec =  shape.algdom.Transpoly([s, 0; 0, s], obj.coefvec) ;
        end
        
		%rotation
		function obj = lt(obj, phi)
			% Redefine the < operator as the rotation of the boundary
			obj = lt@shape.C2boundary(obj,phi);
			obj.phi = obj.phi+phi;
            obj.coefvec =  shape.algdom.Transpoly([cos(phi), -sin(phi); sin(phi), cos(phi)], obj.coefvec) ;
        end
    end
    
    methods (Static)   
        % polynomail
        function poly = genpoly(Coefvec) 
        syms x y 
        % acceptable coefvec sizes
        aclengs = zeros(1,10);
        xexpo = length(Coefvec);
        yexpo = length(Coefvec);
        count = 0 ;
        for iter1 = 2:10
            aclengs(iter1) = iter1 + aclengs(iter1-1);
            for iter2 = 0:(iter1-1)
                count = count + 1 ;
                xexpo(count) = iter1-iter2-1;
                yexpo(count) = iter2 ;
            end
        end

        if ~any(length(Coefvec) == aclengs)|| size(Coefvec,1) ~= 1
            error('Not appropriate coefficient vector size')
        else
            xexpo = xexpo(1:length(Coefvec)) ;
            yexpo = yexpo(1:length(Coefvec)) ;
            poly = sum((x.^xexpo).*(y.^yexpo).*(Coefvec)) ;
        end
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
        
        % transform polynomial
        function newCoef = Transpoly(A,Coefvec)
        deg = (sqrt(1+8*(length(Coefvec)+1))-1)/2 -1 ;
        newCoef = zeros(1,length(Coefvec))  ;

       
        function cval = ceval(j,h,k)
            if j > (k-h)
                cval = 0 ;
            else
                cval = nchoosek(k-h,j)*(A(1,1)^(k-h-j))*(A(1,2)^(j)) ;
            end
        end

        function dval = deval(j,h)
            if j > h
                dval = 0 ;
            else
                dval = nchoosek(h,j)*(A(2,1)^(h-j))*(A(2,2)^(j)) ;
            end
        end

        Aco = cell(1,deg+1)   ;
        count=0 ;
        for k = 2:deg+1
        Aco{k} = zeros(k-1) ;
        for h = 0:k-1    
            for j = 0:k-1
                temp = 0;  
                for p = 0:j
                    temp = temp + ceval(p,h,k-1)*deval(j-p,h) ;
                end        
                Aco{k}(h+1,j+1) = temp ;
            end  
        end
        newCoef(1,count+1:count+k) = (Aco{k}')\Coefvec(1,count+1:count+k)' ;
        count = count + k ;
        end

        end
    end

end