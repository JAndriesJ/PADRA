classdef Zircon< shape.C2boundary
    % Class for regular n-polygon shape 
    properties(SetAccess = protected)
        a= 1; 
		phi = 0; % orientation 
        degree = 8 ;
    end
    
    methods
      function obj = Zircon(a,nbPoints,dspl)  
            if nargin < 3
               dspl = 1 ; 
            end            
            k = 4 ;
            b = a/4 ;
            selector1 = pi*(0:2/(5*nbPoints):(1-2/(5*nbPoints))) ;
            selector2 = 0:1/(5*nbPoints):(1-1/(5*nbPoints)) ;
          
			edge1  = [b*(k-1)*cos(selector1)+b*cos((k-1)*selector1) ; b*(k-1)*sin(selector1)-b*sin((k-1)*selector1) ];          
            edge2  = a*(([-1/2;-sqrt(3)/2] - [-1;0])*selector2 + [-1;0])  ;
            edge3  = a*(([1/2;-sqrt(3)/2] - [-1/2;-sqrt(3)/2])*selector2 + [-1/2;-sqrt(3)/2]) ;
            edge4  = a*(([1;0] - [1/2;-sqrt(3)/2])*selector2 + [1/2;-sqrt(3)/2]) ;
            
            
            points = [edge1, edge2, edge3, edge4] ;            
            points = circshift(points,[0,floor(nbPoints/pi*2)]); % This is so that the begining and end are centered on an edge        
            points = points - points(:,500) ;
            
            theta = 0:(1/length(points)):(1-(1/length(points))) ;            
            [points, tvec, avec, normal] = shape.C2boundary.rescale(points, theta, nbPoints, [], dspl);

            obj = obj@shape.C2boundary(points, tvec, avec, normal, [0;0], 'Zircon');
            
            obj.a = 0;
            obj.phi = 0; 
                       
        end
        function obj = mtimes(obj, s)
			% Overload of the operator *
			obj = mtimes@shape.C2boundary(obj, s);
			obj.radius_r1 = obj.radius_r1*s ;  
            obj.radius_r2 = obj.radius_r2*s ; 
            obj.seperation_a = obj.seperation_a*s ;
		end
		
		function obj = lt(obj, phi)
			% Redefine the < operator as the rotation of the boundary
			obj = lt@shape.C2boundary(obj,phi);
			obj.phi = obj.phi+phi;
		end
            
    end
end      