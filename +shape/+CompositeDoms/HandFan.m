classdef HandFan< shape.C2boundary
    % Class for regular n-polygon shape 
    properties(SetAccess = protected)
        a= 1; 
		phi = 0; % orientation 
        degree = 4 ;
    end
    
    methods
      function obj = HandFan(a,nbPoints,dspl)  
            if nargin < 3
               dspl = 1 ; 
            end
            b = a/2 ;
            selector1 = 0:1/floor(nbPoints/4):1 ;
            selector2 = fliplr(selector1) ;
            selector3 = -(0:1/floor(nbPoints/2):1) ;

            edge1 = [b*(selector1-1); b*selector1] ;
            edge2 = [b*selector1; b*selector2] ;
            edge3 = [b*cos(pi*selector3); b*sin(pi*selector3)] ;
            
            points = [edge1, edge2, edge3] ;
            
            points = circshift(points,[0,floor(nbPoints/4)]); % This is so that the begining and end are centered on an edge        
            points = points - points(:,1) ;
            
            theta = 0:(1/length(points)):(1-(1/length(points))) ;            
            [points, tvec, avec, normal] = shape.C2boundary.rescale(points, theta, nbPoints, [], dspl);

            obj = obj@shape.C2boundary(points, tvec, avec, normal, [0;0], 'Hand Fan');
            
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