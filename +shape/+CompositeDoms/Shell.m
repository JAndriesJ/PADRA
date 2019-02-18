classdef Shell < shape.C2boundary
 
    properties(SetAccess = protected)
		scale = 1;
        phi = 0;
        degree = 5 ;
    end
    
    methods
      function obj = Shell(a,nbPoints,dspl)
            if nargin < 3
                dspl = 1 ;
            end
            selector1  = floor(nbPoints*(2*pi*a/(a*(2*pi+1))) )  ;
            selector2  = floor(nbPoints*(a/(a*(2*pi+1))) )  ;

            
            edge1 = repmat(a*exp(linspace(1,2,selector1)-ones(1,selector1)),2,1).*[cos(linspace(0,pi*2,selector1));sin(linspace(0,pi*2,selector1))] ;
            edge2 = fliplr([(linspace(edge1(1,1),edge1(1,end),selector2));zeros(1,selector2)]) ;
            edge2 = edge2(:,2:end-1) ;
            
            points = [1/4,0;0,1/4]*[edge1,edge2] ;
            points = circshift(points,[0, round(selector1/2)]); 
            points = points - points(:,1) ;
            
            theta = 0:(1/length(points)):(1-(1/length(points))) ;
            [points, tvec, avec, normal] = shape.C2boundary.rescale(points, theta, nbPoints, [], dspl);            


            obj = obj@shape.C2boundary(points, tvec, avec, normal, [0;0], 'trancendental shell');
            
            obj.scale = a;
            obj.phi = 0;
           
        end
      
        function obj = mtimes(obj, s)
			% Overload of the operator *
			obj = mtimes@shape.C2boundary(obj, s);
			obj.size_a = obj.size_a * s;
		
		end
		
		function obj = lt(obj, phi)
			% Redefine the < operator as the rotation of the boundary
			obj = lt@shape.C2boundary(obj,phi);
			obj.phi = obj.phi+phi;
        end
            
    end
   
end