classdef TriloShell < shape.C2boundary
 
    properties(SetAccess = protected)
		scale = 1;
        phi = 0;
        degree = 6 ;
    end
    
    methods
      function obj = TriloShell(a,nbPoints,dspl)
            if nargin < 3
                dspl = 1 ;
            end
            selector1  = floor(nbPoints*(2*pi*a/(a*(2*pi+1))) )  ;            
            selector3  = fliplr(selector1) ;
            
            
            edge1 = repmat(a*exp(linspace(1,2,selector1)-ones(1,selector1)),2,1).*[cos(linspace(0,pi,selector1));sin(linspace(0,pi,selector1))] ;
            edge1 = edge1(:,2:end-1) ;
            edge3 = fliplr(repmat(a/2*exp(linspace(1,2,selector3)-ones(1,selector3)),2,1).*[cos(linspace(0,pi,selector3));sin(linspace(0,pi,selector3))] );
            edge3 = edge3(:,2:end-1) ;
            edge2 = [(linspace(edge1(1,end),edge3(1,1),100));zeros(1,100)] ;
            edge2 = edge2(:,2:end-1) ;
            edge4 = fliplr([(linspace(edge1(1,1),edge3(1,end),100));zeros(1,100)]) ;
            edge4 = edge4(:,2:end-1) ;
            
            points = [1/4,0;0,1/4]*[edge1,edge2,edge3,edge4] ;
            points = circshift(points,[0, round(selector1/pi)]); 
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