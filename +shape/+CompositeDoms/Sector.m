classdef Sector < shape.C2boundary
 
    properties(SetAccess = protected)
		scale = 1 ;
        angle = pi;
        phi = 0;
        degree = 4;
    end
    
    methods
      function obj =Sector(a,angle,nbPoints,dspl)
          
            if nargin < 4
               dspl = 1; 
            end
            
            selector1  = floor(nbPoints*(angle*a/(angle*a+2*a)) )  ;
            selector2  = floor(nbPoints*(a/(angle*a+2*a)) )  ;

    
            edge1 = [linspace(0,a,selector2);zeros(1,selector2)] ;
            edge2 = a*[cos(linspace(0,angle,selector1));sin(linspace(0,angle,selector1))] ;
            edge3 = [edge2(1,end)*fliplr(linspace(0,1,selector2)); edge2(2,end)*fliplr(linspace(0,1,selector2)) ]   ;
            edge2 = edge2(:,2:end-1) ;
 
%             vedge1 = [a;0] ;
%             vedge2 = a*[-sin(linspace(0,angle,selector1));cos(linspace(0,angle,selector1))] ;
%             vedge3 = [-edge2(1,end); -edge2(2,end)]   ;
%             vedge2 = vedge2(:,2:end-1) ;
            
            
%             nedge1 =  [0,-1;1,0]*vedge1 ;
%             nedge2 =  [0,-1;1,0]*vedge2 ;
%             nedge3 =  [0,-1;1,0]*vedge3 ;
% 
%             
%             aedge1 = repmat([0;0],1,selector2);
%             aedge2 = -a*[cos(linspace(0,angle,selector1));sin(linspace(0,angle,selector1))] ;
%             aedge3 = repmat([0;0],1,floor(nbPoints/3)); 
%             aedge2 = aedge2(:,2:end-1) ;
   
            points = [edge1,edge2,edge3] ;
%             tvec   = [vedge1,vedge2,vedge3]   ;
%             avec   = [aedge1,aedge2,aedge3]   ;
%             normal = [nedge1,nedge2,nedge3]   ;
%             normal = [ normal(1,:)./sqrt(sum((normal.^2),1)) ;  normal(2,:)./sqrt(sum((normal.^2),1))] ;

            points = circshift(points,[0, floor(nbPoints/7*pi)]);
            points = points - points(:,1) ;
%             tvec = circshift(tvec,[0, floor(nbPoints*(angle*a/(angle*a+2*a))/2)]); 
%             avec = circshift(avec,[0, floor(nbPoints*(angle*a/(angle*a+2*a))/2)]); 
%             normal = circshift(normal,[0, floor(nbPoints*(angle*a/(angle*a+2*a))/2)]); 
            
            theta = 0:(1/length(points)):(1-(1/length(points))) ;
            [points, tvec, avec, normal] = shape.C2boundary.rescale(points, theta, nbPoints, [], dspl);

            obj = obj@shape.C2boundary(points, tvec, avec, normal, [0;0], 'Sector');
            
            obj.scale = a;
            obj.angle = angle ; 
            obj.phi = 0 ;

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
