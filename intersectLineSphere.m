function point = intersectLineSphere(line, sphere)
%INTERSECTLINESPHERE Return intersection points between a line and a sphere
%
%   GC = intersectLineSphere(LINE, SPHERE);
%   Returns the two points which are the intersection of the given line and
%   sphere. 
%   LINE   : [x0 y0 z0  dx dy dz]
%   SPHERE : [xc yc zc  R]
%   GC     : [x1 y1 z1 ; x2 y2 z2]
%   
%   See also
%   spheres, circles3d, intersectPlaneSphere
%
%   ---------
%   author : David Legland 
%   INRA - TPV URPOI - BIA IMASTE
%   created the 18/02/2005.
%

%   HISTORY
%   2011-06-21 bug for tangent lines, add tolerance

% check if user-defined tolerance is given

% difference between centers

% equation coefficients
a = line(4)^2 + line(5)^2 + line(6)^2;
b = 2*(line(4)*line(1) + line(5)*line(2) + line(6)*line(3));
c = line(1)^2 + line(2)^2 + line(3)^2 - sphere(:,4).*sphere(:,4);

% solve equation
delta = b.*b - 4*a.*c;

if delta > 0
    % delta positive: find two roots of second order equation
    u1 = (-b -sqrt(delta)) / 2 / a;
    u2 = (-b +sqrt(delta)) / 2 / a;
    
    % convert into 3D coordinate
    point = [line(1:3)+u1*line(4:6) ; line(1:3)+u2*line(4:6)];

elseif delta==0
    % delta around zero: find unique root, and convert to 3D coord.
    u = -b/2./a;    
    point = line(1:3) + u*line(4:6);
    
else
    % delta negative: no solution
    point = ones(2, 3);
    point(:) = NaN;
end