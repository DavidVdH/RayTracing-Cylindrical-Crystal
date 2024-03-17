effectStrength = [1e-7];
sizes = [0.9, 2.41, 5];

for effect = effectStrength
for size = sizes

clearvars -except sizes size effectStrength effect

%% Simulation parameters

n1 = 1;     % refractive index surrounding medium
n2 = 1.6;   % refractive index crystals
 
global R L
R = 2.5;    % radius of crystal in mm
L = size;     % length of crystal in mm

global chanceScattering % determines the degree of scattering 
chanceScattering = 1e-5;

NT = 1e7*L;   % number of traps
ET = 0.64;%0.5421;   % trap depth eV
sET = 1.15e7; % freq factor detrapping
s = effect;  % extinction coefficient
opticalDepth = s*NT/L/(3.14*R*R); % extinction coefficient * c

%% Initialization of variables and figures

angleCrit = asin(n1/n2);
nRefl = 0;

k = 8.617e-5;
T = 300;

N_detrappedTherm = [];
n_emitted = [];
N_detrappedOpt = [];

%% Simulation

while NT > 5000
    
    N = floor(sET*NT*exp(-ET/(k*T)));
    n = 0;
    n_emitted(end+1) = 0;
    N_detrappedOpt(end+1) = 0;
    N_detrappedTherm(end+1) = N;

    while n < N && NT > 5000
    
        n = n+1;
        NT = NT-1
        opticalDepth = s * NT / L / (3.14*R*R);

        % Determine path length
        maxDistance = -log(rand(1))/opticalDepth;

        % Create initial point and direction
        p0 = createNewPoint();
        e0 = createNewDirection();

        % Initialize variables
        distance = 0;
        trapped = 1;
        randomDirection = 0;
        nRefl = 0;
        broken = 0;
        nOptDetrap = 0;

        % Start ray tracing
        while trapped
            nRefl = nRefl + 1;
            % Determine location of refection/absorption
            [pNew,randomDirection] = findIntersection(p0, e0, distance, maxDistance);

            % Create new random direction if required
            if randomDirection
                e0 = createNewDirection;
                distance = 0;
                randomDirection  = 0;
                nOptDetrap = nOptDetrap + 1;
                NT = NT-1;
                opticalDepth = s * NT / L / (3.14*R*R);
                maxDistance = -log(rand(1))/opticalDepth; 
            end

            % Determine angle of reflection and reflected direction. 
            [angle,eNew] = calculateReflection(pNew,e0); 

            %plot3([p0(1),pNew(1)],[p0(2),pNew(2)],[p0(3),pNew(3)], 'b-') 
            %plot3(pNew(1),pNew(2),pNew(3), 'go') 

            % Keep track of distance
            distance = distance + norm(pNew - p0);

            % check if photon is still in crystal
            if pNew(1)^2+pNew(2)^2 - R^2 > 0 || pNew(3) < 0 || pNew(3) > L || nRefl > 5000
                NT = NT + 1 + nOptDetrap;
                n = n - 1;
                broken = 1;
                break
            end

            % check if internal reflection or not
            if abs(angle) < angleCrit 
                trapped = 0;
                n_emitted(end) = n_emitted(end)+1;
                N_detrappedOpt(end) = N_detrappedOpt(end) + nOptDetrap;
            elseif pNew(3) < 0.2 || pNew(3) > L-0.2
                if sqrt(pNew(1)^2 + pNew(2)^2) > 2.3 
                    trapped = 0;
                    n_emitted(end) = n_emitted(end)+1;
                    N_detrappedOpt(end) = N_detrappedOpt(end) + nOptDetrap;
                else
                    p0 = pNew;
                    e0 = eNew;
                end
            else
               p0 = pNew;
               e0 = eNew;
            end
        end   

        if broken == 0
            escapePositions(n,:) = pNew;
        else 
            broken = 0;
        end

    end
    
end


dlmwrite(strcat('Output\',num2str(effect),'_', num2str(size),'.txt'), [1:length(n_emitted);n_emitted; N_detrappedTherm; N_detrappedOpt]')

end

end

%% Defining functions

function val = sgn(x)
    if x >= 0
        val = +1;
    else
        val = -1;
    end
end

function p0 = createNewPoint()
    global R L

    % random selection of starting point
    x0 = -R + 2*R*rand(1); 
    y0 = -R + 2*R*rand(1);
    z0 = L*rand(1);
       
    while x0^2+y0^2 >= R^2 || z0 == 0 || z0 == L              
        x0 = -R + 2*R*rand(1); 
        y0 = -R + 2*R*rand(1);
        z0 = L*rand(1);
    end    
    
    p0 = [x0,y0,z0];    
    
end

function e0 = createNewDirection()
    % random selection of direction 
    ff = acos(2*rand(1)-1); % angle phi
    tt = 2*pi*rand(1);      % angle theta
    
    ex = sin(ff)*cos(tt); 
    ey = sin(ff)*sin(tt);
    ez = cos(ff);

    e0 = [ex,ey,ez];    
end


function [pNew,randomDirection] = findIntersection(p,e, distance, maxDistance) 
    global R L
    
    % check at which (x,y)-value line would intersect mantle
    
    if e(1) == 0 && e(2) ~= 0 % special case direction parallel to y-axis
        xNew = p(1);
        if e(2) < 0
            yNew = -sqrt(R^2 -xNew^2);
        else
            yNew = sqrt(R^2 -xNew^2);
        end
        zNew = e(3)/e(2)*(yNew - p(2)+(e(2)/e(3)*p(3)));
        
    elseif e(2) == 0 && e(1) ~= 0 % special case direction parallel to x-axis
        yNew = p(2);
        if e(1) < 0
            xNew = -sqrt(R^2 - yNew^2);
        else
            xNew = sqrt(R^2 - yNew^2);
        end
        zNew = e(3)/e(1)*(xNew - p(1)+(e(1)/e(3)*p(3)));
    elseif e(1) == 0 && e(2) == 0 % special case direction parallel to z-axis
        if e(3)>0
            zNew = L;
        else
            zNew = 0;
        end
            xNew = e(1)/e(3)*(zNew - p(3)+(e(3)/e(1)*p(1)));
            yNew = e(2)/e(3)*(zNew - p(3)+(e(3)/e(2)*p(2)));
    else
        [solx, soly] = ...
            solveNonlinearEquations(1, -e(1)/e(2), (e(1)/e(2)*p(2)-p(1)));
    
        % first solution
        xNew = double(solx(1));
        yNew = double(soly(1));
        zNew = e(3)/e(1)*(xNew - p(1)+(e(1)/e(3)*p(3))); % find corresponding z on line
        
        if e(3)*(zNew-p(3)) <= 0 % check direction of solution
            % second solution
            xNew = double(solx(2));
            yNew = double(soly(2));
            zNew = e(3)/e(1)*(xNew - p(1)+(e(1)/e(3)*p(3)));
        end            
 
        if (zNew < L) && (zNew >= 0)
            %nothing
        else
            if zNew < 0
                zNew = 0;
            else
                zNew = L;
            end
            
            xNew = e(1)/e(3)*(zNew - p(3)+(e(3)/e(1)*p(1)));
            yNew = e(2)/e(3)*(zNew - p(3)+(e(3)/e(2)*p(2)));
            
        end
    end
    
    pNew = [xNew, yNew, zNew];
    
    if (distance + norm(p-pNew)) >= maxDistance
        pNew = findAbsorptionPoint(p, e, maxDistance - distance);
        randomDirection = 1; %SHOULD BE 1
    else
        randomDirection = 0;
    end
end

function [angle,eNew] = calculateReflection(p,e) % point of reflection, direction
    
 	global L chanceScattering
    
    % find normal vector at point of Incidence
    
    if p(3) == 0
        n = [0, 0, -1];
    elseif p(3) == L
        n = [0, 0, 1];
    else
        t = atan(p(2)/p(1));
        n = [sgn(p(1))*cos(t),sgn(p(1))*sin(t), 0];
    end
        
    % find angle from dot product
    angle = acos(dot(n,e)/(norm(n)*norm(e)));
    
    % determine new direction
    if p(3) == 0 || p(3) == L
        eNew = e;
        eNew(3) = -eNew(3);
    else
        angleRef = pi/2 + t;
        reflMatrix = [cos(2*angleRef), sin(2*angleRef), 0;...
                      sin(2*angleRef), -cos(2*angleRef), 0;...
                      0, 0, 1];
               
        eNew = e*reflMatrix;
    end  
    
    if rand < chanceScattering %0.001 up to 5/5/2023
        perturbation = [rand*0.01 , 0 , 0;...
                        0 , rand*0.01 , 0;...
                        0 , 0 , rand*0.01];
        eNew = eNew*perturbation;
    end
   
end

function pNew = findAbsorptionPoint(p, e, distance)
    f = @(x)(abs(norm(x*e)-distance));
    sol = fminsearch(f,1); 
    
    pNew = p + sol*e;
end

function [x,y] = solveNonlinearEquations(a, b, c)

    global R
    
    t = a^2+b^2; %temp variable
    
    x1 = -(sqrt(b^2*(R^2*(t)-c^2))+a*c)/(t);
    x2 = (sqrt(b^2*(R^2*(t)-c^2))-a*c)/(t);
    
    y1 = a*(sqrt(b^2*(R^2*(t)-c^2))-b^2*c)/(b*t);
    y2 = -a*(sqrt(b^2*(R^2*(t)-c^2))+b^2*c)/(b*t);
    
    x = [x1,x2];
    y = [y1,y2];
    
end

