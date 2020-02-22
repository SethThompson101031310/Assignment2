%% ELEC 4700 Assignment 2 Question 2 Part A

% Name: Seth Thompson
% Student Number: 101031310

close all
clear
clc

% Defining the length and width of the box along with Vo

Length = 2;
Width = 1;
Vo = 1;
% Defining the values for sigma both inside and outside of the boxes.

sigmaIN = 1e-2;
sigmaOUT = 1;

% Defining the dimensions of each of the boxes (wb and lb in figure 3)

WidthBox = 0.4;
LengthBox = 0.4;

% Defining the number of elements that will be in each part of the matrix

nx = 100*Length;
ny = 100*Width;

% Defining the conductivity map

conductivity = zeros(ny,nx);

for k = 1:ny
    for l = 1:nx
        % If the element being accesed in the conductivity matrix
        % IS within one of the boxes, set its value to the lower
        % conductivity value
        if(l >= nx*WidthBox && l <= nx-nx*WidthBox && (k >= ny-ny*LengthBox || k <= ny*LengthBox))
            conductivity(k,l) = sigmaIN;
        % Else, put in the higher value
        else
            conductivity(k,l) = sigmaOUT;
        end
    end
end

% Using the surf function to plot the conductivity map
figure(1)
surf(conductivity)
xlabel('Length')
ylabel('Width')
zlabel('Conductivity')
title({'Conductivity Over the Region','Seth Thompson | 101031310'})
% Pausing so the plot will plot before the rest of the code runs
pause (0.001)
% Creating the G matrix and B vector for the GV = B solution

G = sparse(nx*ny,nx*ny);
B = zeros(nx*ny,1);

% Populating the G matrix
for l = 1:nx
    for k = 1:ny
        
        % Node mapping to put entries into the correct place
        n = k + (l - 1)*ny;
        
        % Calculating deltas in the x and y direction
        nxm = k + (l - 2)*ny;
        nxp = k + l*ny;
        nym = (k - 1) + (l - 1)*ny;
        nyp = (k + 1) + (l - 1)*ny;
        
        % Begin inputting all of the correct entires!
        if(l == 1 || l == nx) % Left or right side of the region, set entries to 1
            G(n,n) = 1;
        elseif (k == 1) % We are along the bottom of the region apply finite difference method as needed
            
            entryYup = (conductivity(k,l)+conductivity(k+1,l))/2;
            entryXup = (conductivity(k,l)+conductivity(k,l+1))/2;
            entryXdown = (conductivity(k,l)+conductivity(k,l-1))/2;
            
            G(n,n) = -(entryYup + entryXup + entryXdown);
            G(n,nyp) = entryYup;
            G(n,nxp) = entryXup;
            G(n,nxm) = entryXdown;
            
        elseif (k == ny) % We are along the top of the region
            
            entryYdown = (conductivity(k,l)+conductivity(k-1,l))/2;
            entryXup = (conductivity(k,l)+conductivity(k,l+1))/2;
            entryXdown = (conductivity(k,l)+conductivity(k,l-1))/2;
            
            G(n,n) = -(entryYdown + entryXup + entryXdown);
            G(n,nym) = entryYdown;
            G(n,nxp) = entryXup;
            G(n,nxm) = entryXdown;
        else % else, apply finite differnce as needed without worrying about going out of bounds...
            
            % Storing elements from conductivity matrix that will be mapped
            % into the G matrix
            entryYup = (conductivity(k,l)+conductivity(k+1,l))/2;
            entryYdown = (conductivity(k,l)+conductivity(k-1,l))/2;
            entryXup = (conductivity(k,l)+conductivity(k,l+1))/2;
            entryXdown = (conductivity(k,l)+conductivity(k,l-1))/2;
            
            G(n,n) = -(entryYup + entryYdown + entryXup + entryXdown);
            G(n,nyp) = entryYup;
            G(n,nym) = entryYdown;
            G(n,nxp) = entryXup;
            G(n,nxm) = entryXdown;
            
        end
    end
end

% Looking at the sparsity of G to see if it has the correct form
figure(2)
spy(G)
% Pausing so the plot will plot before the rest of the code runs
pause (0.001)
% Populating the B vector next...
for l = 1:nx
    for k = 1:ny
        % Node mapping to put entries into the correct place
        n = k + (l - 1)*ny;
        
        % Are we along the left side? if so set the value to Vo
        if (l == 1) 
            B(n) = Vo;
        end
        
        % Anywhere also it should be zero, but that was defined by using
        % the zeros function to make the vector.
    end
end

% Obtaining the solution, V, from V = G\B

V = G\B;

% Moing the Solution V from the signle vector to a matrix 

for l = 1:nx
    for k = 1:ny
        % Node mapping to put entries into the correct place
        n = k + (l - 1)*ny;
        
        MatrixV(k,l) = V(n);
    end
end

% plotting the voltage with surf
figure(3)
surf(MatrixV)
xlabel('Length')
ylabel('Width')
zlabel('Voltage (V)')
title({'Voltage Over the Region','Seth Thompson | 101031310'})
% Pausing so the plot will plot before the rest of the code runs
pause (0.001)
% The simulated voltage plot makes sense, it slowly drops before getting
% near the contacts and the drops rapidly once it passes through them, then
% goes back to dropping slowly to 0.

% Using the Gradient Function to obtain the electric field from the voltage
% (E = -del(V))

[Ex,Ey] = gradient(-MatrixV);

% Using the Quiver Function to plot the electric field over the region
figure(4)
quiver(Ex,Ey)
xlabel('Length')
ylabel('Width')
title({'Electric Field over The Region','Seth Thompson | 101031310'})
ylim([0,100])
xlim([0,200])
% Pausing so the plot will plot before the rest of the code runs
pause (0.001)
% The simulated electric field plot makes sense, it is in the direction of
% higher to lower potential.

% If the units of conductivity are S/m = 1/(ohm)m = A/(V*m), and the electircal 
% Field has units of V/m, then multiplying them together gives (A/(V*m))*(V/m)
% Which ends up being A/(m^2), which is current density!

% Calculating x and y components of the current density
CDx = Ex.*conductivity;
CDy = Ey.*conductivity;

% Plotting the current density
figure(5)
quiver(CDx,CDy)
xlabel('Length')
ylabel('Width')
title({'Current Density Over the Region','Seth Thompson | 101031310'})
ylim([0,100])
xlim([0,200])
% Pausing so the plot will plot before the rest of the code runs
pause (0.001)
% The current density plot makes sense, current is traveling from an area
% of higher potential to lower potential and mostly through the bottle
% neck.

% To calculte the current along the two contact, a definite size would have
% been needed, as the current density would have needed to be multiplied by
% the area of the contact. One dimension of the contact has been defined by
% me, but the other has not (I have the length of the contact, ie the
% length, but not the width, or how long it would travel up along the
% z-axis). As a result of this, the current density along the contacts will
% be found.

% Take the conductivity on the one side of the region, and multiply it by
% the x-component of the electric field (y component is small here and is
% neglegibel. After that take the sum of each entry.
LeftsideCD = sum(conductivity(:,1).*Ex(:,1));
RightsideCD = sum(conductivity(:,nx).*Ex(:,nx));

fprintf('The current along the left side (coming in) is %f*(area of contact) Amps.\n',LeftsideCD)
fprintf('The current along the right side (coming in) is %f*(area of contact) Amps.\n',RightsideCD)

%% Part B) Graphing The Current at the Contacts Versus the Mesh Size

% To plot the current versus the mesh size, The previous analysis must be
% done again but for different mesh sizes, and the results must be stored
% in their own vectors. This is done next

MeshSizeCurrentLEFT = LeftsideCD;
MeshSizeCurrentRIGHT = RightsideCD;

% Creating variables for current versus box lengths/widths

BoxSizeLEFT = LeftsideCD;
BoxSizeRIGHT = RightsideCD;

% Creating a vector for when the conductivity is incremented.
ConducCurLEFT = LeftsideCD;
ConductCurRIGHT = RightsideCD;

% Copy and Paste previous code in this for loop
for SolNumber = 2:5
    % Both the size of nx and ny must be altered by the multiplier:
    % Defining the number of elements that will be in each part of the matrix

    nx = 100*Length*(1+(0.1*SolNumber));
    ny = 100*Width*(1+(0.1*SolNumber));

    % Defining the conductivity map

    conductivity = zeros(ny,nx);

    for k = 1:ny
        for l = 1:nx
            % If the element being accesed in the conductivity matrix
            % IS within one of the boxes, set its value to the lower
            % conductivity value
            if(l >= nx*WidthBox && l <= nx-nx*WidthBox && (k >= ny-ny*LengthBox || k <= ny*LengthBox))
                conductivity(k,l) = sigmaIN;
            % Else, put in the higher value
            else
                conductivity(k,l) = sigmaOUT;
            end
        end
    end
    
    % Creating the G matrix and B vector for the GV = B solution

    G = sparse(nx*ny,nx*ny);
    B = zeros(nx*ny,1);

    % Populating the G matrix
    for l = 1:nx
        for k = 1:ny

            % Node mapping to put entries into the correct place
            n = k + (l - 1)*ny;

            % Calculating deltas in the x and y direction
            nxm = k + (l - 2)*ny;
            nxp = k + l*ny;
            nym = (k - 1) + (l - 1)*ny;
            nyp = (k + 1) + (l - 1)*ny;

            % Begin inputting all of the correct entires!
            if(l == 1 || l == nx) % Left or right side of the region, set entries to 1
                G(n,n) = 1;
            elseif (k == 1) % We are along the bottom of the region apply finite difference method as needed

                entryYup = (conductivity(k,l)+conductivity(k+1,l))/2;
                entryXup = (conductivity(k,l)+conductivity(k,l+1))/2;
                entryXdown = (conductivity(k,l)+conductivity(k,l-1))/2;

                G(n,n) = -(entryYup + entryXup + entryXdown);
                G(n,nyp) = entryYup;
                G(n,nxp) = entryXup;
                G(n,nxm) = entryXdown;

            elseif (k == ny) % We are along the top of the region

                entryYdown = (conductivity(k,l)+conductivity(k-1,l))/2;
                entryXup = (conductivity(k,l)+conductivity(k,l+1))/2;
                entryXdown = (conductivity(k,l)+conductivity(k,l-1))/2;

                G(n,n) = -(entryYdown + entryXup + entryXdown);
                G(n,nym) = entryYdown;
                G(n,nxp) = entryXup;
                G(n,nxm) = entryXdown;
            else % else, apply finite differnce as needed without worrying about going out of bounds...

                % Storing elements from conductivity matrix that will be mapped
                % into the G matrix
                entryYup = (conductivity(k,l)+conductivity(k+1,l))/2;
                entryYdown = (conductivity(k,l)+conductivity(k-1,l))/2;
                entryXup = (conductivity(k,l)+conductivity(k,l+1))/2;
                entryXdown = (conductivity(k,l)+conductivity(k,l-1))/2;

                G(n,n) = -(entryYup + entryYdown + entryXup + entryXdown);
                G(n,nyp) = entryYup;
                G(n,nym) = entryYdown;
                G(n,nxp) = entryXup;
                G(n,nxm) = entryXdown;

            end
        end
    end
    
    % Populating the B vector next...
    for l = 1:nx
        for k = 1:ny
            % Node mapping to put entries into the correct place
            n = k + (l - 1)*ny;

            % Are we along the left side? if so set the value to Vo
            if (l == 1) 
                B(n) = Vo;
            end

            % Anywhere also it should be zero, but that was defined by using
            % the zeros function to make the vector.
        end
    end

    % Obtaining the solution, V, from V = G\B

    V = G\B;

    % Moing the Solution V from the signle vector to a matrix 

    for l = 1:nx
        for k = 1:ny
            % Node mapping to put entries into the correct place
            n = k + (l - 1)*ny;

            MatrixV(k,l) = V(n);
        end
    end
    
    % Using the Gradient Function to obtain the electric field from the voltage
    % (E = -del(V))

    [Ex,Ey] = gradient(-MatrixV);
    
    % Apply same logic as before, yet store results in a vector
    MeshSizeCurrentLEFT(SolNumber) = sum(conductivity(:,1).*Ex(:,1));
    MeshSizeCurrentRIGHT(SolNumber) = sum(conductivity(:,nx).*Ex(:,nx));
end

% With the resulting currents calculated versus the mesh size, a plot can
% be made!

figure(8)
plot(linspace(1,5,5),MeshSizeCurrentLEFT,'r:')
hold on
plot(linspace(1,5,5),MeshSizeCurrentRIGHT,'k--')
hold off
title({'Plot of Mesh Size Versus Current of Left and Right Contacts','Seth Thompson | 101031310'})
xlabel('Increasing Mesh Size')
ylabel({'Currents at Contacts', '(Current Density Multiplied By Area of Contact)'})
legend('Current Through Left Contact','Current Through Right Contact')

% Pausing so the plot will plot before the rest of the code runs
pause (0.001)

% As can be seen in the previous plot, as the mesh size is increased, the
% slope of the plot gets smaller and smaller. This is indicating that it is
% converging to a value, and that value can be considered to be the actual
% current.

%% Part C) Plotting the Current At the Contacts Versus Different Sized Bottle Necks

% To obtain a plot of the current versus various bottle nexk sizes is next.
% The same approach from before must be applied again to obtain a solution.

% Copy and Paste previous code in this for loop
for SolNumber = 2:4
    % Both the size of nx and ny must be altered by the multiplier:
    % Defining the number of elements that will be in each part of the matrix

    nx = 100*Length;
    ny = 100*Width;
    
    % Editing the length and width of the boxes for this part of the
    % assignment
    
    WidthBox = 0.4*(1+(SolNumber/20));
    LengthBox = 0.4*(1+(SolNumber/20));
    
    % Defining the conductivity map

    conductivity = zeros(ny,nx);

    for k = 1:ny
        for l = 1:nx
            % If the element being accesed in the conductivity matrix
            % IS within one of the boxes, set its value to the lower
            % conductivity value
            if(l >= nx*WidthBox && l <= nx-nx*WidthBox && (k >= ny-ny*LengthBox || k <= ny*LengthBox))
                conductivity(k,l) = sigmaIN;
            % Else, put in the higher value
            else
                conductivity(k,l) = sigmaOUT;
            end
        end
    end
    
    % Creating the G matrix and B vector for the GV = B solution

    G = sparse(nx*ny,nx*ny);
    B = zeros(nx*ny,1);

    % Populating the G matrix
    for l = 1:nx
        for k = 1:ny

            % Node mapping to put entries into the correct place
            n = k + (l - 1)*ny;

            % Calculating deltas in the x and y direction
            nxm = k + (l - 2)*ny;
            nxp = k + l*ny;
            nym = (k - 1) + (l - 1)*ny;
            nyp = (k + 1) + (l - 1)*ny;

            % Begin inputting all of the correct entires!
            if(l == 1 || l == nx) % Left or right side of the region, set entries to 1
                G(n,n) = 1;
            elseif (k == 1) % We are along the bottom of the region apply finite difference method as needed

                entryYup = (conductivity(k,l)+conductivity(k+1,l))/2;
                entryXup = (conductivity(k,l)+conductivity(k,l+1))/2;
                entryXdown = (conductivity(k,l)+conductivity(k,l-1))/2;

                G(n,n) = -(entryYup + entryXup + entryXdown);
                G(n,nyp) = entryYup;
                G(n,nxp) = entryXup;
                G(n,nxm) = entryXdown;

            elseif (k == ny) % We are along the top of the region

                entryYdown = (conductivity(k,l)+conductivity(k-1,l))/2;
                entryXup = (conductivity(k,l)+conductivity(k,l+1))/2;
                entryXdown = (conductivity(k,l)+conductivity(k,l-1))/2;

                G(n,n) = -(entryYdown + entryXup + entryXdown);
                G(n,nym) = entryYdown;
                G(n,nxp) = entryXup;
                G(n,nxm) = entryXdown;
            else % else, apply finite differnce as needed without worrying about going out of bounds...

                % Storing elements from conductivity matrix that will be mapped
                % into the G matrix
                entryYup = (conductivity(k,l)+conductivity(k+1,l))/2;
                entryYdown = (conductivity(k,l)+conductivity(k-1,l))/2;
                entryXup = (conductivity(k,l)+conductivity(k,l+1))/2;
                entryXdown = (conductivity(k,l)+conductivity(k,l-1))/2;

                G(n,n) = -(entryYup + entryYdown + entryXup + entryXdown);
                G(n,nyp) = entryYup;
                G(n,nym) = entryYdown;
                G(n,nxp) = entryXup;
                G(n,nxm) = entryXdown;

            end
        end
    end
    
    % Populating the B vector next...
    for l = 1:nx
        for k = 1:ny
            % Node mapping to put entries into the correct place
            n = k + (l - 1)*ny;

            % Are we along the left side? if so set the value to Vo
            if (l == 1) 
                B(n) = Vo;
            end

            % Anywhere also it should be zero, but that was defined by using
            % the zeros function to make the vector.
        end
    end

    % Obtaining the solution, V, from V = G\B

    V = G\B;

    % Moing the Solution V from the signle vector to a matrix 

    % Set MatrixV = 0 to reset its size
    MatrixV = 0;
    for l = 1:nx
        for k = 1:ny
            % Node mapping to put entries into the correct place
            n = k + (l - 1)*ny;

            MatrixV(k,l) = V(n);
        end
    end
    
    % Using the Gradient Function to obtain the electric field from the voltage
    % (E = -del(V))

    [Ex,Ey] = gradient(-MatrixV);
    
    % Apply same logic as before, yet store results in a vector
    BoxSizeLEFT(SolNumber) = sum(conductivity(:,1).*Ex(:,1));
    BoxSizeRIGHT(SolNumber) = sum(conductivity(:,nx).*Ex(:,nx));
end

% With all of the new vectors made, plot the results!

figure(9)
plot(linspace(1,4,4),BoxSizeLEFT,'g:')
hold on
plot(linspace(1,4,4),BoxSizeRIGHT,'b--')
hold off
title({'Plot of Box Size Versus Current of Left and Right Contacts','Seth Thompson | 101031310'})
xlabel('Increasing Box Size')
ylabel({'Currents at Contacts', '(Current Density Multiplied By Area of Contact)'})
legend('Current Through Left Contact','Current Through Right Contact')

% The results from the plot above make sense. Increasing the size of the
% boxes can be though of as increasing the resistance of the device, and if
% the resistance is increased than less current will go through, as shown
% in the plot.

%% Part D) Graph of the Current at the Contacts Versus 

% For the final part of the assignment, the whole soultion must be iterated
% through AGAIN

% Copy and Paste previous code in this for loop
for SolNumber = 2:4
    % Both the size of nx and ny must be altered by the multiplier:
    % Defining the number of elements that will be in each part of the matrix

    nx = 100*Length;
    ny = 100*Width;
    
    % Editing the length and width of the boxes for this part of the
    % assignment
    
    WidthBox = 0.4;
    LengthBox = 0.4;
    
    % Defining the conductivity map

    conductivity = zeros(ny,nx);

    for k = 1:ny
        for l = 1:nx
            % If the element being accesed in the conductivity matrix
            % IS within one of the boxes, set its value to the lower
            % conductivity value
            if(l >= nx*WidthBox && l <= nx-nx*WidthBox && (k >= ny-ny*LengthBox || k <= ny*LengthBox))
                conductivity(k,l) = sigmaIN;
            % Else, put in the higher value
            else
                conductivity(k,l) = sigmaOUT;
            end
        end
    end
    
    % Increasing the size of each element of the conductivity matrix
    conductivity = conductivity*(1+(0.1*SolNumber));
    
    % Creating the G matrix and B vector for the GV = B solution

    G = sparse(nx*ny,nx*ny);
    B = zeros(nx*ny,1);

    % Populating the G matrix
    for l = 1:nx
        for k = 1:ny

            % Node mapping to put entries into the correct place
            n = k + (l - 1)*ny;

            % Calculating deltas in the x and y direction
            nxm = k + (l - 2)*ny;
            nxp = k + l*ny;
            nym = (k - 1) + (l - 1)*ny;
            nyp = (k + 1) + (l - 1)*ny;

            % Begin inputting all of the correct entires!
            if(l == 1 || l == nx) % Left or right side of the region, set entries to 1
                G(n,n) = 1;
            elseif (k == 1) % We are along the bottom of the region apply finite difference method as needed

                entryYup = (conductivity(k,l)+conductivity(k+1,l))/2;
                entryXup = (conductivity(k,l)+conductivity(k,l+1))/2;
                entryXdown = (conductivity(k,l)+conductivity(k,l-1))/2;

                G(n,n) = -(entryYup + entryXup + entryXdown);
                G(n,nyp) = entryYup;
                G(n,nxp) = entryXup;
                G(n,nxm) = entryXdown;

            elseif (k == ny) % We are along the top of the region

                entryYdown = (conductivity(k,l)+conductivity(k-1,l))/2;
                entryXup = (conductivity(k,l)+conductivity(k,l+1))/2;
                entryXdown = (conductivity(k,l)+conductivity(k,l-1))/2;

                G(n,n) = -(entryYdown + entryXup + entryXdown);
                G(n,nym) = entryYdown;
                G(n,nxp) = entryXup;
                G(n,nxm) = entryXdown;
            else % else, apply finite differnce as needed without worrying about going out of bounds...

                % Storing elements from conductivity matrix that will be mapped
                % into the G matrix
                entryYup = (conductivity(k,l)+conductivity(k+1,l))/2;
                entryYdown = (conductivity(k,l)+conductivity(k-1,l))/2;
                entryXup = (conductivity(k,l)+conductivity(k,l+1))/2;
                entryXdown = (conductivity(k,l)+conductivity(k,l-1))/2;

                G(n,n) = -(entryYup + entryYdown + entryXup + entryXdown);
                G(n,nyp) = entryYup;
                G(n,nym) = entryYdown;
                G(n,nxp) = entryXup;
                G(n,nxm) = entryXdown;

            end
        end
    end
    
    % Populating the B vector next...
    for l = 1:nx
        for k = 1:ny
            % Node mapping to put entries into the correct place
            n = k + (l - 1)*ny;

            % Are we along the left side? if so set the value to Vo
            if (l == 1) 
                B(n) = Vo;
            end

            % Anywhere also it should be zero, but that was defined by using
            % the zeros function to make the vector.
        end
    end

    % Obtaining the solution, V, from V = G\B

    V = G\B;

    % Moing the Solution V from the signle vector to a matrix 

    % Set MatrixV = 0 to reset its size
    MatrixV = 0;
    for l = 1:nx
        for k = 1:ny
            % Node mapping to put entries into the correct place
            n = k + (l - 1)*ny;

            MatrixV(k,l) = V(n);
        end
    end
    
    % Using the Gradient Function to obtain the electric field from the voltage
    % (E = -del(V))

    [Ex,Ey] = gradient(-MatrixV);

    % Apply same logic as before, yet store results in a vector
    ConducCurLEFT(SolNumber) = sum(conductivity(:,1).*Ex(:,1));
    ConductCurRIGHT(SolNumber) = sum(conductivity(:,nx).*Ex(:,nx));
end

% Once again, another plot can be made to show these results
figure(10)
plot(linspace(1,4,4),ConducCurLEFT,'k:')
hold on
plot(linspace(1,4,4),ConductCurRIGHT,'c--')
hold off
title({'Plot of Conductivity Size Versus Current of Left and Right Contacts','Seth Thompson | 101031310'})
xlabel('Increasing Conductivity Size')
ylabel({'Currents at Contacts', '(Current Density Multiplied By Area of Contact)'})
legend('Current Through Left Contact','Current Through Right Contact')

% Once again, these results make sense. Conductivity is the inverse of
% resistance, so by increasing the conductivity we are decreasing the
% resistance, and by decreasing the resistance we are allowing more current
% to go through the device, as shown in the previous plot!