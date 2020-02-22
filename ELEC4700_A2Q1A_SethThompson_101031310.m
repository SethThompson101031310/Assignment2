%% ELEC 4700 Assignment 2 Question 1 PART A

% Name: Seth Thompson
% Student Number: 101031310
close all
clear
clc

% While working on this assignment, using matricies in the form GV = F
% would be benificial, as mentioned in the assignment handout.

% Defining the length and width of the plot.
Width = 1;
Length = (3/2)*Width;

% Creating a nx by ny matrix to hold the solution to the system and a B
% vector
nx = 20;
ny = 20;
Va = zeros(nx*ny,nx*ny);
Ba = zeros(nx*ny,1);

% Defining Vo, as mentioned in the assignment handout
Vo = 1;

% Populating the sparse matrix Va

for i = 1 : nx
    for j = 1 : ny
        % Node mapping to put entries into the correct place
        n = j + (i - 1)*ny;
        % Calculating deltas in the x and y direction
        nxm = j + (i - 2)*ny;
        nxp = j + i*ny;
        nym = (j - 1) + (i - 1)*ny;
        nyp = (j + 1) + (i - 1)*ny;
        % Begin inputting all of the correct entires!
        if (i == 1 || i == nx || j == 1 || j == ny)
            Va(n,n) = 1;
        else
            Va(n,n) = -4;
            Va(n,nxm) = 1;
            Va(n,nxp) = 1;
            Va(n,nym) = 1;
            Va(n,nyp) = 1;
        end
        
        % Add a seperate if statement to account for boudaries along the
        % y-axis. This part accounts for the derivatives at the boundaries
        % onthe y-axis
        if i > (nx-(nx-1)) && i < nx
            Va(i,i + nx) = -1;
            Va(i + (ny - 1)*nx,i + (ny - 2)*nx) = -1;
        end
    end
end

% Looking at the sparsity of matrix Va
figure(1)
spy(Va)
title({'Sparsity of matrix Va (partA)','Seth Thompson | 101031310'})

% A seperate loop must be made to populate the B vector with its contents.

for k = 1: ny
    Ba(1 + (k - 1)*ny,1) = Vo;
    % No need to set other boundary equal to zero since it is already
    % initialized to that!
end

% If Ax = B, then x = A\B
x = Va\Ba;
for i = 1:nx
    for j = 1:ny
        % Populating an array with the Z values of the final solution.
        n = j + (i - 1)*ny;
        xArray(i,j) = x(n);
    end
end

% Plotting the output with the correct aspect ratio (3/2)
xDOM = linspace(0,Length,nx);
yDOM = linspace(0,Width,ny);
figure(2)
surf(xDOM,yDOM,xArray)
title({'Solution Obtained from part A','Seth Thompson | 101031310'})
xlabel('x-Axis')
ylabel('y-Axis')
zlabel('z-Axis')

% The solution obtained from this code give a flat plane at an angle. After
% some inspection , it looks as if the simulation has worked correctly
% since the final results starts at the value set for Vo and ends at 0,
% just like how the boundary conditions described!