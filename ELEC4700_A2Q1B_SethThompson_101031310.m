%% ELEC 4700 Assignment 2 Question 1 PART B

% Name: Seth Thompson
% Student Number: 101031310
close all
clear
clc

% To solve this part of the assignment, the GV = F method will be used
% again.

% Defining the length and width of the plot.
Width = 1;
Length = (3/2)*Width;

% Creating a nx by ny matrix to hold the solution to the system and a B
% vector
nx = 20;
ny = 20;
Vb = zeros(nx*ny,nx*ny);
Bb = zeros(nx*ny,1);

% Defining Vo, as mentioned in the assignment handout
Vo = 1;

for i = 1:nx
    for j = 1:ny
        % Node mapping to put entries into the correct place
        n = j + (i - 1)*ny;
        % Calculating deltas in the x and y direction
        nxm = j + (i - 2)*ny;
        nxp = j + i*ny;
        nym = (j - 1) + (i - 1)*ny;
        nyp = (j + 1) + (i - 1)*ny;
        % Begin inputting all of the correct entires!
        if (i == 1 || i == nx || j == 1 || j == ny)
            Vb(n,n) = 1;
        else
            Vb(n,n) = -4;
            Vb(n,nxm) = 1;
            Vb(n,nxp) = 1;
            Vb(n,nym) = 1;
            Vb(n,nyp) = 1;
        end
        % No need for extra bounddary statement condition this time, it is
        % accounted for in the first part of the if-else statement!
    end
end

% Looking at the sparsity of matrix Va
figure(1)
spy(Vb)
title({'Sparsity of matrix Vb (partB)','Seth Thompson | 101031310'})

% A seperate loop must be made to populate the B vector with its contents.

for k = 2: (ny-1)
    Bb(1 + (k - 1)*ny,1) = Vo;
    % This time, the opposite end of the boundary must be set to Vo as well
    Bb(nx + (k - 1)*ny,1) = Vo;
end

% If Ax = B, then x = A\B
x = Vb\Bb;
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
title({'Solution Obtained from part B','Seth Thompson | 101031310'})
xlabel('x-Axis')
ylabel('y-Axis')
zlabel('z-Axis')

% This next part of the script will calculate the analytical solution of
% the system and compare it to the results found previously.

% To start, a matrix will be defined to hold the analytical solution along
% with a maximum number of iterations.

AnalyticSOL = zeros(nx,ny);
nMAX = 100;

% The solution each iteration of the analytical method will be calculated
% all at once using element-wise operations. The niital x and y values
% Will be put on a grid so element wise operations can be done to complete
% This in one loop as opposed to two.

[x,y] = meshgrid(linspace(-Length/2,Length/2,nx), linspace(0,Width,ny));

% Also defining a and b before entering the loop. a spans from 0 to the
% width of the box, and b spans from -L/2 to L/2
a = Width;
b = Length/2;


for iterNUM = 1:nMAX
     % We can only use odd numbers for n, so use the formula shown below
     t = 2*iterNUM - 1;
     % The code to update the solution matrix is shown next
     AnalyticSOL = AnalyticSOL + (4*Vo/pi) .* (1/t) .* (cosh((t*pi).*x./a) ...
         ./ cosh((t*pi).*b./a)) .* sin((t*pi).*y./a);
     figure(3)
     surf(xDOM,yDOM,AnalyticSOL)
     pause(0.001)
end


% Plotting the results from the analytical solution
figure(3)
surf(xDOM,yDOM,AnalyticSOL)
title({'Solution Obtained from Analytical Method','Seth Thompson | 101031310'})
xlabel('x-Axis')
ylabel('y-Axis')
zlabel('z-Axis')

% Comparing the analytical solution to the numerical one shows some very
% similar results. It can be seen that the analytical solution approches
% that of the final one obtained from the first technique, however both
% methods have their pros and cons

% The numerical solution came to an answer very quickly and appears to be
% much more accurate than the analytical solution after 100 iterations, but
% it required a very large matrix to solve and as a result more data will
% be taken up to solve this. This can be partially resolved though by
% seting the elements of the numerical method to be sparse since most of
% the elements are zero.

% The analytical solution took a more brute force approch to the problem
% and in the end was much easier to follow, however it will take more time
% to solve depending on how accurate you want your final solution to be.
% Also, being able to see the final solution as a movie did help while
% debugging the script for this part of the lab.

% In conclusion, both methods are good ways to solve the system, however
% each has their own advanteges and disiadvantages. The numerical method
% comes up with a more acurate solution faster, but the process is harder
% to follow and will take up more memory while running. The analytical
% solution is easier to follow, yet it'll take longer to reach the final 
% solution depending on how accurate you want it to be.
