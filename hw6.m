clear all
  %% Initial condition
  h = 0.1;
  % Define the domain
  x = 0: h: 2;
  y = 0: h: 1;
  M = length(x);
  N = length(y);
  dx = x(2) - x(1);
  dy = y(2) - y(1); %new

  % Set initial condition for u
u = zeros(M,N);

for i=1:N
    u(1,i) = 1;  % u(x = 0,y) = 0
    u(M,i) = i;  % u(x = 2,y) = y
end

u_Point_Jacobi = u;
u_Point_Jacobi_old = u;
u_Gauss_Seidel = u;
u_Gauss_Seidel_old = u;
u_SOR = u;
u_SOR_old = u;
  
  % Define maximum number of iterations
  MAX_NUMBER_ITERS = 1000;
  
  for iter = 1: MAX_NUMBER_ITERS
    % Using Point Jacobi method to solve for the phi;  
    % For all points inside the loop
    for i = 2: M-1
      for j = 2: N-1
        % New    =  Old value
        u_Point_Jacobi(i,j) = ((1/(dx*dx)) * ( u_Point_Jacobi_old(i-1,j) + u_Point_Jacobi_old(i+1,j)) + (1/(dy*dy)) * ( u_Point_Jacobi_old(i,j-1) + u_Point_Jacobi_old(i,j+1)))/(2*((1/(dx*dx))+1/(dy*dy)));
        u_Gauss_Seidel(i,j) = ((1/(dx*dx)) * ( u_Gauss_Seidel(i-1,j) + u_Gauss_Seidel_old(i+1,j)) + (1/(dy*dy)) * ( u_Gauss_Seidel(i,j-1) + u_Gauss_Seidel_old(i,j+1)))/(2*((1/(dx*dx))+1/(dy*dy)));
      end
    end
       % Assign old value back
    u_Point_Jacobi_old = u_Point_Jacobi;
   
  
  % Plot the value of phi in the x direction going through the center
  % center_location_x = floor(M/2);
  % center_location_y = floor(N/2);
  
  % u_Point_Jacobi_y = u_Point_Jacobi(center_location_x,:);
  % u_Point_Jacobi_x = u_Point_Jacobi(:,center_location_y);

 
  
% Exact solution
  x_one = ones(1,11);
  n = 1;
  u_exact = 0;
  if n<Inf
    u_exact = u_exact-4*(1/((n*pi)*(n*pi)*sinh(2*n*pi)))*sinh(n*pi*x)'.*cosh(n*pi*y);
    n = n + 2;
  end
  u_exact =  (x'/4)*x_one-u_exact;


figure(1)
mesh(x,y,u_Point_Jacobi);
xlabel('x','FontSize',12,...
       'FontWeight','bold','Color','k')
ylabel('y','FontSize',12,...
       'FontWeight','bold','Color','k')
zlabel('u_Point_Jacobi','FontSize',12,...
       'FontWeight','bold','Color','k')

figure(2)
mesh(x,y,u_Gauss_Seidel);
xlabel('x','FontSize',12,...
       'FontWeight','bold','Color','k')
ylabel('y','FontSize',12,...
       'FontWeight','bold','Color','k')
zlabel('u_Gauss_Seidel','FontSize',12,...
       'FontWeight','bold','Color','k')

figure(3)
mesh(x,y,u_exact);
xlabel('x','FontSize',12,...
       'FontWeight','bold','Color','k')
ylabel('y','FontSize',12,...
       'FontWeight','bold','Color','k')
zlabel('u_exact)','FontSize',12,...
       'FontWeight','bold','Color','k')
