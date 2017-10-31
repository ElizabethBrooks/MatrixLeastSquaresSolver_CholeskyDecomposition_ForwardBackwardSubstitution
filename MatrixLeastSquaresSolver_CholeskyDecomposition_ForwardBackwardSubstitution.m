%{
  Name: Elizabeth Brooks
  File: MatrixLeastSquaresSolver_CholeskyDecomposition_ForwardBackwardSubstitution
  Modified: 31 October 2017
%}
%Main function of script
function MatrixLeastSquaresSolver_CholeskyDecomposition_ForwardBackwardSubstitution
%%Sample matrices to test accuracy in solving least squares problems
g = 10^8+0.3 %Initialize gamma
m = solveLeastSquaresB(g)
%
%Cholesky factorization A = LL' at O((n^3)/3) efficiency
function A = choleskyA(A,n)
for j = 1:n %Columnwise on A
  for i = 1:n %Rowwise on A
    if(i < j) %Above diagonal
      A(i,j) = 0;
    elseif(i == j) %Diagonal
      if(A(i,j) <= 0) %Check if SPD
        Err = 'Error: not SPD matrix.' %Display error message
        return; %Exit matrix computations
      end
      A(i,j) = sqrt(A(i,j)); %Corner
    else %Below diagonal
      A(i,j) = A(i,j)/A(j,j); %l
    end
  end
  for z = j+1:n %Columnwise on submatrix of A
    for i = z:n %Rowwise on submatrix of A
      A(i,z) = A(i,z) - (A(i,j)*A(z,j)); %a = ll'
    end
  end
end
end %End function choleskyA
%
%Solve systems of linear equations with forward and backward substitution
%Forward substitution
function y = forwardSubstitutionA(A,b,n)
y = zeros(n,1); %Matrix to be solved, y
for j = 1:n %Columnwise on L
  y(j) = b(j)/A(j,j); %Matrix is SPD, no need to check A(j,j) > 0
  for i = j+1:n %Rowwise on L
    if(A(i,j) ~= 0) %Skip computations with zero multiplication
      b(i) = b(i) - (A(i,j)*y(j));
    end
  end
end
end %End function forwardSubstitutionA
%
%Backward substitution
function x = backwardSubstitutionA(A,y,n)
x = zeros(n,1); %Matrix to be solved, x
for j = n:-1:1 %Decrement columnwise index on L'
  x(j) = y(j)/A(j,j); %Matrix is SPD, no need to check A(j,j) > 0
  for i = 1:j-1 %Rowwise on L'
    if(A(j,i) ~= 0) %Skip computations with zero multiplication
      y(i) = y(i) - (A(j,i)*x(j));
    end
  end
end
end %End function backwardSubstitutionA
%
%Compute the solution to the least squares problem which satisfies B'Bx = B'b
function m = solveLeastSquaresB(g)
  B = [0 2 0;g g 0;g 0 g;0 1 1]; %Initialize 4x3 matrix B
  b = [2;2*g;2*g;2]; %Initialize 4 row vector
  x = ones(3,1); %Set exact solution vector
  E = zeros(3,3); %Initialize computed intermediary solution matrix
  e = zeros(3,1); %Initialize solution vector with zeros
  %Compute B'*B
  for j = 1:3
    for i = 1:3
      for z = 1:4
        E(i,j) = E(i,j) + (B(z,i)*B(z,j));
      end
    end
  end
  %Compute B'*b
  for i = 1:3
    for z = 1:4
      e(i) = e(i) + (B(z,i)*b(z));
    end
  end
  %Compute B'*B*x = B'*b
  E = choleskyA(E,3);
  m = zeros(3,1);
  for j = 1:3
    for i = 1:3
      y = forwardSubstitutionA(E,e,3);
      m = backwardSubstitutionA(E,y,3);
    end
  end
  %B\b
end %End function solveLeastSquaresB
%
end %End main function of script