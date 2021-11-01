function [ A ] = generateArandom( m,n )
%% Generates A = {v1, ..., vn} randomly according to a uniform distribution on the unit ball of Rm
%
% Syntax: 
%       [ A ] = generateArandom( m,n )
%
% Input: 
%       m: Dimension 
%       n: Amount of points
%
% Outside:
%       A: Matrix A randomly according to a uniform distribution on the unit ball of Rm
%
    j=1;
    A=zeros(m,n);
    while j<=n
        vj=randn(m,1);                    %Generate random vectors(m,1)
        vj=vj/norm(vj);                   %Normalizes these vectors to be on the border of unit ball
        vj=(rand)^(1/m)*vj;               %Multiplies by the m-th root of a scalar between 0 and 1 with uniform distribution
        A(:,j)=vj;                        %This vector is the j-th column of A
        j=j+1;
    end
end

