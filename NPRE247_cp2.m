function [M,F,eig_value,eig_vector,eig_max_value,eig_max_vector] = NPRE247_cp2(group)
% M: Migration matrix
% F: Fission matrix
% eig_value: eigenvalues
% eig_vector: eigenvectors, column vector corresponding to the eigenvalues
% eig_max_value: largest eigenvalue
% eig_max_vector: eigenvector corresponding to the largest eigenvalue
%
% by Hao Xiong, Nov. 17, 2013

if nargin < 1
  group = 2;
end

% Read data and init matrix -----------------------------------------------
file_input = ['cross_section_',int2str(group),'g.txt'];
data = load(file_input,'-ascii');
N = size(data,1);
A    = zeros(N);
Sout = zeros(N);
Sin  = data(:,5:end);
F    = data(:,4)*data(:,3)';

for i = 1:N
  A(i,i) = data(i,2);
  Sout(i,i) = sum(data(:,i+4))-data(i,i+4);
  Sin(i,i) = 0;
end

% Numerical calculation of eigenvalues/eigenvectors -----------------------
M = A+Sout-Sin;
B = inv(M)*F;
[V,D] = eig(B);
eig_value = 1:N;
eig_vector = zeros(N);
NN = 0;
for i = 1:N
  if isreal(D(i,i))
    NN = NN+1;
    eig_value(NN) = D(i,i);
    eig_vector(:,NN) = V(:,i);
  end
end
eig_value(NN+1:end) = [];
eig_vector(:,NN+1:end) = [];

% Power iteration ---------------------------------------------------------
MAXIT = 2;
k = 1;
phi = ones(N,1);
for i = 1:MAXIT
  phi = B*phi/norm(B*phi);
  k = (B*phi)'*phi/(phi'*phi);
end
eig_max_value = k;
eig_max_vector = phi;

% Output ------------------------------------------------------------------
file_output = ['output_',int2str(group),'g.txt'];
fn = fopen(file_output,'w');
fprintf(fn,'%% Input data\r\n');
fprintf(fn,[repmat('%f ',1,size(data,2)),'\r\n'],data');
fprintf(fn,'\r\n%% Migration matrix\r\n');
fprintf(fn,[repmat('%f ',1,N),'\r\n'],M');
fprintf(fn,'\r\n%% Fission matrix\r\n');
fprintf(fn,[repmat('%f ',1,N),'\r\n'],F');
fprintf(fn,'\r\n%% Eigenvalues\r\n');
fprintf(fn,[repmat('%f ',1,NN),'\r\n'],eig_value');
fprintf(fn,'\r\n%% Eigenvectors\r\n');
fprintf(fn,[repmat('%f ',1,NN),'\r\n'],eig_vector');
fprintf(fn,'\r\n%% Largest eigenvalue\r\n');
fprintf(fn,'%f\r\n',eig_max_value);
fprintf(fn,'\r\n%% Eigenvector corresponding to the largest eigenvalue\r\n');
fprintf(fn,'%f\r\n',eig_max_vector);
fclose(fn);
end