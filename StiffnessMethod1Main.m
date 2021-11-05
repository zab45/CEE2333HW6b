%%
% This code is an adaptation of the code from
% http://www.ce.memphis.edu/7117/FEM_codes/Camp's_code/Camp_Matlab_code.html
%
%% Sample FEM code for CEE 2333 Sec 1020
% ergersg
%Clear variables and MATLAB workspace
clear all;
clc;
format shortEng;

%%hey im adding thiss line

%%zach uploading

%% Pre-processing step 
%
% Read inputs from external file
% Input data from Excel spreadsheets 
%
% Input general information

disp("Reading input data");
delimiterIn = ',';
%
%General inputs about the problem
%
filename = 'generalInput.csv';
headerlinesIn = 0;
General = importdata(filename,delimiterIn);
 NumNod = General.data(1,1);   % read the number of nodes.
NumDofPerNode = General.data(2,1);  % read the number of DOFs per node.
NumElem = General.data(3,1);  %read the number of elements.
%
% Compute the total number of DOFs in the problem - we'll need this to read the BCs, so define it here.
%
NumDof = NumDofPerNode*NumNod; 
%
%Node coordinates
%
filename = 'coordinates.csv';
delimiterIn = ',';
headerlinesIn = 2;
geom = importdata(filename,delimiterIn,headerlinesIn)

xx(1:NumNod) = geom.data(1:NumNod,2);

%Element information
filename = 'elements.csv';
delimiterIn = ',';
headerlinesIn = 2;
element = importdata(filename,delimiterIn,headerlinesIn)
C  = element.data(1:NumElem,5:6); % read the connectivity information that defines the nodes of each element. Similar to the comment above, the number of rows being read is equal to NumEle.
t  = element.data(1:NumElem,2); % read the area or thickness (depending on the element type) of each element. Similar to the comment above, the number of rows being read is equal to NumEle.
E  = element.data(1:NumElem,3); % read the Young modulus of each element. Similar to the comment above, the number of rows being read is equal to NumEle.
v  = element.data(1:NumElem,4); % read Poisson's ratio of each element. Similar to the comment above, the number of rows being read is equal to NumEle.
%
%Boundary conditions
filename = 'BC.csv';
delimiterIn = ',';
headerlinesIn = 1;
boundary = importdata(filename,delimiterIn,headerlinesIn)%F  = xlsread(filename,'BC',strcat('C2:C',string(2+NumDof-1))); %From sheet 'BC,' read the vector of prescribed forces.
bcType(1:NumDof) = boundary.data(1:NumDof,2)
bcValue(1:NumDof) = boundary.data(1:NumDof,3)
%
%% Processing step
%
%Define global stiffness matrix by calling the corresponding stiffness
%function based on element type.
%
disp("Forming stiffness matrix");

        K = StiffnessSpring(E,C,xx,t,v,NumElem,NumDof);

disp("Applying boundary conditions");
%
%Eliminate rows and columns to reduce the problem to only free DOFs
%
locFree     = find(bcType==0); %Find DOFs that are free, represented by a value of 0 in boundary. See 'Boundary' sheet in the input file.
locFIX  = find(bcType==1); %Find DOFs that are fixed, corresponding to a value of 1. See comment above.
DofFree = length(locFree); %Number of free DOFs

Kelim = zeros(DofFree); %Pre-allocate Kelim, the stiffness matrix corresponding to free DOFs. The function zeros(n) creates an nxn matrix of zeros.
Felim = zeros(DofFree,1); %Pre-allocate Felim, the force vector corresponding to free DOFs. The function zeros(n,m) creates an nxm matrix of zeros.
%
%Copy values of K corresponding to the free DOFs into Kelim.
%Note that the nested loop below could've been replaced with Kelim = K(loc,loc), which would also be faster.
%
for i = 1:DofFree %Loop through rows
    for j=1:DofFree %Loop through columns
        Kelim(i,j) = K(locFree(i),locFree(j));  
    end
end
%
%Copy values of F corresponding to the free DOFs into Felim.
%Note that the loop below could've been replaced with Felim(:,1) = F(loc),which would also be faster.
%
for i = 1:DofFree
    Felim(i,1) = bcValue(locFree(i));
end
%
% Solve for free DOFs and store results in d. Note the backslash operator
% will find the fastest way to compute the inverse of Kelim. Typically, it
% will be much faster than calling inv(Kelim) directly.
disp("Finding unknown displacements");
d = Kelim\Felim;
%
%Write full vector of displacements
%
displacement = zeros(NumDof,1); %Pre-allocated a zero-length vector for displacements.
displacement(locFree,1) = d; %Write displacements of free DOFs into the displacement vector - the fixed DOFs have zero displacement.
%

%% Post-processing step
% Evaluate element stresses by calling the appropriate stress function
% corresponding to the element type.

        Se =StressSpring(E,C,xx,t,v,displacement,NumElem);


%==============================================================

% Write displacements to Excel file
D = zeros(NumNod,3); %Pre-allocate a matrix D to store the displacements. The number of rows is the number of nodes, and there are three columns, as defined below.
%D(i,:) = [ith node number,ith node x-displacement,ith node y-displacement]
for ip = 1:NumNod %Loop through each node
    node = NumDofPerNode*ip; %'node' variable is basically the y-direction DOF at the node ip, this is a not-so-elegant way of expressing DOFs in terms of number of nodes.
    D(ip,1) = ip; %Column 1: Node number
    D(ip,2) = displacement(node,1); %Column 2: x-displacement
%    D(ip,3) = displacement(node  ,1); %Column 3: y-displacement
end

disp("Writing the displacements in the Excel file");
%
% 
%
header = {'Displacements'};

% output results
fileName = fullfile(pwd, 'ResultsStiffnessMethod.csv');
fid = fopen(fileName, 'wt');
% Write headers
fprintf(fid, 'Displacements\n');
fprintf(fid,'Node \n');

% Write data.
for i=1:NumNod
    fprintf(fid, '%f , %f\n',D(i,1:2));
end

fprintf(fid, 'Stress\n');
fprintf(fid,'Element\n');

% Write data.
for i=1:NumElem
    fprintf(fid,'%i , %f\n',i,Se(i,1));
end

fclose(fid)

 