% Supplementary Code for Rabiee, M., Aslani, B., & Rezaei, J. (2021). A decision support system for detecting and handling biased decision-makers in multi criteria group decision-making problems.
% Expert Systems with Applications, 171, 114597.
% https://doi.org/10.1016/j.eswa.2021.114597

% Implementation By Babak Aslani
% 03/07/2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Initilaization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameter setting
alpha=0.95;
B=2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main Loop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Provide the input normalized score matrix in which each column is related to a DM
%%Option 1
% Reading the matrix from an excel file
%a=xlsread('Data.xlsx','A1:B4');

%%Option 2
% Insert the normalized performance matrix directly
a=[
0.486916938	0.570432583	0.670772345	0.593694246	0.36060752;
0.415863049	0.75375489	0.577712187	0.516114074	0.346950812;
0.545364152	0.30970587	0.49204504	0.602286828	0.215574985;
0.450330663	0.671935065	0.724090987	0.604087627	0.350131365;
0.527723246	0.592345116	0.699563641	0.584261375	0.287776678;
0.348133048	0.834233765	0.541779161	0.623329748	0.174433718;
0.488940901	0.612744821	0.750149603	0.681816313	0.161714054;
0.328910609	0.72740824	0.555213158	0.561438752	0.428360647;
0.488665382	0.713862762	0.713332171	0.497396116	0.307851593;
0.37971174	0.432293511	0.628134545	0.497956208	0.398180392;
0.376480982	0.797080933	0.630331695	0.63249759	0.262438643;
0.427774852	0.539190855	0.762354326	0.689662928	0.349487986;
0.331888083	0.276852884	0.705510511	0.516618333	0.207914837;
0.431120358	0.569413742	0.561758829	0.576548989	0.30156787;
0.336988283	0.772262718	0.724249105	0.503000727	0.460469966;
0.545315101	0.426020552	0.642560371	0.590175479	0.249987491;
0.527547882	0.776185946	0.504685218	0.635343066	0.408269273;
0.403658157	0.856619371	0.415158759	0.583495501	0.466041013
];

NO=size(a,1);       % Number of Options*Criteria (required score for each DM)
NDM=size(a,2);      % Number of DMs

% Selecting the desired version of method
choice = menu('Choose the version','EABM','MABM','SABM');
   
%%%%%%%%%%%%%%%%%%%%%%%%%%% Algorithm 1- EABM Version %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% First phase
% Column elimination of biased DMs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% calculate the CI of DMs
for i=1:NDM
SEM(i) = std(a(:,i))/(sqrt(length(a(:,i)))-1);          % Standard Error
ts(i) = tinv(alpha,length(a(:,i))-1);                   % T-Score
CILow(i) = mean(a(:,i)) - ts(i)*SEM(i);                 % Lower bound of CI
CIupp(i) = mean(a(:,i)) + ts(i)*SEM(i);                 % Upper bound of CI
end

% Find the overlaps among DMs
NOR=ones(NDM);
for i=1:NDM
Inter1(i)= fixed.DataSpecification('single', 'Intervals',{CILow(i),CIupp(i)});    
end
N=ones(NDM, NDM);
for i=1:NDM-1
    for j=i+1:NDM
        if overlaps(Inter1(i).Intervals, Inter1(j).Intervals)==0
        NOR(i,j)=0;  
        NOR(j,i)=0; 
        end
    end
end

NNOR=zeros(NDM,1);
for i=1:NDM
    % Reverse the NOR by substracting the toal number of DMs form the
    % calculated value 
NNOR(i)=NDM-sum(NOR(i,:));    
end

% Eliminating the biased DMs by comparing to the predefined threshold 
c=a;
e=0;
for i=1:NDM    
  if NNOR(i)>B
     c(:,i)=0; 
    e=e+1;
  else
     e=e+0; 
  end  
       
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Second phase
% Assign weight based on the combination of overlap ratio and the relative
% CI to the Total CI
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Remove the columns for biased DMs
columnsWithAllZeros = all(c == 0);
d = c(:, ~columnsWithAllZeros);
NNDM=size(d,2);
if NNDM>1
% Calculate the CI Total
TSEM = std(d(:))/sqrt(length(d(:)));             % Standard Error
tts = tinv(alpha,length(d(:))-1);                % T-Score
TCIl = mean(d(:)) - tts*TSEM;                    % Lower bound of CI
TCIu = mean(d(:)) + tts*TSEM;                    % Upper bound of CI
TCI=TCIu-TCIl;
% Calculate the Indivdual CIs
for i=1:NNDM
SE(i) = std(d(:,i))/(sqrt(length(d(:,i)))-1);           % Standard Error
tss(i) = tinv(alpha,length(d(:,i))-1);                   % T-Score
Cil(i) = mean(d(:,i)) - tss(i)*SE(i);                    % Lower bound of Confidence Intervals
Ciu(i) = mean(d(:,i)) + tss(i)*SE(i);                    % Upper bound of Confidence Intervals
l(i)=Ciu(i)-Cil(i);                                     % Length of the confidence Interval
end 
% Calculate the Overlap ratio
for i=1:NNDM
Inter(i)= fixed.DataSpecification('single', 'Intervals', {Cil(i), Ciu(i)});    
end
N=ones(NNDM, NNDM);
for i=1:NNDM-1
    for j=i+1:NNDM
        if overlaps(Inter(i).Intervals, Inter(j).Intervals)==1
ax=intersect(Inter(i).Intervals, Inter(j).Intervals);
ub=ax.LeftEnd;
UB=quantize(ub, numerictype);
lb=ax.RightEnd;
LB=quantize(lb, numerictype);
g=lb-ub;
f = sfi(g);
% %convert to double
h= double(f);
N(i,j)=h;
        else
         N(i,j)=0;   
        end
    end
end
% Mirror the upper triangle values to the lower one
X=N;
for i=1:NNDM-1
    for j=i+1:NNDM
        X(j,i)=X(i,j); 
    end
end
% Change the main diagonal to zero
N=X-eye([NNDM NNDM]);

% Calculate the total Overlap with other DMs
sr=zeros(1,NNDM);
or=zeros(1,NNDM);
for i=1:size(N,1)
sr(i)=sum(N(i,:));  
or(i)=sr(i)/((NNDM-1)*l(i));
end

% Calculate the fraction of CI(i) to the total CI
g=zeros(1,NNDM);
for i=1:NNDM
    g(i)=l(i)/TCI;
end

% Calculate the Product of OR and Relative CI
w=zeros(1,NNDM);
for i=1:NNDM
    w(i)=g(i)*or(i);
end

% Print the final Results
Sum=sum(w(:));
z=zeros(NNDM,1);
for i=1:NNDM
z(i)=w(i)/Sum;    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Third phase
% Form the final results and cacluate the performance measures
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Final Matrix
columnsWithAllZeros = all(c == 0);
fm = c(:, ~columnsWithAllZeros);
x=zeros(1,NDM);
ww=zeros(1,NDM);
for i=1:NDM
    for j=1:NNDM
        if a(:,i)==fm(:,j)
        x(i)=1;
        ww(i)=z(j);
         break;
         else
        x(i)=0;
        ww(i)=0;
        end
    end

end

% Assign the equally distributed weights as the initial ones to DMs and
% Calculate the absolute difference in weights for all DMs
v=0;
for t=1:NDM
    v=v+abs(ww(t)-(1/NDM));
end

%Calculate the weight deviation from equallly distributed scheme for the
%remaining DMs in the pool
WD=0;
for t=1:NDM
    if x(t)==1
    WD=WD+abs(ww(t)-(1/NNDM));
    end
end
AWD=WD/NNDM;

% Calculate MAD and AWD values
if nnz(N) == 0
   MAD=0;
    AWD=0;
else
MAD=v/NDM;
end

else
   MAD=0; 
   AWD=0;
end

% Print the NEDM
NEDM=e;

% Print the Ratio of eliminated DMs during the process
RNEDM=e/NDM;

%%%%%%%%%%%%%%%%%%%%%%%%%%% Algorithm 2- MABM Version %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Second phase
% Assign weight based on the combination of overlap ratio and the relative
% CI to the Total CI
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if NNDM>1
% Calculate the Product of OR and Relative CI
w2=zeros(1,NNDM);
for i=1:NNDM
    w2(i)=((g(i)*or(i))/2);
end

% Print the final Results
Sum2=sum(w2(:));
z2=zeros(NNDM,1);
for i=1:NNDM
z2(i)=w2(i)/Sum2;    
end
weight2=zeros(1,NNDM);
for i=1:NNDM
    weight2(i)=0.5*z2(i)+(1/(2*NNDM));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Third phase
% Form the final results and cacluate the performance measures
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Final Matrix
columnsWithAllZeros2 = all(c == 0);
fm2 = c(:, ~columnsWithAllZeros2);
x2=zeros(1,NDM);
ww2=zeros(1,NDM);
for i=1:NDM
    for j=1:NNDM
        if a(:,i)==fm2(:,j)
        x2(i)=1;
        ww2(i)=weight2(j);
         break;
         else
        x2(i)=0;
        ww2(i)=0;
        end
    end

end

% Assign the equally distributed weights as the initial ones to DMs and
% Calculate the absolute difference in weights for all DMs
v2=0;
for t=1:NDM
    v2=v2+abs(ww2(t)-(1/NDM));
end

%Calculate the weight deviation from equallly distributed scheme for the
%remaining DMs in the pool
WD2=0;
for t=1:NDM
    if x2(t)==1
    WD2=WD2+abs(ww2(t)-(1/NNDM));
    end
end
AWD2=WD2/NNDM;

% Calculate MAD and AWD values
if nnz(N) == 0
   MAD2=0;
    AWD2=0;
else
MAD2=v2/NDM;
end

else
   MAD2=0; 
   AWD2=0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%% Algorithm 3- SABM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Remove the columns for biased DMs
NNDM=size(a,2);
if NNDM>1
% Calculate the CI Total
TSEM1 = std(a(:))/sqrt(length(a(:)));               % Standard Error
tts1 = tinv(alpha,length(a(:))-1);                    % T-Score
TCIl1 = mean(a(:)) - tts1*TSEM1;                    % Confidence Intervals
TCIu1 = mean(a(:)) + tts1*TSEM1;                    % Confidence Intervals
TCI1=TCIu1-TCIl1;
% Calculate the Indivdual CIs
for i=1:NNDM
SE1(i) = std(a(:,i))/(sqrt(length(a(:,i)))-1);         % Standard Error
tss1(i) = tinv(alpha,length(a(:,i))-1);                % T-Score
Cil1(i) = mean(a(:,i)) - tss1(i)*SE1(i);               % Lower bound of CI
Ciu1(i) = mean(a(:,i)) + tss1(i)*SE1(i);               % Upper bound of CI
l1(i)=Ciu1(i)-Cil1(i);                                 % Length of the CI
end 
% Calculate the Overlap ratio
for i=1:NNDM
Inter1(i)= fixed.DataSpecification('single', 'Intervals', {Cil1(i), Ciu1(i)});    
end
N1=ones(NNDM, NNDM);
for i=1:NNDM-1
    for j=i+1:NNDM
        if overlaps(Inter1(i).Intervals, Inter1(j).Intervals)==1
ax1=intersect(Inter1(i).Intervals, Inter1(j).Intervals);
ub1=ax1.LeftEnd;
UB1=quantize(ub1, numerictype);
lb1=ax1.RightEnd;
LB1=quantize(lb1, numerictype);
g1=lb1-ub1;
f1 = sfi(g1);
% %convert to double
h1= double(f1);
N1(i,j)=h1;
        else
         N1(i,j)=0;   
        end
    end
end
% Mirror the upper triangle values to the lower one
X1=N1;
for i=1:NNDM-1
    for j=i+1:NNDM
        X1(j,i)=X1(i,j); 
    end
end
% Change the main diagonal to zero
N1=X1-eye([NNDM NNDM]);

% Calculate the total Overlap with other DMs
sr1=zeros(1,NNDM);
or1=zeros(1,NNDM);
for i=1:size(N1,1)
sr1(i)=sum(N1(i,:));  
or1(i)=sr1(i)/((NNDM-1)*l1(i));
end

% Calculate the fraction of CI(i) to the total CI
g1=zeros(1,NNDM);
for i=1:NNDM
    g1(i)=l1(i)/TCI1;
end

% Calculate the Product of OR and Relative CI
w3=zeros(1,NNDM);
for i=1:NNDM
    w3(i)=(g1(i)*or1(i));
end

% Print the final Results
Sum3=sum(w3(:));
z3=zeros(NNDM,1);
for i=1:NNDM
z3(i)=w3(i)/Sum3;    
end
weight3=zeros(1,NNDM);
for i=1:NNDM
    weight3(i)=0.5*z3(i)+(1/(2*NDM));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Third phase
% Form the final results and cacluate the performance measures
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Final Matrix
columnsWithAllZeros3 = all(a == 0);
fm3 = a(:, ~columnsWithAllZeros3);
x3=zeros(1,NDM);
ww3=zeros(1,NDM);
for i=1:NDM
    for j=1:NNDM
        if a(:,i)==fm3(:,j)
        x3(i)=1;
        ww3(i)=z3(j);
         break;
         else
        x3(i)=0;
        ww3(i)=0;
        end
    end

end

% Assign the equally distributed weights as the initial ones to DMs and
% Calculate the absolute difference in weights for all DMs
v3=0;
for t=1:NDM
    v3=v3+abs(ww3(t)-(1/NDM));
end

%Calculate the weight deviation from equallly distributed scheme for the
%remaining DMs in the pool
WD3=0;
for t=1:NDM
    if x3(t)==1
    WD3=WD3+abs(ww3(t)-(1/NNDM));
    end
end
AWD3=WD3/NNDM;

% Calculate MAD and AWD values
if nnz(N) == 0
   MAD3=0;
    AWD3=0;
else
MAD3=v/NDM;
end

else
   MAD3=0; 
   AWD3=0;
end

if choice==1
    fprintf('You selected the EABM version\n');  
fprintf('Results for EABM Version \n');

fprintf('********************************\n');
fprintf('*Elimination information*\n');
fprintf('The Number of eliminated DMs is %d \n', NEDM);
fprintf('The Ratio of eliminated DMs is %f \n', RNEDM);
fprintf('********************************\n');


fprintf('*Performance Measures*\n');
fprintf('The MAD value is %f \n', MAD);
fprintf('The AWD value is %f \n', AWD);
fprintf('********************************\n');


fprintf('*Final Weights*\n');
wwt=transpose(ww);
B = [(1:size(wwt,1))' wwt];
VarNames = {'DM', 'Weight'};
fprintf(1, '  %s\t\t%s\t\t\t\n', VarNames{:})
fprintf(1, '\t%d\t%06f\n', B')
fprintf('********************************\n');
end

if choice==2
fprintf('You selected the MABM version\n');  
fprintf('Results for MABM Version \n');
fprintf('********************************\n');
fprintf('*Elimination information*\n');
fprintf('The Number of eliminated DMs is %d \n', NEDM);
fprintf('The Ratio of eliminated DMs is %f \n', RNEDM);
fprintf('********************************\n');

fprintf('*Performance Measures*\n');
fprintf('The MAD value is %f \n', MAD2);
fprintf('The AWD value is %f \n', AWD2);
fprintf('********************************\n');

fprintf('*Final Weights*\n');
wweight2=transpose(ww2);
B2 = [(1:size(wweight2,1))' wweight2];
VarNames2 = {'DM', 'Weight'};
fprintf(1, '  %s\t\t%s\t\t\t\n', VarNames2{:})
fprintf(1, '\t%d\t%06f\n', B2')
fprintf('********************************\n');
end


if choice==3
    fprintf('You selected the SABM version\n');  
fprintf('Results for SABM Version \n');
fprintf('********************************\n');

fprintf('*Performance Measures*\n');
fprintf('The MAD value is %f \n', MAD3);
fprintf('The AWD value is %f \n', AWD3);
fprintf('********************************\n');

fprintf('*Final Weights*\n');
wweight3=transpose(weight3);
B3 = [(1:size(wweight3,1))' wweight3];
VarNames3 = {'DM', 'Weight'};
fprintf(1, '  %s\t\t%s\t\t\t\n', VarNames3{:})
fprintf(1, '\t%d\t%06f\n', B3')
fprintf('********************************\n');
end
