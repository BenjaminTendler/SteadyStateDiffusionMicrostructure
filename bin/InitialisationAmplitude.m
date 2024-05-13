function [Signal,State,Gwave,Trans,b] = InitialisationAmplitude(alpha,phi,E1,E2,Scenario)
%%
if strcmp(Scenario,'First')
    %Initialise First Dictionary
    %Obtain initialisation amplitude
    Lrec=1-E1;
    LtoL=E1.*cosd(alpha);
    A=Lrec/(1-LtoL);
    %Calculate rotation matrix
    [rotMat] = EPG_Operator(alpha,phi,E1,E2,0);
    %Initialise parameters
    Signal{1}=transpose(A*rotMat*[0;0;1]);Signal{1}(3)=0;
    State{1}=transpose([1;-1;0]);
    Gwave{1}{1}=1;Gwave{1}{2}=-1;Gwave{1}{3}=[];
    Trans{1}=transpose([1;1;0]);
    b{1}{1}=1;b{1}{2}=1;b{1}{3}=[];
elseif strcmp(Scenario,'Higher')
    %Initialise Second Dictionary
    Signal{1}=[0,1,0];
    State{1}=[0,0,0];
    Gwave{1}{1}=[];Gwave{1}{2}=[];Gwave{1}{3}=[];
    Trans{1}=transpose([0;0;0]);
    b{1}{1}=[];b{1}{2}=[];b{1}{3}=[];
else
    print('Scenario not recognised');
    exit
end