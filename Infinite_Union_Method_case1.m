%% Y and Y1 %%
Y=zeros(123,123);
Y(1,31)=1;Y(2,59)=2;Y(3,98)=1;Y(4,87)=2;Y(5,71)=1;Y(6,123)=2;Y(7,78)=5;Y(8,80)=5;Y(9,56)=5;Y(10,49)=2;Y(11,35)=1;Y(12,41)=2;Y(13,43)=2;Y(14,44)=2;Y(15,47)=10;Y(16,53)=2;Y(17,65)=1;
Y(18,68)=5;Y(19,70)=4;Y(20,72)=5;Y(21,91)=2;Y(22,84)=1;Y(23,96)=5;Y(24,100)=2;Y(25,103)=1;Y(26,108)=1;Y(27,112)=5;Y(28,115)=1;Y(29,118)=4;Y(30,122)=0.2;Y(31,34)=5;Y(32,33)=5;
Y(33,34)=2;Y(34,38)=1;Y(35,36)=4;Y(36,37)=5;Y(37,38)=2;Y(38,39)=5;Y(39,40)=1;Y(39,41)=5;Y(41,42)=10;Y(42,43)=5;Y(42,58)=2;Y(42,54)=2;Y(44,45)=1;Y(44,46)=5;Y(46,47)=2;Y(47,48)=2;Y(47,54)=2;Y(48,49)=5;Y(49,50)=5;Y(49,51)=1;Y(51,52)=5;Y(51,53)=1;Y(54,55)=1;Y(54,74)=10;Y(55,56)=4;Y(56,57)=1;Y(58,69)=1;Y(59,60)=1;Y(59,64)=10;Y(60,61)=5;Y(61,62)=4;Y(62,63)=1;Y(63,93)=0.5;Y(64,65)=5;Y(64,66)=2;Y(66,67)=5;Y(66,68)=0.5;Y(68,69)=10;Y(69,70)=2;Y(70,71)=1;Y(72,73)=4;Y(73,77)=2;Y(73,91)=2;Y(74,75)=1;Y(75,76)=5;Y(76,77)=1;Y(77,78)=2;Y(77,79)=1;Y(79,81)=5;Y(80,81)=4;Y(81,82)=4;Y(82,83)=2;Y(82,84)=1;Y(84,85)=4;Y(85,86)=5;Y(85,119)=2;Y(87,88)=1;Y(88,89)=2;
Y(89,90)=5;Y(90,91)=1;Y(91,92)=1;Y(92,109)=1;Y(93,94)=1;Y(94,95)=4;Y(94,99)=1;Y(95,96)=2;Y(96,97)=1;Y(97,98)=5;Y(99,100)=4;Y(99,101)=5;Y(101,102)=4;Y(101,104)=5;Y(102,103)=5;Y(104,105)=2;Y(105,106)=2;Y(105,109)=1;Y(106,107)=5;Y(107,108)=2;Y(109,110)=2;Y(109,113)=1;Y(110,111)=4;Y(111,112)=5;Y(113,114)=4;Y(113,116)=1;Y(114,115)=4;Y(116,117)=4;Y(116,119)=1;Y(117,118)=1;Y(118,120)=2;Y(120,121)=4;Y(121,122)=1;Y(121,123)=1;
Y=diag((Y+Y')*ones(123,1))-Y-Y'; %% Y is obtained %%

YSS=Y(1:10,1:10); YSL=Y(1:10, 11:123); YLS=YSL';YLL=Y(11:123,11:123);
K=diag([0.5 0.5 0.5 0.5 0.5 1 1 1 1 1]);
c=ones(1,10)*K^(-1)*YSS^(-1)*K^(-1)*ones(10,1);alpha=(ones(1,10)*K^(-1)*YSS^(-1)*YSL)';
Y1=YLL-YLS*YSS^(-1)*YSL+alpha*alpha'*c^(-1); %% Y1 is obtained %%
uref=3000;

%% The feasible regions proposed in this paper under %%

fplot(@(x)0.2+x/1.2,[0, 1.2/2.2]); hold on
fplot(@(x)x/0.8-0.2,[0.16, 0.8/1.8]);
fplot(@(x)x.^2./(1-x),[0.8/1.8, 1.2/2.2]);  %% the feasible proposed in this paper %%

%% T1, T2, s and r 
lambda=1.2; %% load factor (Fig. 5 can be obtained by $\lambda=1.2$ and $\lambda=1.4$)
pL=[3*ones(31,1);6*ones(31,1);2*ones(31,1)]*10000*lambda;
pM=[1; 1; 1; 1; 1; 1; 1; 2; 2; 2; 2; 2; 2; 3; 3; 3; 3; 3; 3; 3;]*20000;
T=Y1^(-1)*diag([-pM;pL])/uref/uref;  %% T is obtained 
T1=T.*(T>0); %% T1 is obtained 
T2=-T.*(T<0); %% T2 is obtained  
s=T1*ones(113,1);r=T2*ones(113,1); %% s and r are obtained 

%% Solvability condition: (s_i, r_i) should belongs to the proposed feasible region
plot(s,r,'*')


%% power low calculation (Newton-Raphson method)
x=ones(113,1);
for i=1:5
    x=x-(eye(113)-T*diag(x)^(-2))^(-1)*(x+T*diag(x)^(-1)*ones(113,1)-ones(113,1)); 
end
error=x+T*diag(x)^(-1)*ones(113,1)-ones(113,1);
max(abs(error))  %% if the error is less than 1e-10, the equation is solvable 
max(x)  %% the highest voltage among CPD's voltage
min(x)  %% the lowest voltage among CPD's voltage
u_L=uref*x; %% CPD's voltage is obtained 

i_S=ones(1,113)*([-pM;pL]./u_L)*(ones(1,10)*K^(-1)*ones(10,1))^(-1)*K^(-1)*ones(10,1); %% Currents of DCC-DGs
i_max=[200*ones(5,1);100*ones(5,1)];

tau=ones(1,113)*([-pM;pL]./[1.2*ones(20,1);0.8*ones(93,1)])/(uref*ones(1,10)*i_max) %% if tau is less than 1, then i_S <i_max. 
