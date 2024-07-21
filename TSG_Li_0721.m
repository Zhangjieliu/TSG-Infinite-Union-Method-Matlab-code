%% Y and Y1
Y=zeros(123,123);
Y(1,21)=1; Y(2,53)=2; Y(3,91)=1; Y(4,86)=2; Y(5,63)=1; Y(6,123)=2; Y(7,69)=5; Y(8,72)=5; Y(9,47)=5; Y(10,40)=2; Y(11,25)=1; Y(12,54)=1; Y(13,89)=5; Y(14,105)=1; Y(15,121)=0.2; Y(16,45)=2; Y(17,67)=5; Y(18,62)=4; Y(19,35)=2; Y(20,31)=2; Y(21,24)=5; Y(22,23)=5; Y(23,24)=2; Y(24,28)=1; Y(25,26)=4; Y(26,27)=5; Y(27,28)=2; Y(28,29)=5; Y(29,30)=1; Y(29,31)=5; Y(31,32)=10; Y(32,33)=5; Y(32,49)=2; Y(32,60)=2; Y(33,34)=2; Y(35,36)=1; Y(35,37)=5; Y(37,38)=2; Y(38,39)=2; Y(38,42)=10; Y(38,49)=2; Y(39,40)=5; Y(40,41)=5; Y(40,43)=1; Y(43,44)=5; Y(43,45)=1; Y(46,47)=1; Y(47,48)=4; Y(48,49)=1; Y(49,64)=10; Y(50,51)=4; Y(50,87)=1; Y(51,52)=5; Y(52,53)=1; Y(53,55)=10; Y(54,55)=5; Y(55,57)=2; Y(56,57)=5; Y(57,59)=0.5; Y(58,59)=5; Y(59,61)=10; Y(60,61)=1; Y(61,62)=2; Y(62,63)=1; Y(64,65)=1; Y(65,66)=5; Y(66,70)=1; Y(67,68)=4; Y(68,70)=2; Y(68,81)=2; Y(69,70)=2; Y(70,71)=1; Y(71,73)=5; Y(72,73)=4; Y(73,74)=4; Y(74,75)=2; Y(74,76)=1; Y(76,77)=1; Y(76,78)=4; Y(78,79)=5; Y(78,122)=2; Y(80,81)=2; Y(81,82)=1; Y(81,83)=1; Y(82,106)=1; Y(83,84)=5; Y(84,85)=2; Y(85,86)=1; Y(87,88)=0.5; Y(88,93)=1; Y(89,90)=1; Y(89,92)=2; Y(90,91)=5; Y(92,93)=4; Y(93,94)=1; Y(94,95)=4; Y(94,97)=5; Y(95,96)=2; Y(97,98)=4; Y(97,101)=5; Y(98,99)=5; Y(99,100)=1; Y(101,102)=2; Y(102,103)=2; Y(102,106)=1; Y(103,104)=5; Y(104,105)=2; Y(106,107)=2; Y(106,111)=1; Y(107,108)=4; Y(108,109)=5; Y(109,110)=5; Y(111,112)=4; Y(111,115)=1; Y(112,113)=4; Y(113,114)=1; Y(115,116)=4; Y(115,122)=1; Y(116,117)=1; Y(117,118)=4; Y(117,119)=2; Y(119,120)=4; Y(120,121)=1; Y(120,123)=1;
Y=diag((Y+Y')*ones(123,1))-Y-Y'; %% Y is obtained

YSS=Y(1:10,1:10); YSL=Y(1:10, 11:123); YLS=YSL';YLL=Y(11:123,11:123);
K=diag([0.5 0.5 0.5 0.5 0.5 1 1 1 1 1]);
c=ones(1,10)*K^(-1)*YSS^(-1)*K^(-1)*ones(10,1);alpha=(ones(1,10)*K^(-1)*YSS^(-1)*YSL)';
Y1=YLL-YLS*YSS^(-1)*YSL+alpha*alpha'*c^(-1); %% Y1 is obtained
uref=1500;

%% The feasible regions proposed in [31]-[37] and this paper
line([1/4 0],[0 1/4]); hold on;  %% the feasible region proposed in [31]-[33]

f1 = @(x, y) x + y + abs(x - y) + 2 * sqrt((x + y).* abs(x - y)) - 1;
fimplicit(f1, [0 2/3 0 4/3], 'MeshDensity', 1000); %% the feasible region proposed in [34]-[35]

f2=@(x,y)x+y+3*(abs(x-y)*0.5).^(2/3)-1;
fimplicit(f2, [0 2/3 0 4/3], 'MeshDensity', 1000); %% the feasible region proposed in [37] 

line([0 2/3],[1/2 4/3]); fplot(@(x)2*sqrt(x)-1,[1/4 (3-sqrt(5))/2]); 
fplot(@(x)x.^2./(1-x),[(3-sqrt(5))/2 2/3]);  %% the feasible proposed in this paper

%% T1, T2, s and r
lambda=1.5; %% load factor
pL=[4*ones(35,1);5*ones(35,1);6*ones(33,1)]*1000*lambda;
pM=-[1; 1; 1; 2; 2; 2; 3; 3; 3; 3]*20000;
T=Y1^(-1)*diag([pM;pL])/uref/uref;  %% T is obtained
T1=T.*(T>0); %% T1 is obtained
T2=-T.*(T<0); %%T2 is obtained
s=T1*ones(113,1);r=T2*ones(113,1); %% s and r are obtained 

%% Solvability condition: (s_i, r_i) should belongs to the proposed feasible region
plot(s,r,'*')


%% power low calculation (Newton-Raphson method)
x=ones(113,1);
for i=1:5
    x=x-(eye(113)-T*diag(x)^(-2))^(-1)*(x+T*diag(x)^(-1)*ones(113,1)-ones(113,1)); 
end
error=x+T*diag(x)^(-1)*ones(113,1)-ones(113,1);  
u_L=uref*x; %% CPD's voltage is obtained

i_S=ones(1,113)*([pM;pL]./u_L)*(ones(1,10)*K*ones(10,1))^(-1)*K*ones(10,1) %% Currents of DCC-DGs
