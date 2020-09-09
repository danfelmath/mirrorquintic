 %% The next code computes the orbit of the vector v, by the monodromy matrices M0, M1 of the mirror CY threefold, modulo p.
 %% This new version includes the new vector computed in the actual list, therefore the code is faster than the original code.
clc; clear;
 tic
v=[0 1 0 0]; %Initial vector:  [0,1,0,0] is the 3-torus delta2. [0,0,0,1] is the 3-sphere delta4.
d=1; k=3;    %(d,k) for 14 examples of monodromy for Mirror threefold.  (5,5) is the mirro quintic.
p=23;         %Prime number to consider Zp  

M0=[1 1 0 0;0 1 0 0;d d 1 0;0 -k -1 1]; %Monodromy around 0
M1=[1 0 0 0;0 1 0 1;0 0 1 0;0 0 0 1];   %Monodromy around 1

if p>5  %In order to do more efficient the algorithm: for p>6 the identity mod(M0^(p),p)=Id holds.
    P=p;
else
    P=p^2;
end

Wp=v(1)*p^3+v(2)*p^2+v(3)*p+v(4);  %It is the list Wou of vectors, written as list of number in p-decimal
W0=0;

Cond=0; %% Condition to finish
Cont=1; %% Counting steps
while Cond==0
    fprintf('%d\n',Cont)
    Wd=setdiff(Wp,W0); %Difference between new and old list
    W0=Wp;
    for i=1:length(Wd)
        Wn=orbita0(Wd(i),p,p,P,d,k); %% Function explained below.
        Waux=[Wp Wn];
        Wp=unique(Waux); 
        Cond=(length(Wp)==length(W0)) || length(Wp)==p^4-1; %% Condition to stop: The new orbit is equal to previuos, or Orbp=H1
    end
    Cont=Cont+1;
end

for i=1:length(Wp)
Orbp1(i,:)=[rem(floor(Wp(i)*p.^(-3:0)),p),0,Wp(i)];
Orbp(i,:)=[rem(floor(Wp(i)*p.^(-3:0)),p)];
end
Lorbit=length(Wp); % Size of the orbit


fprintf('Orbit of v=%d\n');
fprintf('Size of the orbit=%d\n',Lorbit);
fprintf('Size diference between the whole H1(X,Zp) and the orbit=%d\n',p^4-Lorbit);
toc
%Orbp


%% Function which compute the orbit of the vector associated to M=v(1)p^3+v(2)p^2+v(3)+v(4)
%% With monodromy Matrix M0 with values d,k. 
%% The powers of M0^k for k=1...L0, and M1^l for l=1..L1
function N=orbita0(M,p,L1,L0,d,k)
%Help
M0=[1 1 0 0;0 1 0 0;d d 1 0;0 -k -1 1]; %Monodromy around 0
M1=[1 0 0 0;0 1 0 1;0 0 1 0;0 0 0 1];   %Monodromy around 1
vi=rem(floor(round(M*p.^(-3:0),5)),p); %Decimal to p-esimal 
for l=1:L1
    for m=1:L0
        w= mod(vi*M1^l*M0^m,p);
        n=w(1)*p^3+w(2)*p^2+w(3)*p+w(4);
        N(n+1)=n;
    end
end
N=nonzeros(N)';
end
