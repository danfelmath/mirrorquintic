clc
clear
 %% The next code computes the orbit of the vector v, by the monodromy matrices M0, M1 of the mirror CY threefold, modulo p.

v=[0 1 0 0]; %Initial vector:  [0,1,0,0] is the 3-torus delta2. [0,0,0,1] is the 3-sphere delta4.
%v(2,:)=[0 1 0 0];
d=5; k=5;    %(d,k) for 14 examples of monodromy for Mirror threefold.  (5,5) is the mirro quintic.
p=5;         %Prime number to consider Zp  

M0=[1 1 0 0;0 1 0 0;d d 1 0;0 -k -1 1]; %Monodromy around 0
%M0=M0^5;
%M0=[51 90 -25 0;0 1 0 0;100 175 -49 0;-75 -125 35 1]; % Candelas

M1=[1 0 0 0;0 1 0 1;0 0 1 0;0 0 0 1];   %Monodromy around 1

if p>6  %In order to do more efficient the algorithm: for p>6 the identity mod(M0^(5p),p)=Id holds.
    P=p;
else
    P=p^2;
end

Orbp=mod(v,p);
Wp=pscale(Orbp,p);  %It is the list Wou of vectors, written as list of number in p-decimal

norm=1;
while norm>0
    W=Orbp; L=size(Orbp,1);
    l=1;   c=1;
    while l<L+1
        j=0;
        while j<p+1
                i=0;
            while i<P+1
                %vaux=mod(M1^j*(M0^5)^i*W(:,l),p);
                vaux=mod(W(l,:)*M1^j*(M0)^i,p);
                vp=pscale(vaux,p); %It is the list vaux of vectors, written as list of number in p-decimal
                if length(find(Wp==vp))==0
                    Orbp(L+c,:)=vaux;
                    Wp(L+c)=vp;
                    c=c+1; i=i+1;
                else
                    i=i+1;
                end
            end
                j=j+1;
        end
            l=l+1;
    end
        norm=size(Orbp,1)-size(W,1);
end
Orbp=sortrows([Orbp,Wp'],5); %% Orbit of v. In the last column there are the vectors in p-decimal form.
Lorbit=length(Wp); % Size of the orbit

fprintf('Orbit of v=%d\n');
Orbp(:,1:end-1)
fprintf('Size of the orbit=%d\n',Lorbit);
fprintf('Size diference between the whole H1(X,Zp) and the orbit=%d\n',p^4-Lorbit);






%% The next function write any vector [a b c d] in a list as a number in p-decimal, i.e. a*p^3+b*p^2+c*p+d.
function wn=pscale(vect,p)
    for j=1:size(vect,1)
        wn(j)=vect(j,1)*p^3+vect(j,2)*p^2+vect(j,3)*p^1+vect(j,4)*p^0;
    end
end
