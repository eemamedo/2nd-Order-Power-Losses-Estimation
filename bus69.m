format short

linedata = xlsread('linedata.xlsx', 1, 'A3:F70');
busdata = xlsread('busdata.xlsx', 1, 'A3:J71');

Zb = (12.66E3)^2/100E6;
fb = linedata(:,1); %from
tb = linedata(:,2); %to
R =  linedata(:,3)/Zb; %R
X =  linedata(:,4)/Zb; %X
B =  linedata(:,5);
a =  linedata(:,6);

z =  R+1i*X;
y = 1./z;
b= 1i*B;
nbus = max(max(fb),max(tb)); 
nbranch = length(fb); 

Y = zeros(nbus,nbus);
%formation of off-diagonal element
for k=1:nbranch
    Y(fb(k),tb(k))= Y(fb(k),tb(k))-y(k);
    Y(tb(k),fb(k))=Y(fb(k),tb(k));
end
%formation of diagonal element
for m=1:nbus
    for n = 1:nbranch
        if fb(n)== m
            Y(m,m)= Y(m,m)+y(n);
        elseif tb(n)== m
            Y(m,m)= Y(m,m)+y(n);
        end
    end 
end
    
Y;