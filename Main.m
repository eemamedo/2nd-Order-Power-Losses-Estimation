clc
clear all
format short
beep off
clear all

dgbus = 18; %DG is connected at this bus. Bus is arbitrary
bus69;
Y;
busd = busdata;
BMva = 100E3;
bus = busd(:,1);
V = busd(:,2);
Vsp = busd(:,2);
del = busd(:,3);
Pg = busd(:,4)/BMva;
Qg = busd(:,5)/BMva;
Pl = busd(:,6)/BMva;
Ql = busd(:,7)/BMva;
Qmin = busd(:,8)/BMva;
Qmax = busd(:,9)/BMva;
type = busd(:,10);
P = Pg - Pl;
Q = Qg - Ql;
Psp = P;
Qsp = Q;

G = real(Y);
B = imag(Y);

pv = find(type == 1);
pq = find(type == 3);
npv = length(pv);
npq = length(pq);

Tolerance = 1;
Iteration = 1;
while (Tolerance > 1e-5)

    P = zeros(nbus,1);
    Q = zeros(nbus,1);
    for i = 1:nbus
        for k = 1:nbus
            P(i) = P(i) + V(i)* V(k)*(G(i,k)*cos(del(i)-del(k)) + B(i,k)*sin(del(i)-del(k)));
            Q(i) = Q(i) + V(i)* V(k)*(G(i,k)*sin(del(i)-del(k)) - B(i,k)*cos(del(i)-del(k)));
        end
    end

        for n = 2:nbus
            if type(n) == 2
                QG = Q(n)+Ql(n);
                if QG < Qmin(n)
                    V(n) = V(n) + 0.01;
                elseif QG > Qmax(n)
                    V(n) = V(n) - 0.01;
                end
            end
         end


    dPa = Psp-P;
    dQa = Qsp-Q;
    k = 1;
    dQ = zeros(npq,1);
    for i = 1:nbus
        if type(i) == 3
            dQ(k,1) = dQa(i);
            k = k+1;
        end
    end
    dP = dPa(2:nbus);
    M = [dP; dQ];

    %J1 - Derivative of Real Power Injections with Angles
    J11 = zeros(nbus-1,nbus-1);
    J22 = zeros(nbus-1,npq);
    J21 = zeros(npq,nbus-1);
    J22 = zeros(npq,npq);

    J11 = zeros(nbus-1,nbus-1);
    for i = 1:(nbus-1)
        m = i+1;
        for k = 1:(nbus-1)
            n = k+1;
            if n == m
                for n = 1:nbus
                    J11(i,k) = J11(i,k) + V(m)* V(n)*(-G(m,n)*sin(del(m)-del(n)) + B(m,n)*cos(del(m)-del(n)));
                end
                J11(i,k) = J11(i,k) - V(m)^2*B(m,m);
            else
                J11(i,k) = V(m)* V(n)*(G(m,n)*sin(del(m)-del(n)) - B(m,n)*cos(del(m)-del(n)));
            end
        end
    end

    J12 = zeros(nbus-1,npq);
    for i = 1:(nbus-1)
        m = i+1;
        for k = 1:npq
            n = pq(k);
            if n == m
                for n = 1:nbus
                    J12(i,k) = J12(i,k) + V(n)*(G(m,n)*cos(del(m)-del(n)) + B(m,n)*sin(del(m)-del(n)));
                end
                J12(i,k) = J12(i,k) + V(m)*G(m,m);
            else
                J12(i,k) = V(m)*(G(m,n)*cos(del(m)-del(n)) + B(m,n)*sin(del(m)-del(n)));
            end
        end
    end

    J21 = zeros(npq,nbus-1);
    for i = 1:npq
        m = pq(i);
        for k = 1:(nbus-1)
            n = k+1;
            if n == m
                for n = 1:nbus
                    J21(i,k) = J21(i,k) + V(m)* V(n)*(G(m,n)*cos(del(m)-del(n)) + B(m,n)*sin(del(m)-del(n)));
                end
                J21(i,k) = J21(i,k) - V(m)^2*G(m,m);
            else
                J21(i,k) = V(m)* V(n)*(-G(m,n)*cos(del(m)-del(n)) - B(m,n)*sin(del(m)-del(n)));
            end
        end
    end

    J22 = zeros(npq,npq);
    for i = 1:npq
        m = pq(i);
        for k = 1:npq
            n = pq(k);
            if n == m
                for n = 1:nbus
                    J22(i,k) = J22(i,k) + V(n)*(G(m,n)*sin(del(m)-del(n)) - B(m,n)*cos(del(m)-del(n)));
                end
                J22(i,k) = J22(i,k) - V(m)*B(m,m);
            else
                J22(i,k) = V(m)*(G(m,n)*sin(del(m)-del(n)) - B(m,n)*cos(del(m)-del(n)));
            end
        end
    end

    J = [J11 J12; J21 J22];

    X = inv(J)*M;
    dTh = X(1:nbus-1);
    dV = X(nbus:end);

    del(2:nbus) = dTh + del(2:nbus);
    k = 1;
    for i = 1:nbus
        if type(i) == 3
            V(i) = dV(k) + V(i);
            k = k+1;
        else
            V(i)=Vsp(i);
        end
    end

    Iteration = Iteration + 1;
    Tolerance = max(abs(M));

end

P;
Q;
V;
del;


%Loss Calculation: Base power flow case
g = real(Y);
b = imag(Y);
nbranch = length (fb);
Ploss_b = zeros(nbus-1,1);
for m = 1:nbranch
    i = fb(m);
    j = tb(m);
    Ploss_b(m) = (-g(i,j)*((V(i))^2+V(j)^2-2*V(i)*V(j)*cos(del(i)-del(j))));
end

PLoss_b = sum(Ploss_b)*BMva; %Ploss-b is the base loss from the pf equations
PLoss_b

Qloss_b = zeros(nbus-1,1);
for m = 1:nbranch
    i = fb(m);
    j = tb(m);
    Qloss_b(m) = (b(i,j)*((V(i))^2+V(j)^2-2*V(i)*V(j)*cos(del(i)-del(j))));
end
QLoss_b = sum(Qloss_b)*BMva; %Qloss-b is the base loss from the pf equations
QLoss_b


Ploss_del = zeros(nbus,1);
for i = 1:(nbus)
    for j = 1:(nbus)
        if i~=j
            for n = 1:nbus
                Ploss_del(n) = Ploss_del(n) + (2*g(i,j)*V(i)*V(j)*sin(del(i)-del(j)));
            end
        end
    end
end

Qloss_del = zeros(nbus,1);
for i = 1:(nbus)
    for j = 1:(nbus)
        if i~=j
            for n = 1:nbus
                    Qloss_del(n) = Qloss_del(n) + (-2*b(i,j)*V(i)*V(j)*sin(del(i)-del(j)));
            end
        end
    end
end

Ploss_V = zeros(nbus,1); %derivative of P wrt to voltage
for i = 1:(nbus)
    for j = 1:(nbus)
        if i~=j
            for n = 1:nbus
                Ploss_V(n) = Ploss_V(n) + (2*g(i,j)*(V(i)-V(j)*cos(del(i)-del(j))));
            end
        end
    end
end

Qloss_V = zeros(nbus,1);
for i = 1:(nbus)
    for j = 1:(nbus)
        if i~=j
            for n = 1:nbus
                Qloss_V(n) = Qloss_V(n) + (-2*b(i,j)*(V(i)-V(j)*cos(del(i)-del(j))));
            end
        end
    end
end


Pld = Ploss_del(1:end-1,1); Plv = Ploss_V(1:end-1,1);

Qld = Qloss_del(1:end-1,1); Qlv = Qloss_V(1:end-1,1);

N = [Pld; Plv];
M = [Qld; Qlv];

J_apq = inv(J')*N;
J_rpq = inv(J')*M;



for k=1:nbus-1
    J_ap(k) = J_apq(k);
end

for k=1:nbus-1
    J_aq(k) = J_apq((nbus-1)+k);
end


for k=1:nbus-1
    J_rp(k) = J_rpq(k);
end

for k=1:nbus-1
    J_rq(k) = J_rpq((nbus-1)+k);
end


Pg = Pg(1:end-1);
Qg = Qg(1:end-1);
deltaPloss = Pg*J_ap+Qg*J_aq;

deltaPloss=deltaPloss*BMva;
fprintf ('deltaPloss first order %0.2f\n', deltaPloss(dgbus,dgbus));
Ploss = (PLoss_b + deltaPloss);
Error_fo_p = abs ((Ploss-PLoss_b)/PLoss_b)*100;
fprintf ('Error first order P %0.2f\n', Error_fo_p(dgbus,dgbus));

deltaQloss = Pg*J_rp + Qg*J_rq;
deltaQloss=deltaQloss*BMva;
fprintf ('deltaQloss first order %0.2f\n', deltaQloss(dgbus,dgbus));
Qloss = (QLoss_b + deltaQloss);
Error_fo_q = abs ((Qloss-QLoss_b)/QLoss_b)*100;
fprintf ('Error first order Q %0.2f\n', Error_fo_q(dgbus,dgbus));

% 2nd
Plossdd = zeros(nbus,nbus);
for i = 1:nbus
    for j = 1:nbus
        if (i == j)
            for n = 1:nbus
                if (i~=n)
                    Plossdd(i,j) =  Plossdd(i,j)+(2*V(i)*V(n)*g(i,n)*cos(del(i)-del(n)));
                end
            end
        else
            Plossdd(i,j) = -2*g(i,j)*V(i)*V(j)*cos(del(i)-del(j));
        end
    end
end

PlossVV = zeros(nbus,nbus);
for i = 1:nbus
    for j = 1:nbus
        if (i == j)
            for n = 1:nbus
                if (i~=n)
                    PlossVV(i,j) =  PlossVV(i,j) + 2*g(i,n);
                end
            end
        else
            PlossVV(i,j) = -2*g(i,j)*cos(del(i)-del(j));
        end
    end
end

clear k m j i n
PlossdV = zeros(nbus,nbus);
for i = 1:nbus
    for j = 1:nbus
        if (i == j)
            for n = 1:nbus
                if (i~=n)
                    PlossdV(i,j) = PlossdV(i,j)+2*g(i,n)*V(n)*sin(del(i)-del(n));
                end
            end
        else
            PlossdV(i,j) = -2*g(i,j)*V(j)*sin(del(i)-del(j));
        end
    end
end

PlossVd = zeros(nbus,nbus);
for i = 1:nbus
    for j = 1:nbus
        if (i == j)
            for n = 1:nbus
                if (i~=n)
                    PlossVd(i,j) = PlossVd(i,j)+2*g(i,n)*V(n)*sin(del(i)-del(n));
                end
            end
        else
            PlossVd(i,j) = 2*g(i,j)*V(i)*sin(del(i)-del(j));
        end
    end
end

Pdd = Plossdd(2:end,2:end);
PdV = PlossdV(2:end,2:end);
PVd = PlossVd(2:end,2:end);
Pvv = PlossVV(2:end,2:end);


H = [Pdd, PdV; PVd, Pvv];
K = transpose(inv(J))*H*inv(J);


% % k_pp
for i = 1:nbus-1
    for j=1:nbus-1
        k_pp(i,j)=K(i,j);
    end
end

% % k_pq
for i = 1:nbus-1
    for j=1:nbus-1
        k_pq(i,j)=K(i,(nbus-1)+j);
    end
end

% % k_qp
for i = 1:nbus-1
    for j=1:nbus-1
        k_qp(i,j)=K((nbus-1)+i,j);
    end
end

% % k_qq
for i = 1:nbus-1
    for j=1:nbus-1
        k_qq(i,j)=K((nbus-1)+i,(nbus-1)+j);
    end
end

clear Pg Qg
Pg = busd(:,4)/BMva;        % PGi..
Qg = busd(:,5)/BMva;        % QGi..
Pg = Pg(2:end);
Qg = Qg(2:end);
Pg = diag(Pg);
Qg = diag(Qg);
dPloss_so = 0.5*[Pg Qg]*[k_pp,k_pq;k_qp,k_qq]*[Pg;Qg];
dPloss_so=dPloss_so*BMva;
fprintf ('deltaPloss second order %0.2f\n', dPloss_so(dgbus-1,dgbus-1));
Ploss_2 = PLoss_b+deltaPloss(dgbus,dgbus)+abs(dPloss_so(dgbus-1,dgbus-1));
Error_so = abs ((Ploss_2-PLoss_b)/PLoss_b)*100;
fprintf ('Error second order P %0.2f\n', Error_so);

% Second order for Q
Qlossdd = zeros(nbus,nbus);
for i = 1:nbus
    for j = 1:nbus
        if (i == j)
            for n = 1:nbus
                if (i~=n)
                    Qlossdd(i,j) =  Qlossdd(i,j)+(2*V(i)*V(n)*(-b(i,n))*cos(del(i)-del(n)));
                end
            end
        else
            Qlossdd(i,j) = 2*b(i,j)*V(i)*V(j)*cos(del(i)-del(j));
        end
    end
end

QlossVV = zeros(nbus,nbus);
for i = 1:nbus
    for j = 1:nbus
        if (i == j)
            for n = 1:nbus
                if (i~=n)
                    QlossVV(i,j) =  QlossVV(i,j) + (-2*b(i,n));
                end
            end
        else
            QlossVV(i,j) = 2*b(i,j)*cos(del(i)-del(j));
        end
    end
end

QlossdV = zeros(nbus,nbus);
for i = 1:nbus
    for j = 1:nbus
        if (i == j)
            for n = 1:nbus
                if (i~=n)
                    QlossdV(i,j) = QlossdV(i,j)+(-2*b(i,n)*V(n)*sin(del(i)-del(n)));
                end
            end
        else
            QlossdV(i,j) = 2*b(i,j)*V(j)*sin(del(i)-del(j));
        end
    end
end

QlossVd = zeros(nbus,nbus);
for i = 1:nbus
    for j = 1:nbus
        if (i == j)
            for n = 1:nbus
                if (i~=n)
                    QlossVd(i,j) = QlossVd(i,j)+(-2*b(i,n)*V(n)*sin(del(i)-del(n)));
                end
            end
        else
            QlossVd(i,j) = -2*b(i,j)*V(i)*sin(del(i)-del(j));
        end
    end
end

Qdd = Qlossdd(2:end,2:end);
QdV = QlossdV(2:end,2:end);
QVd = QlossVd(2:end,2:end);
Qvv = QlossVV(2:end,2:end);
H_q = [Qdd, QdV; QVd, Qvv];
L = transpose(inv(J))*H_q*inv(J);

% % l_pp
for i = 1:nbus-1
    for j=1:nbus-1
        l_pp(i,j)=L(i,j);
    end
end

% % l_pq
for i = 1:nbus-1
    for j=1:nbus-1
        l_pq(i,j)=L(i,(nbus-1)+j);
    end
end

% % l_qp
for i = 1:nbus-1
    for j=1:nbus-1
        l_qp(i,j)=L((nbus-1)+i,j);
    end
end

% % l_qq
for i = 1:nbus-1
    for j=1:nbus-1
        l_qq(i,j)=L((nbus-1)+i,(nbus-1)+j);
    end
end

clear Pg Qg
Pg = busd(:,4)/BMva;
Qg = busd(:,5)/BMva;
Pg = Pg(1:end-1);
Qg = Qg(1:end-1);
Pg = diag(Pg);
Qg = diag(Qg);
dQloss_so = 0.5*[Pg Qg]*[l_pp,l_pq;l_qp,l_qq]*[Pg;Qg];

dQloss_so=dQloss_so*BMva;
fprintf ('deltaQloss second order %0.2f\n', dQloss_so(dgbus,dgbus));
Qloss_2 = QLoss_b+deltaQloss(dgbus,dgbus)+abs(dQloss_so(dgbus,dgbus));

Error_q = abs ((Qloss_2-QLoss_b)/QLoss_b)*100;
fprintf ('Error second order Q %0.2f\n', Error_q);
