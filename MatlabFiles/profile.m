function [s, sdot ,sdotdot] = profile (pi,pf,t1,t2,tf,is_circ,is_via,rho,alpha,delta)

si=0;
dt=t2-t1;
t=0:0.001:dt;

s_primo = zeros(1,length(t));
sdot_primo = zeros(1,length(t));
sdotdot_primo = zeros(1,length(t));

num = tf*1000 + 1;
s = zeros(1,length(num));
sdot = zeros(1,length(num));
sdotdot = zeros(1,length(num));

%traiettoria circolare
if is_circ == 1
    sf=alpha*rho;
else
    sf=norm(pf-pi);
end

sdot_circ=1.5*abs(sf-si)/dt; 
tc=(si-sf+sdot_circ*dt)/sdot_circ;
sdotdot_circ=sdot_circ/tc;

% trapezoidal profile
for k=1:length(t)
    if t(k)<=tc
        s_primo(k)=si+0.5*sdotdot_circ*t(k).^2;
        sdot_primo(k)=sdotdot_circ*t(k);
        sdotdot_primo(k)=sdotdot_circ;
    elseif (t(k)>tc && t(k)<=dt-tc)
        s_primo(k)=si+sdot_circ*(t(k)-tc/2);
        sdot_primo(k)=sdot_circ;
        sdotdot_primo(k)=0;
    elseif (t(k)>dt-tc && t(k)<=dt)
        s_primo(k)=sf-0.5*sdotdot_circ*(dt-t(k))^2;
        sdot_primo(k)=sdotdot_circ*(dt-t(k));
        sdotdot_primo(k)=-sdotdot_circ;
    end
end

ta=0:0.001:t1;
tb=t1+0.001:0.001:t2;
tc=t2+0.001:0.001:tf;

%via point
if is_via == 0 
    s(1:length(ta))=0;
    s(length(ta)+1:length(ta)+length(tb)+1)=s_primo;
    s(length(ta)+length(tb):length(ta)+length(tb)+length(tc))=sf;
    
    sdot(1:length(ta))=0;
    sdot(length(ta)+1:length(ta)+length(tb)+1)=sdot_primo;
    sdot(length(ta)+length(tb):length(ta)+length(tb)+length(tc))=0;
    
    sdotdot(1:length(ta))=0;
    sdotdot(length(ta)+1:length(ta)+length(tb)+1)=sdotdot_primo;
    sdotdot(length(ta)+length(tb):length(ta)+length(tb)+length(tc))=0;
else
    inc=delta/0.001;
    
    s(1:(length(ta))-inc)=0;
    s(length(ta)+1-inc:length(ta)+length(tb)+1-inc)=s_primo;
    s(length(ta)+length(tb)-inc:length(ta)+length(tb)+length(tc))=sf;
    
    sdot(1:length(ta)-inc)=0;
    sdot(length(ta)+1-inc:length(ta)+length(tb)+1-inc)=sdot_primo;
    sdot(length(ta)+length(tb)-inc:length(ta)+length(tb)+length(tc))=0;
    
    sdotdot(1:length(ta)-inc)=0;
    sdotdot(length(ta)+1-inc:length(ta)+length(tb)+1-inc)=sdotdot_primo;
    sdotdot(length(ta)+length(tb)-inc:length(ta)+length(tb)+length(tc))=0;
end

end