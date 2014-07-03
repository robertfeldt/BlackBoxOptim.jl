p = 2
n = 10000

function ftrefethen(p)
  x = p[1]; y = p[2]
  return exp(sin(50 * x)) + sin(60 * exp(y)) + sin(70 * sin(x)) + 
             sin(sin(80 * y)) - sin(10 * (x + y)) + (x^2 + y^2)/4
end
f = ftrefethen

alpha = 0.602
gamma = 0.101

theta_lo = -100*ones(p,1)
theta_hi = 100*ones(p,1)
theta = 0.40 * (theta_hi - theta_lo) * rand(p, 1)
theta=min(theta,theta_hi)
theta=max(theta,theta_lo)

# Input to gain settings calculations
initial_step_length = 0.40 * maximum(theta_hi - theta_lo)
gavg = 2 # How many gradient estimates per iteration
n_for_gain_calc = n/20
estimated_std_dev_measurement_noise = 0.0

A = 0.10 * n / (2*gavg)
c = max(estimated_std_dev_measurement_noise/(gavg^0.5), 0.0001)

gbar=zeros(p, 1)
for i=1:int(n_for_gain_calc/(2*gavg))
   ghat=zeros(p,1)
   for j=1:gavg
      delta=2*round(rand(p,1))-1;
      thetaplus=theta+c*delta;
      thetaminus=theta-c*delta;
      yplus=f(thetaplus);
      yminus=f(thetaminus);
      ghat=(yplus-yminus)./(2*c*delta)+ghat;
   end
   gbar=gbar+abs(ghat/gavg);
end
meangbar=mean(gbar)/(n_for_gain_calc/(2*gavg))
a = initial_step_length * ((A+1)^alpha)/meangbar;

println("a = $(a), c = $(c), A = $(A)")

#a = 0.0017
#c = 0.00001 # No noise so we can select c small
#A = 5

min_f = Inf
best = 1

for k=0:n-1

    ak = a/(k+1+A)^alpha
    ck = c/(k+1)^gamma

    delta = 2*round(rand(p,1))-1

    thetaplus = theta + ck*delta
    thetaminus = theta - ck*delta

    yplus=f(thetaplus)
    yminus=f(thetaminus)

    ghat = (yplus - yminus)./(2*ck*delta)
    theta = theta - ak * ghat

    theta=min(theta,theta_hi)
    theta=max(theta,theta_lo)

    tf = f(theta)

    if tf < min_f
      min_f = tf
      best = theta
    end

    println("theta = ", theta, ", fitness = ", )

end
println("a = $(a), c = $(c), A = $(A)")
println("best = ", best, ", fitness = ", min_f)