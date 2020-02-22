K = QQ[a1,a2,a3,a4,a5,a6,a7,a8,c1,c2]
R = K[t1,t2,t3]

I=ideal(
a1*a5*t1*t2+a2*a5*t1*t2*t3+a3*a6*t1^2*t2*t3-a1*a5*c1,
a1*a5*t2+a3*a8*t1^2*t2*t3+a4*a5*t1*t2*t3-a1*a5*c2,
a1*a5*t3-a1*a5*t1*t2+a3*a5*t1*t3+a3*a7*t1^2*t3-a1*a5*t2-a1*a5*(1-2*c1-2*c2)
)

