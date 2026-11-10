import sympy as s
x,t=s.symbols('x t');A=[-8,-7,1,5,9]
P=s.prod(x-a for a in A);Q=x*s.prod(x+a for a in A)
F=s.Poly(s.expand(P*(Q+t)),x)
assert s.expand(P)==x**5-110*x**3+2629*x-2520
assert s.expand(Q)==x**6-110*x**4+2629*x**2+2520*x
assert all(F.nth(i)==0 for i in [2,4,6,8,10]) and F.nth(0)==-2520*t
h=729*t**5-6403320*t**4-930770808176*t**3-3787455316296960*t**2+27536385315056910336*t-15173754195676653158400
Delta=s.discriminant(Q+t,x)
assert s.expand(Delta+64*h)==0
assert s.degree(s.gcd(h,s.diff(h,t)),t)==0
assert -h.subs(t,0)==123181793280**2
collisions={int(-Q.subs(x,a)) for a in A}
assert collisions=={40320,35280,-5040,-25200,-45360}
for prime,degrees in [(11,[5]),(19,[2,3])]:
    f=s.Poly(h,t,modulus=prime)
    assert f.degree()==5 and s.degree(s.gcd(f,f.diff()))==0
    factors=s.factor_list(f)[1]
    assert sorted(g.degree() for g,m in factors)==degrees
    assert all(m==1 for g,m in factors)
print('All exact genus-2 pencil assertions passed.')
print('F=',F.as_expr())
print('Correct curve: y^2 =',-h)
print('No nontrivial split fiber or rational-point completeness result claimed.')
