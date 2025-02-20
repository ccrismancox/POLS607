

from sympy import *
init_printing(pretty_print=True)
Ti, s2a, s2e = symbols("T_i, sigma^2_\\alpha, sigma^2_\\varepsilon", positive=True)
y,x,theta,ybar,xbar=symbols("y_i,X_i,theta,\\bar{y}_i, \\bar{X}_i",real=True)


ln2_det_i = log(Ti*s2a + s2e)/2 + (Ti-1)/2 * log(s2e)
omega_i = 1 -sqrt(s2e)/sqrt(Ti*s2a+s2e)
ehat = y-x*theta
ebar = ybar-xbar*theta
etilde = (ehat-omega_i*ebar)/sqrt(s2e)

s2=symbols("sigma^2", real=True)


LL = -ln2_det_i - (etilde**2)/2

Dtheta = LL.diff(theta)
Ds2a = LL.diff(s2a)
Ds2e = LL.diff(s2e)

Dtheta = Dtheta.subs([(omega_i, "omega_i"), 
                                     (ehat, symbols("varepsilon_i")),
                                      (ebar, symbols("\\bar{\\varepsilon}_i"))])
Ds2a = Ds2a.subs([(omega_i, "omega_i"), 
                                     (ehat, symbols("varepsilon_i")),
                                     (ebar, symbols("\\bar{\\varepsilon}_i")),
                                     (Ti*s2a+s2e, "sigma^2")])
Ds2e = Ds2e.subs([(omega_i, "omega_i"), 
                                     (ehat, symbols("varepsilon_i")),
                                     (ebar, symbols("\\bar{\\varepsilon}_i")),
                                     (Ti*s2a+s2e, "sigma^2")])


Matrix([[Dtheta], [Ds2a], [Ds2e]])

print(latex(Dtheta))
print(latex(Ds2a))
print(latex(Ds2e))

'''
This section is for if you wanted to explore the
second derivatives (Hessian)

LL.diff(theta).diff(theta).subs([(omega_i, "omega_i"), 
                                     (ehat, symbols("epsilon")),
                                      (ebar, symbols("\\bar{\\epsilon}"))])

LL.diff(theta).diff(s2a).subs([(omega_i, "omega_i"), 
                                     (ehat, symbols("epsilon")),
                                      (ebar, symbols("\\bar{\\epsilon}")),
                                      (Ti*s2a+s2e, "sigma^2")])

LL.diff(theta).diff(s2e).subs([(omega_i, "omega_i"), 
                                     (ehat, symbols("epsilon")),
                                      (ebar, symbols("\\bar{\\epsilon}")),
                                      (Ti*s2a+s2e, "sigma^2")])




LL.diff(s2a).diff(s2a).subs([(omega_i, "omega_i"), 
                                     (ehat, symbols("epsilon")),
                                     (ebar, symbols("\\bar{\\epsilon}")),
                                     (Ti*s2a+s2e, "sigma^2"),
                                     ((sqrt(s2e)/(2*s2**Rational(3/2))-1/(2*sqrt(s2)*sqrt(s2e))) , "a")])



LL.diff(s2a).diff(s2e).subs([(omega_i, "omega_i"), 
                                     (ehat, symbols("epsilon")),
                                     (ebar, symbols("\\bar{\\epsilon}")),
                                     (Ti*s2a+s2e, s2),
                                     ((sqrt(s2e)/(2*s2**Rational(3/2))-1/(2*sqrt(s2)*sqrt(s2e))) , "a")])
LL.diff(s2e).diff(s2e).subs([(omega_i, "omega_i"), 
                                     (ehat, symbols("epsilon")),
                                     (ebar, symbols("\\bar{\\epsilon}")),
                                     (Ti*s2a+s2e,s2 ),
                                     ((sqrt(s2e)/(2*s2**Rational(3/2))-1/(2*sqrt(s2)*sqrt(s2e))) , "a")])
'''
