# Notes

## 2/16/26

Now onto the fun stuff.

Group $G$ with $a^\nu, b^\nu\dots \in G$, with $\nu = 1, \dots \nu$ group dimension.

For the bilinear operation with this $\nu$ term we'll use $a \times b = c$, but $c^\nu = \phi^\nu(a,b)$

We say that $\phi$ is a set of functions with closure.

1) Unit element $a_0$

$$
\phi^\nu (a_0,a)=\phi^\nu(a,a_0) = a^\nu
$$ 

Previously, we proved the left unit is the same as the right unit. Actually, we'll use $a_0 = 0$ going forward.

2) Inverse element

$$
\phi^\nu (\bar a, a) = \phi^\nu(a, \bar a) = 0^\nu
$$

3) Associativity

$$
\phi(\phi(a,b),c)=\phi(a,\phi(b,c))
$$

We want the Lie Algebra structure constants to tell us more about the group. 

Lie algebra picture: $\delta a^\nu$ near the unit element, $d a^\nu$ near the element $a^\nu$. This yields

$$
\phi^\nu (a, \delta a) = a^\nu + d a^\nu % and hope da^\nu isn't zero
$$

This comes from $\delta a$ being a tiny bit on top of the unit element, so it only changes $a$ by a small $da^\nu$

Generators idea:

$$
\phi^\nu(a,\delta a) = \phi^\nu(a,0) + \left. \partial \phi^\nu (a,b) \over \partial b^\nu \right|_{b=0} \partial a^\beta
$$

vierbein / tetrad time. Kinda a tangent. This is an alternative formalism to $g$.

$$
    \mu_{(\beta)}^\nu (a) = \left. \partial \phi^\nu (a,b) \over \partial b^\beta \right|_{b=0}
\\  da^\nu = \mu_{(\beta)}^\nu(a) \delta a^{(\beta)}
$$

Don't freak out, $\beta$ is a dummy index. Also, we use parenthesis to indicate that we're near $0$.

---

Fermions are stuck in special relativity, so we multiply by $\mu$ to get generalized out of flat minkowski space.

$$
    \mu_{(\alpha)}^\nu \psi^{(\alpha)}
\\  g^{\mu \nu} = \xi^{(\alpha)(\beta)} \mu_{(\alpha)}^\mu (a) \mu_{(\beta)}^{\nu} (a)
$$

Tetradized. Okay, how about any general function $F(a)$

$$
      dF(a) 
    = F(a+da) - F(a)
\\  = \left.\partial F \over \partial a^\alpha\right|_{\alpha=0} da^\alpha
\\  = {\partial F \over \partial a^\alpha} \mu^\alpha_{(\beta)} (a) \delta a^{(\beta)}
\\  = \delta a^{(\beta)} X_{(\beta)}(a)F
\\    X_{(\beta)}(a) = \mu_{(\beta)}^\alpha {\partial \over \partial a^\alpha}
$$

Mathematicians don't care about $i$, so $X_{(\beta)}$ is our generator. 

$$
\left[ X_{(\alpha)}, X_{(\beta)} \right] = C_{(\alpha)(\beta)}^{(\gamma)} X_{(\gamma)}
$$

---

Near $0$ with a nice function $\phi$

$$
\begin{align*}
      \phi^\nu(a,b) 
   &= \phi^\nu ( 0,0) 
    + n^\nu_\kappa a^\kappa 
    + m^\nu_\kappa b^\kappa
    + {f'}_{\kappa\lambda}^\nu a^\kappa a^\lambda 
    + {f''}_{\kappa\lambda}^\nu b^\kappa b^\lambda
    + f_{\kappa \lambda}^\nu a^\kappa b^\lambda
\\ &+ {g'}^{\nu}_{\kappa \lambda \sigma} a^\kappa a^\lambda a^\sigma 
    + g^\nu_{\kappa \lambda \sigma}a^\kappa a^\lambda b^\sigma
    + h^\nu_{\kappa \lambda \sigma} b^\kappa b^\lambda b^\sigma
    + h^\nu_{\kappa \lambda \sigma} b^\kappa b^\lambda a^\sigma 
    + \dots
\end{align*}
$$
Multivariable calculus. Let's get to it.

$$
    \phi(0,0) = 0
\\  \phi^\nu(0,0) = 0
$$
which takes out half the terms.

$$
    \phi^\nu (a,0)=a^\nu
\\  \phi^\nu(0,b)=b^\nu
$$
which takes out an additional
$$
    n^\nu_\kappa = \delta^\nu_\kappa, \quad m_\kappa^\nu = \delta_\kappa^\nu
\\  f' = f'' = g' = h' = 0
$$
yielding
$$
      \phi^\nu(a,b) 
    = a^\nu 
    + b^\nu 
    + f^\nu a^\kappa b^\lambda
    + g^\nu_{\kappa \lambda \sigma} a^\kappa a^\lambda b^\sigma
    + h^\nu_{\kappa \lambda \sigma} b^\kappa b^\lambda a^\sigma
$$ 

Next up: $\phi^\nu(a,\phi(b,c))-\phi^\nu(\phi(a,b),c)=0$. Two pages later, we have

$$
      f_{\rho \xi}^\nu f_{\kappa \lambda}^\xi - f_{\xi \lambda}^\nu f_{\rho \kappa}^\xi
    = g_{\rho \kappa \lambda}6\nu + g_{\kappa \rho \lambda}^\nu - \left(h_{\kappa \lambda \rho}^\nu + h_{\lambda \kappa \rho}^\nu \right)
$$
And we can use jacobi to get
$$
    C_{\kappa \lambda}^\rho f_{\kappa \lambda}^\rho-f_{\lambda \kappa}^\rho
\\  C_{\rho \xi}^\nu C_{\kappa \lambda}^\xi + C_{\lambda \xi}^\nu C_{\rho \kappa}^\xi + C_{\kappa \xi}^\nu C_{\lambda\rho}^\xi = 0
$$

This is all General Relativity too btw.

## 2/18/26