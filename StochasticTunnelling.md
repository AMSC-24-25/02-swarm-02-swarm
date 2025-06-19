#  Stochastic Tunnelling Approach
Stochastic Tunnelling (ST) is an optimization algorithm designed for determining the global minima of complex and rugged energy landscapes. It is a generic physically motivated generalization of
simulated annealing. This approach circumvents the freezing problem which arises when the energy difference between “adjacent” local minima on the energy surface is much smaller than the energy of 
intervening transition states separating them.The physical idea behind the stochastic tunneling method is to allow the particle to “tunnel” forbidden regions, once it has been determined
that they are irrelevant for the low-energy properties of the problem. This can be accomplished by applying the transformation:

$$
f_{\text{STUN}}(x) = 1 - \exp[-\gamma(f(x) - f_0)],
$$

where $f_0$ is the lowest minimum encountered thus far. The effective potential preserves the locations of all minima, but maps the entire energy space from $f_0$ to the maximum
of the potential onto the interval [0, 1]. The degree of steepness of the cutoff of the high-energy regions is controlled by the tunneling parameter $\gamma$.

