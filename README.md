# Keyfitz-Kranzer
The Local Lax-Friedrichs (LLF) scheme applied to a system of Keyfitz-Kranzer type hyperbolic conservation laws.

We consider the Riemann problem for the following non-linear system comprising a conservation law and a balance law with time-dependent source term
```math
      \begin{cases}
          \rho_t + \left( \rho u \left(1 - \left(\frac{\, \rho \,}{\overline{\rho}}\right)^a \right) \right)_x = 0, \\ \\
          \left(\rho u\right)_t + \left(\rho u^2 \left( 1 - \left(\frac{\, \rho \,}{\overline{\rho}}\right)^a \right) \right)_x = a(t) \rho,
      \end{cases}
```
with initial conditions
```math
(\rho, u)(x,0) = \begin{cases}
        (\rho_L,u_L)& x < 0, \\
        (\rho_R,u_R)& x > 0,
    \end{cases}
```
where $\rho$ is a density variable and $u$ a velocity. Under a suitable change of variables, the above system becomes 
```math
\begin{cases}
\rho_t + \left(\rho\left(\tilde{u} + \int_0^ta(s)\text{d}s\right) \left(1 - \left(\frac{\, \rho \,}{\overline{\rho}}\right)^a \right) \right)_x=0, \\ \\
\left(\rho\tilde{u}\right)_t + \left(\rho \tilde{u}\left(\tilde{u} + \int_0^ta(s)\text{d}s\right) \left(1 - \left(\frac{\, \rho \,}{\overline{\rho}}\right)^a \right) \right)_x=0,
\end{cases}
```
under the same initial conditions. The system can be shortened to $\partial_t H + \partial_x G = 0$, where $H$ represents the conserved quantities and $G$ represents the flux. Given a left state $\left(\rho_L,\tilde{u}_L\right)$ and a right state $\left(\rho_R,\tilde{u}_R\right)$, we seek a solution curve in the $\rho\tilde{u}$-plane that shows the connection between the states across the one-dimensional domain. 

The present code is primarily concerned with self-similar solutions to the Riemann problem obtained through the LLF numerical flux-splitting scheme. The Local Lax-Friedrichs (LLF) scheme is a commonly used numerical method for approximating solutions to hyperbolic conservation laws. The scheme uses a discretization of the spatio-temporal domain into grid points at which the conserved quantities $H=\left( \rho\enspace, \rho\tilde{u} \right)^T$ are reconstructed. Letting $\Delta x$ and $\Delta t$ be the corresponding cell dimensions, we can represent any point in our grid as $\left( x_i, t_j \right)=\left( i\Delta x, j\Delta t \right)$. The solution at point $\left(i,j\right)$ can be written similarly as $H_i^j=H(x_i,t_j)$. By splitting the flux $G$ across the left and right spatial neighbors, we can numerically reconstruct the solution using the following formula:
```math
H_{i}^{j+1} = \frac{1}{2}(H_{i-1}^j + H_{i+1}^j)+\frac{CFL}{2\lambda}(G_{i+1}^{j}-G_{i-1}^{j}),
```
where $\lambda := \max_{i} |\lambda_i|$ is the greater of the two characteristic speeds. The Courant number $CFL=\lambda \frac{\Delta t}{\Delta x}$ measures the numerical stability of the scheme and must be chosen such that 
```math
\lambda \frac{\Delta x}{\Delta t} \leq \frac{1}{2}.
```
We will keep $\Delta x=1$ constant so that the CFL condition simply entails requiring $\lambda \Delta t \leq \frac{1}{2}$ throughout the procedure. Note, then, that our spatial grid size is always fixed; however, the time steps between each iteration of the procedure may differ significantly depending on the eigenvalues calculated at each middle state. By obeying the CFL condition, we guarantee that the scheme converges to the physically correct weak solution satisfying the Lax Entropy Condition.

To run this code, make any changes to the initial conditions in the initvars.m file and then run autogen.m. The procedure will produce a series of graphical solutions to the Riemann problem for both $\rho$ and $\tilde{u}$ in terms of the self-similar $x/t$. Additionally, a solution curve will be plotted in $\rho\tilde{u}$-space. 
