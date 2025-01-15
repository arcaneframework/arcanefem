# Solving Acoustics with ArcaneFEM #

**WARNING: ill posed problem**

<img width="300" align="left" src="https://github.com/user-attachments/assets/430570ff-5185-4331-992f-efebbb23bb91"/>
Single Frequncey acoustic probelm reads: 

$$\(\frac{k}{c}\)^2 u + \nabla u = 0 \quad \forall (x,y,z)\in \Omega$$

$$\frac{\partial u}{\partial\textbf{n}}  = g  \quad \forall (x,y,z)\in \partial{\Omega}_{\text{N}}$$

here $u$ is infact solution of Helmholtz’s problem, for certain values of $\frac{k}{c}$ this problem is ill posed, _resonance_
