---
layout: default
title: Home
---


# Théo André - Welcome to my Personal Page

<!-- <center>
<img src="./images/photo-theo.jpg" width="300" height="300">
</center> -->


Hi there! I'm a 3rd year PhD student at Heidelberg University. I work on Reaction-Diffusion-ODE systems and mathematical modeling of stem cells dynamics, under the supervision of [Prof. Dr. Anna Marciniak-Czochra](https://biostruct.iwr.uni-heidelberg.de/compactseminar_res_bib.php).

## Research projects

### Reaction-Diffusion-ODE

We study the phenomenon of pattern formation and symmetry breaking in coupled dynamics involving both diffusing and non-diffusing components. Such problems have the following strong form

$$\begin{aligned} \dfrac{\partial}{\partial t} u(x, t) &= F \bigl(u(x, t), v(x, t) \bigr) & x \in \Omega, t \ge 0, \\ \dfrac{\partial}{\partial t} v(x, t) &= D\Delta v(x, t) + G \bigl( u(x,t), v(x, t) \bigr) & x \in \Omega, t \ge 0, \\ \bigl( u(\mathbf{\cdot}, 0), v(\mathbf{\cdot}, 0) \bigr) &\in \left( L^\infty(\bar \Omega) \right)^{\dim(u)} \times \left(W^{2, p}(\Omega)\right)^{\dim(v)}, \end{aligned},$$

endowed with either Neumann or Periodic boundary conditions on the domain $\Omega = (0, L)$. In more detais, we investigate the existence of stable patterns in a Turing-like sense, and the behavior of solutions near steady-states exhibiting diffusion-driven instability, in a bistable setting. [(see this work)](https://arxiv.org/abs/2511.15648)

### Neural Stem Cells (NSC) dynamics modeling 

Neural Stem Cells (NSC) are a specialized population of cells with the ability to self-renew and differentiate into various neural lineages, such as neurons and glial cells. 
In most adult mammals (e.g. the mouse) it is observed that the number of NSCs declines with age. In adult zebrafish, NSCs maintain lifelong neurogenic activity. Experimental evidence suggests that the preservation of that homeostatic state is the result of a finely tuned regulation mechanism of stem cell fate dynamics.
We investigate what shape these regulation feedbacks may take by means of non-linear ODE models.

## Publications

- [Spatial scale separation and emergent patterns in coupled diffusive-nondiffusive systems](https://arxiv.org/abs/2511.15648) (preprint) -- Théo André, Szymon Cygan, Anna Marciniak-Czochra, Finn Münnich. Link: 
