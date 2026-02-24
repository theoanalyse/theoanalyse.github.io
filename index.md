---
layout: default
title: Home
---

## About Me
Hi there! I'm a 3rd year PhD student at Heidelberg University. I work on Reaction-Diffusion-ODE systems and mathematical modeling of neural stem cells dynamics, under the supervision of [Prof. Dr. Anna Marciniak-Czochra.](https://biostruct.iwr.uni-heidelberg.de/compactseminar_res_bib.php)

## Research Projects

### Turing-like patterns in Reaction-Diffusion-ODE systems

We study the phenomenon of pattern formation and symmetry breaking in coupled dynamics involving both diffusing and non-diffusing components. Such problems have the following strong form

$$\begin{aligned} \dfrac{\partial}{\partial t} u(x, t) &= F \bigl(u(x, t), v(x, t) \bigr) & x \in \Omega, t \ge 0, \\ \dfrac{\partial}{\partial t} v(x, t) &= D\Delta v(x, t) + G \bigl( u(x,t), v(x, t) \bigr) & x \in \Omega, t \ge 0, \\ \bigl( u_0(\mathbf{\cdot}), v_0(\mathbf{\cdot}) \bigr) &\in \left( L^\infty(\bar \Omega) \right)^{\dim(u)} \times \left(W^{2, p}(\Omega)\right)^{\dim(v)}, \end{aligned},$$

endowed with either Neumann or Periodic boundary conditions on the domain $\Omega = (0, L)$. In more detais, we investigate the existence of stable patterns in a Turing-like sense, and the behavior of solutions near steady-states exhibiting diffusion-driven instability, in a bistable setting. [(see this work)](https://arxiv.org/abs/2511.15648)

### Neural Stem Cells (NSC) dynamics modeling in adult zebrafish

Neural Stem Cells (NSC) are a specialized population of cells with the ability to self-renew and differentiate into various neural lineages, such as neurons and glial cells. 
In most adult mammals (e.g. the mouse) it is observed that the number of NSCs declines with age. In adult zebrafish, NSCs maintain lifelong neurogenic activity. Experimental evidence suggests that the preservation of that homeostatic state is the result of a finely tuned regulation mechanism of stem cell fate dynamics.
We investigate what shape these regulation feedbacks may take by means of non-linear ODE models.

## Publications

- **Théo André**, Szymon Cygan, Anna Marciniak-Czochra, Finn Münnich:
   “Spatial scale separation and emergent patterns in coupled diffusive-nondiffusive systems”  
   *Preprint*, (2025) [[Paper]](https://arxiv.org/abs/2511.15648)

- Marco David, **Théo André**, Mathis Bouverot-Dupuis, Eva
Brenner, ..., Jonas Bayer:  
   “Universal Pairs for Diophantine Equations”  
   *Isabelle Library*, (2026) [[Library Paper]](https://www.isa-afp.org/browser_info/current/AFP/Diophantine_Universal_Pairs/document.pdf) and [[Associated Paper]](https://arxiv.org/abs/2505.16963)


## List of Talks and Posters

- Poster presentation "Reduction and Diffusion-Driven Instability Analysis In a Receptor-Based Model for Hydra Morphogenesis" at MIMUW, 05/2023, Warsaw. [(Poster)](./talks_and_posters/PosterMIMUW052023.pdf)
- Invited talk "" at I2M, 10/2024, Marseille. [(slides)](./talks_and_posters/TalkI2M122024%20.pdf)
- Talk "Modeling Zebrafish’s Adult Neural Stem Cell
Dynamics with Feedback", 10/2025, Heidelberg. [(slides)](./talks_and_posters/TalkUH102025.pdf)
- Invited talk "Comparing mathematical models of NSC Dynamics in Zebrafish against the Mouse" at DKFZ, 02/2026, Heidelberg. [(slides)](./talks_and_posters/TalkDKFZ022026.pdf)


