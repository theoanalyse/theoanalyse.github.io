---
layout: post
use_math: true
title:  "On the study of Turing instabilities in the $p$ cells reaction-diffusion model with two morphogens"
date:   2022-10-09 22:53:16 +0200
categories: project
---


## Some biological context

Nature is known to present patterns, resulting from a perfect balance between regularity and chaos. Probably the most classical example one can give is the strangely regular distribution of petals and pistils in flowers. Based on this observations, mathematicians of the world started to try finding models that could describe best the emergence of these patterns. Some people see a close link with the Fibonacci sequence and the golden ratio $\phi \approx 1.618$, some see a connection to prime numbers' distribution but the existence or not of these links is beyond the scope of this report. We also observe the appearance of patterns in the animal. Whether it is the skin of the zebrafish, the stripes of the zebra, the spots of the giraffe, the back of the ladybug, or a mixture of spots and stripes in the leopard... It would seem that the model introduced in 1952 by Alan Turing, in his paper "The chemical basis of Morphogenesis", covers all these pigmentation processes involved in the generation of all this fauna of shapes. The goal of this project is, thus, to present and study the 1-dimensional, 2-morphogens, reaction-diffusion equation across $p$-cells arranged on a torus of length $L$, as proposed in 1952.

$$
\begin{equation}
\begin{cases}
\frac{\partial}{\partial t} X (t, s) = D_X \frac{\partial^2}{\partial s^2} X (t, s) + f \bigl( X(t,s), Y(t,s) \bigr), \quad (t, s) \in \Omega
\\ 
\frac{\partial}{\partial t} Y (t, s) = D_Y \frac{\partial^2}{\partial s^2} Y (t, s) + g \bigl( X(t,s), Y(t,s) \bigr), \quad (t, s) \in \Omega
\end{cases}
\end{equation}
$$

where $\Omega = \mathbb{R}^+ \times [0, L]$ and <font color='blue'>$X$</font>, <font color='red'>$Y$</font> denote the quantity of <font color='blue'>activators</font> and <font color='red'>inhibitors</font> respectively. $D_X$ and $D_Y$ are positive real numbers representing diffusion coefficients of the model. Finally, $f$ and $g$ are two functions that defines the reaction between $X$ and $Y$. In our case for Turing's model, $f$ and $g$ are chosen as follows

$$
\begin{align}
f(X, Y) &= 5X - 6Y + 1 \\ g(X, Y) &= 6X - 7Y + 1
\end{align}
$$

## Towards a discretization of the equation using  finite differences

Currently, the model is described using a system of PDE, which we cannot simply simulate. Hence, our way to go is by descretizing the diffusion term using finite differences. We start off by splitting our periodic domain $[0, L]$ into $p$ intervals that are going to define our $p$ cells. for $1 \le r \le p$, we let $X_r$, $Y_r$ denote the quantity of activator and inhibitor in cell $r$.

**Remark:** Since our domain has periodic boundary condition, we have that $X_{p+1} = X_1$, and $Y_{p+1} = Y_1$.

It is only left to discretize our diffusion term using two, third-order Taylor expansion around $s$ (one above $s$, one below). This results in the common expression of the 1D discrete Laplacian operator (modulo some term in $\mathcal o (h^4)$)

$$\frac{\partial^2 X_r}{\partial s^2} = \frac{X_{r-1} - 2X_r - X_{r+1}}{h^2}$$

Assuming that $h = 1$ for now, we thus land on the final system of $2p$ equations at the heart of our study.


$$
\begin{equation}
\begin{cases}
\frac{\partial}{\partial t} X_r = D_X (X_{r-1} - 2X_r + X_{r+1}) + 1 + 5X_r - 6Y_r, \quad 1 \le r \le p
\\ 
\frac{\partial}{\partial t} Y_r = D_Y (Y_{r-1} - 2Y_r + Y_{r+1}) + 1 + 6X_r - 7Y_r, \quad 1 \le r \le p
\end{cases}
\end{equation}
$$


## Parameters and building stones of the program

The code presented in this report has a very linear structure and will be assembled step by step through each chunk of code. The first natural step is to import packages that will be mandatory for the well-being of our program


{% highlight python %}
import numpy as np # obviously.
from numpy import exp, cos, sin, pi # some shortcuts

import matplotlib.pyplot as plt # plotting and visuals
from matplotlib.lines import Line2D # custom legend for curves
from scipy.integrate import odeint # numerical integration for ODE simulation
plt.style.use("ggplot") # visual style that I personaly like :)

# restore classic colormap
from cycler import cycler
plt.rcParams['axes.prop_cycle'] = cycler(color='bgrcmyk')

# adjust figures' size
plt.rcParams['figure.figsize'] = (16, 9)
{% endhighlight %}

