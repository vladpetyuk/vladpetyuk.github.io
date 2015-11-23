---
layout: post
title: Lotka-Volterra Model
subtitle: forward and reverse ODE problem solving in R
categories: category1
body-class: categoryclass
tags: [R, ODE, modeling]
---

#Intro
[Lotka-Volterra](https://en.wikipedia.org/wiki/Lotka–Volterra_equations) is a small model that describes a number of biological processes. Perhaps the most know example is describing the populations of prey and predator specie. Although the model is (as low as 2 ODE) it possesses some interesting properties that makes it challending and interesting to model.

#Mode Specifications
Let's consider a system of chemical reaction where a substance `A` converts to `B` through a couple of intermediate steps. The chemical reactions equations are as follows. Essentially it maps to a more common prey/predator as `X` is a prey and `Y` is a predator. The specie `A` would be something like grass or whatever the prey feeds on and `B` is a dead predator. I don't think there is an <u>actual</u> system of chemical reactions that is described by Lotka-Volterra. There are plenty of more complicated ones like [Belousov–Zhabotinsky reaction](https://en.wikipedia.org/wiki/Belousov–Zhabotinsky_reaction). Nonetheless, the relative simplicity of Lotka-Volterra makes it a great toy example.

\\( A + X \xrightarrow{k_1} 2 X \\)

\\( X + Y \xrightarrow{k_2} 2 Y \\)

\\( Y \xrightarrow{k_3} B \\)

<!-- can't control size this way
![img1](/blog/figs/2015-11-21-lotka_volterra/lv.png)
-->
<!-- this way has more handle on figure position and size
     Jekyll will look for it in
     http://localhost:4000/blog/figs/2015-11-21-lotka_volterra/lv.png
     Knitr, however does not compile the document.  I need to sort out 
     the paths. -->

<!-- http://localhost:4000/figs/2015-11-21-lotka_volterra/lv.png -->

<left>
<img src="/blog/figs/2015-11-21-lotka_volterra/lv.png" alt="None" width="300">
<!-- <img src="../figs/2015-11-21-lotka_volterra/lv.png" alt="None" width="300"> -->
<br>
<em>Figure made in Systems Biology Workbench</em>
</left>

#ODEs
The set of ordinary differential equations describing the system can be written as follows.

\\( \frac{dA}{dt} = -k_{1}AX \\)

\\( \frac{dX}{dt} = k\_{1}AX - k_{2}XY \\)

\\( \frac{dY}{dt} = k\_{2}XY - k_{3}Y \\)

\\( \frac{dB}{dt} = k_{3}Y \\)

