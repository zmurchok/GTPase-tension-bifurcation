# GTPase-tension-bifurcation

This repo contains code for a bifurcation analysis of the GTPase-tension model in Zmurchok et al. 2018... see my webpage for links.

## Requirements

Matcont (we used version 7p2, availble at https://sourceforge.net/projects/matcont/files/matcont/matcont7p2/) and MATLAB.

XPPAUT is available at http://www.math.pitt.edu/~bard/xpp/xpp.html.

## Contents

- Fig1.m generates Fig 1A and 1B
- Fig2X.m generates Fig 2X, X = A,B,C
- bettercolors contains some better colors
- SingleGTPaseTension.m contains the functions
- GTPase-tension.ode for use with XPPAUT
- phaseplanes/
    - PPMaker6.m generates the panels for Fig 4
    - dXX/ are folders with .dat (.txt files) that contain information obtained from XPPAUT regarding the stable/unstable manifolds and the periodic solutions that exist in each region

## Contributors

Cole Zmurchok https://zmurchok.github.io
Matthew Sahota
Wayne Nagata https://www.math.ubc.ca/~nagata/
Eric N Cytrynbaum https://www.math.ubc.ca/~cytryn/
