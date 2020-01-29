---
title: 'PDP: a C++ package for data processing at Scale'
tags:
  - C++
  - Fluid Dynamics
  - Linear Algebra
  - MPI 
  - ScaLAPACK
authors:
  - name: Nicholas Arnold-Medabalimi
    orcid: 0000-0003-2810-7996
    affiliation: 1
affiliations:
 - name: Nicholas Arnold-Medabalimi, Ph.D. Candidate, University of Michigan
   index: 1=
date: 28 January 2020
bibliography: paper.bib

# Optional fields if submitting to a AAS journal too, see this blog post:
# https://blog.joss.theoj.org/2018/12/a-new-collaboration-with-aas-publishing
#aas-doi: 10.3847/xxxxx <- update this with the DOI from AAS once you know it.
#aas-journal: Astrophysical Journal <- The name of the AAS journal.
---

# Summary

In the scientific model reduction community a variety of 
decomposition techniques are used for analysis of complex
spatio-temporal systems as well as reduced-order model 
development. Application of these methods to more complex
problems require memory usage that often exceeds that of 
single machines. Parallel Data Processor (PDP) is a MPI 
C++ package built upon the ScaLAPACK toolset of distributed
Linear Algebra. 



# Statement of need

In the current development of projection based reduced-order models a 
basis is generated from a "data matrix". This matrix is sized as 
(degrees of freedom x number of timesteps). For development problems 
this is often small enough to be done in the scritping environement. 
Additionally in the cases where one of the given dimensions of small
method of snapshots() can be utilized. (however at wall time costs
due to file I/O) However as attempts are made to apply these methods 
to large problems other tools must be used to overcome the memory 
and I/O problem. PDP aims to make overcoming these issues easier
without users needing to understand the depths of distributed
linear algebra computing.




# Mathematics

Single dollars ($) are required for inline mathematics e.g. $f(x) = e^{\pi/x}$

Double dollars make self-standing equations:

$$\Theta(x) = \left\{\begin{array}{l}
0\textrm{ if } x < 0\cr
1\textrm{ else}
\end{array}\right.$$


# Citations

Citations to entries in paper.bib should be in
[rMarkdown](http://rmarkdown.rstudio.com/authoring_bibliographies_and_citations.html)
format.

For a quick reference, the following citation commands can be used:
- `@author:2001`  ->  "Author et al. (2001)"
- `[@author:2001]` -> "(Author et al., 2001)"
- `[@author1:2001; @author2:2001]` -> "(Author1 et al., 2001; Author2 et al., 2002)"

# Figures

Figures can be included like this: ![Example figure.](figure.png)

# Acknowledgements

We acknowledge contributions from Brigitta Sipocz, Syrtis Major, and Semyeong
Oh, and support from Kathryn Johnston during the genesis of this project.

# References