---
title: '`Peacock.jl`: Photonic crystals in Julia'
tags:
  - Julia
  - physics
  - optics
  - photonics
  - photonic crystals
  - photonic topological insulators
  - Maxwell's equations
  - metamaterials
authors:
  - name: Samuel J. Palmer
    orcid: 0000-0003-0485-7047
    affiliation: "1"
  - name: Vincenzo Giannini
    orcid: 0000-0001-8025-4964
    affiliation: "2"
affiliations:
 - name: The Blackett Laboratory, Imperial College London, London, SW7 2AZ, UK
   index: 1
 - name: Instituto de Estructura de la Materia (IEM), Consejo Superior de Investigaciones Científicas (CSIC), Serrano 121, 28006, Madrid, Spain
   index: 2
date: 21 May 2020
bibliography: paper.bib
---

# Summary

![The irridescent colours of peacocks arise from the nanoscale 'photonic
crystal' structure of the feathers, rather than from pigmentation. The image of
the peacock's photonic crystal structure is reproduced from @zi2003coloration,
Copyright National Academy of Sciences. \label{fig:zoom}](../docs/src/assets/peacock_feathers_zoom.png)

The **P**lane-wave **E**xpansion **A**pproach to **C**haracterising **O**ptical
**C**rystals in **k**-space - otherwise known as
`Peacock.jl` - is a Julia package for studying
photonic crystals using the Plane Wave Expansion Method [@rumpf2006design].
Photonic crystals are materials whose optical properties arise from the
structure of the material [@yablonovitch1987inhibited; @john1987strong], and
`Peacock.jl` is named for the irridescent colours of peacock feathers which
arise not from pigmentation but from their photonic crystal structure
[@zi2003coloration], as shown in \autoref{fig:zoom}.

The response of a photonic crystal is strongest
when the periodicity of the structure is comparable to the wavelength of light.
For visible light, photonic crystals are built from components that are just a
few hundred nanometers in size. Advances in nanofabrication mean that 'designer'
photonic crystals can now be manufactured for unprecedented control over the
flow of light, with applications ranging from optical fibers
[@knight2003photonic] to photonic circuitry [@joannopoulos2008molding].

`Peacock.jl` provides a user-friendly interface to analyse 2D photonic crystals
with support for non-orthogonal unit cells and inhomogeneous permittivity and/or
permeability. The package focuses on ease of use and includes a `Zoo` submodule
of crystals from published works. We hope that in time the `Zoo` will grow and
aid others to reproduce and extend existing work.

Photonic crystals are also a promising platform for exotic materials known as
photonic topological insulators [@wu2015scheme; @blanco2019engineering;
@rider2019perspective], the photonic analogue of the electronic topological
insulator [@kane2005z; @kane2005quantum; @bernevig2006quantum] for which the
2016 Nobel Prize in Physics was awarded. `Peacock.jl` includes built-in methods
to identify topological photonic crystals using Wilson loops
[@blanco2020tutorial], for which open source implementations have been lacking.

![Examples. (a) ... (b) ..... \label{fig:examples}](examples.svg)

# Acknowledgements

S.J.P. acknowledges his studentship from the Centre for Doctoral Training on
Theory and Simulation of Materials at Imperial College London funded
by EPSRC Grant No. EP/L015579/1.
​
V.G. acknowledges the Spanish Ministerio de Economia y Competitividad for
financial support through the grant NANOTOPO (FIS2017-91413-EXP) and also the
Ministerio de Ciencia, Innovació n y Universidades through the grant MELODIA
(PGC2018-095777-B-C21).

# References