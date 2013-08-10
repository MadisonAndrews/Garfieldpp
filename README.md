About Garfield++
==========

In a nutshell
-------------

Garfield++ is a toolkit for the detailed simulation of detectors which use gases or semi-conductors as sensitive medium. The main area of application is currently in micropattern gaseous detectors.

Garfield++ shares functionality with Garfield. The main differences are the more up-to-date treatment of electron transport in gases and the user interface, which is derived from ROOT.

Elements
--------

**Ionisation**

The Heed program generates ionisation patterns of fast charged particles. The core of Heed is a photo-absorption and ionisation model. Heed in addition provides atomic relaxation processes and dissipation of high-energy electrons.

**Electric fields**

Garfield++ currently offers the following techniques for calculating electric fields:

solutions in the thin-wire limit for devices made of wires and planes;
interfaces with the finite element programs
Ansys,
Elmer, and
CST,
which can compute approximate fields in nearly arbitrary 3-dimensional configurations with dielectrics and conductors;
an interface with the Synopsys Sentaurus device simulation program.
In the future, the program should become interfaced with neBEM, as is already the case for Garfield.

**Transport of electrons**

Magboltz is used to compute electron transport and avalanches in nearly arbitrary gas mixtures.

Authors
-------

Garfield++ includes contributions from

Kim Baraka,
Joshua Renner,
Heinrich Schindler,
Nicholi Shiell,
Rob Veenhof,
Klaus Zenker.
Magboltz is developed by Steve Biagi.

Heed was written by Igor Smirnov.

