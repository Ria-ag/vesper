Overview

This project computationally generates Lewis dot structures for given molecular formulas. It uses Weighted Constraint Satisfaction Problems (WCSPs) to model chemical bonding rules and then 
solves them using Toulbar2, a constraint optimization solver. The goal is to automate the traditionally manual process of drawing Lewis structures, ensuring correct valence, bonding, and 
charge distribution for both covalent and ionic compounds.

Features

Converts molecular formulas into WCSP problem files.
Encodes bonding and valence constraints (including octet rule, formal charges, and resonance handling).
Interfaces with Toulbar2 to find optimized solutions.
Outputs Lewis structures with 2d and 3d rendering by using Cytoscape.
Supports a variety of molecules, including: Simple covalent molecules (e.g., CH₄, H₂O), complex covalent molecules (e.g., CH₃OH), ionic compounds (e.g., NaCl).
Currently, some complex molecules do not always produce valid structures as resonance is not fully implemented. 

