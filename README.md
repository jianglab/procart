# ProCart
ProCart is a Web app that plots a cartoon to illustrate the residue properties of amyloid proteins in the context of their atomic structures.

Click the link (<a href="https://protein-structure-procart.herokuapp.com">ProCart</a>) to plot a cartoon of your amyloid protein structure!</a>

<a href="https://protein-structure-procart.herokuapp.com"><img src="./procart.png" style='width: 100%; object-fit: contain'></a>


## Features
* The backbone is represented as thin line segments connected at the Cα atoms
* The β-strands are represented by thicker lines and an arrow ended at the last residue of the strand
* You can choose a coloring scheme from [three choices (Cinema, Lesk, and Clustal)](https://www.bioinformatics.nl/~berndb/aacolour.html)
* To accurately depict the positions of the residues, the residues are presented as circles positioned at the center-of-mass of the entire residue (including backbone and side-chain atoms)
* To quantitatively illustrate the packing of the residues, the sizes of the circles are proportional to the radius-of-gyration computed for each of the residues
* Mouse hovering over the residues will display the coordinates and the amino acid identity in the tooltips
* The plots can be shared/reproduced via a URL displayed in the browser address bar or a QR code displayed below the plots