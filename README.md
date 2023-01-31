# ProCart
ProCart is a Web app that plots a cartoon to illustrate the residue properties of amyloid proteins in the context of their atomic structures.

Click the link ([ProCart](https://procart.streamlit.app)) to plot a cartoon of your amyloid protein structure!

<a href="https://procart.streamlit.app/?chain_ids=A&chain_ids=B&pdb_id=7MKH"><img src="./procart.png" style='width: 100%; object-fit: contain'></a>

## Why ProCart?
While the atomic structures of amyloid fibrils, for example, [the structure of Tau filaments from Alzheimer’s disease](https://www.nature.com/articles/nature23002#Fig3) and [the structure of Tau filaments in Prion protein amyloidoses](https://link.springer.com/article/10.1007%2Fs00401-021-02336-w#Fig4), are often presented nicely in XY plane as colored circles around line segments through the Cα atoms, the figures were manually drawn due to the lack of an automated graphic display tool. To eliminate the tedious manual process and to be more quantitative with the positions and sizes of different residues for our new structures of TMEM106B filaments from MSTD, we developed the ProCart method originally as [a Google Colab notebook](https://colab.research.google.com/drive/10q6qJyccMtJlNWn-7I49g6sRwVf_IkWX?usp=sharing). To make the method more convenient and accessible to all users, we then further developed it into a Web app with some nice features listed below:
## Features
* Multiple methods to input an atomic model: upload from your local disk, a URL, or a PDB ID
* Automatically rotate the model around Z-axis to orient its longest direction along X-axis
* The backbone is represented as thin line segments connected at the Cα atoms
* The β-strands are represented by thicker lines and an arrow ended at the last residue of the strand
* You can choose a coloring scheme from five presets ([Charge](https://www.sigmaaldrich.com/US/en/technical-documents/technical-article/protein-biology/protein-structural-analysis/amino-acid-reference-chart), [Hydrophobicity](https://www.cgl.ucsf.edu/chimera/docs/UsersGuide/midas/hydrophob.html), [Cinema, Lesk, Clustal](https://www.bioinformatics.nl/~berndb/aacolour.html)) or specify your own coloring scheme to map any chain, sequence range or amino acid to any color
* To accurately illustrate the packing of the residues, the residues are presented as the union outline of all atoms in the residue with the size of each atom equal to its van der waals radius.
* If the user chooses to represent each residue as a circle, the circles will be positioned at the center-of-mass of the entire residue (including backbone and side-chain atoms) with the sizes of the circles proportional to the radius-of-gyration computed for each of the residues
* The undulation of the chains perpendicular to the XY plane can be displayed in the Z-plot
* Mouse hovering over the residues will display the coordinates and the amino acid identity in the tooltips
* The plots can be shared/reproduced via a URL displayed in the browser address bar or a QR code displayed below the plots