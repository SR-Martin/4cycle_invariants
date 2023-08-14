# 4cycle_invariants
Files for inference of 4-leaf 4-cycle phylogenetic networks from aligned DNA sequence data. Simulation scripts for simulating aligned sequence data along a 4-leaf 4-cycle network under the Jukes-Cantor (JC) and Kimura 2-parameter (K2P) substitution models.

## Dependencies
Scripts are written in Python3 and use the libraries sci-kit bio (requires Python v3.8 or later) and mp-math. 
```console
pip install scikit-bio
pip install mpmath
```

Code to calculate phylogenetic invariants is written in Macaulay2, available from https://www.macaulay2.com/.

The 4-leaf, 4-cycle network from which data is simulated is displayed below. Leaves are labelled by the integers 0,1,2, and 3 in some order, e.g. the leaf-labelled network below is denoted (0,1,2,3).

<img src="/docs/images/4cycle.png" alt="4-cycle network" width="250" height="250">

### Installation

The software does not require installation, simply install the dependencies above, clone the repository, and run the scripts.

```console
git clone https://github.com/SR-Martin/4cycle_invariants.git
cd 4cycle_invariants
```

There are two example datasets with the repository. The first is a simulated alignment from the network (0,1,2,3) under the K2P substitution model. Running

```console
python evaluate.py -i invariants/4LeafK2P_GB_deg3.txt -a example_data/network_0123_K2P_100000.phylip
```

should give (after a couple of minutes)

```console
1: (Taxon0,Taxon1,Taxon2,Taxon3)	7.27211e-9
2: (Taxon3,Taxon0,Taxon1,Taxon2)	1.11965e-8
3: (Taxon1,Taxon0,Taxon3,Taxon2)	1.36818e-8
4: (Taxon2,Taxon1,Taxon0,Taxon3)	1.70732e-8
5: (Taxon0,Taxon2,Taxon1,Taxon3)	1.66039e-7
6: (Taxon2,Taxon0,Taxon3,Taxon1)	1.79623e-7
7: (Taxon3,Taxon0,Taxon2,Taxon1)	1.80426e-7
8: (Taxon1,Taxon2,Taxon0,Taxon3)	1.8475e-7
9: (Taxon0,Taxon1,Taxon3,Taxon2)	2.64079e-7
10: (Taxon3,Taxon1,Taxon0,Taxon2)	2.77587e-7
11: (Taxon2,Taxon0,Taxon1,Taxon3)	2.92029e-7
12: (Taxon1,Taxon0,Taxon2,Taxon3)	3.09538e-7
```
Note that the score for the network that generated the alignment is an order of magnitude lower than the others.

Evolution that has been tree like is detected automatically. The command

```console
python evaluate.py -i invariants/4LeafK2P_GB_deg3.txt -a example_data/tree_01-23_100000.phylip
```
gives

```console
1: (Taxon1,Taxon0,Taxon3,Taxon2)	2.72771e-9
2: (Taxon0,Taxon1,Taxon2,Taxon3)	2.90057e-9
3: (Taxon3,Taxon0,Taxon1,Taxon2)	3.13039e-9
4: (Taxon1,Taxon0,Taxon2,Taxon3)	3.14916e-9
5: (Taxon0,Taxon1,Taxon3,Taxon2)	3.23117e-9
6: (Taxon2,Taxon1,Taxon0,Taxon3)	3.31187e-9
7: (Taxon2,Taxon0,Taxon1,Taxon3)	3.48568e-9
8: (Taxon3,Taxon1,Taxon0,Taxon2)	4.52525e-9
9: (Taxon2,Taxon0,Taxon3,Taxon1)	2.34451e-7
10: (Taxon1,Taxon2,Taxon0,Taxon3)	2.41453e-7
11: (Taxon0,Taxon2,Taxon1,Taxon3)	2.56563e-7
12: (Taxon3,Taxon0,Taxon2,Taxon1)	2.63093e-7
Tree-like evolution detected with tree ((Taxon2,Taxon3),(Taxon0,Taxon1)).
```

## evaluate.py

Evaluates the supplied invariants on all 12 possible 4-leaf 4-cycle networks and writes a score for each to standard output, as above. The script attempts to detect evolution that has been tree-like, and identifies the 4-leaf tree.
```console
python evaluate.py -i invariants/4LeafK2P_GB_deg3.txt -a /path/to/alignment.phylip
```

## simulateJC.py

Simulates multiple sequence alignment data in phylip format, according to the JC msubstition model, on the specified 4-leaf 4-cycle network. Substitution parameters and the tree-ratio can either be specified or randomly generated, and are written to standard output. The `-r` flag accepts "low" (0.1% chance of substitution along edge) or "medium" (1% chance of substitution along edge) when randomly generating substitution parameters. 

```console
python simulate.py -l <MSA length> -s <seed> -p <edge parameters> -g <gamma parameter> -o <output filename>
-l <MSA length>		 Length of MSA to output
-s <seed>		 Integer value for seeding random numbers
-p <edge parameters>	 Comma-separated list of values giving the probability of a substitution along each edge.
-g <gamma parameter>	 Floating point value giving the probability of a site evolving along edge e in the network.
-o <output filename>	 Filename for output MSA in phylip format
```

### Examples
To simulate an alignment of 10,000 base-pairs with tree-ratio 0.5 and randomly generate substitution parameters with random seed 1:
```console
python simulateJC.py -n '(0,1,2,3)' -g 0.5 -r medium -l 10000 -s 1 -o alignment_1.phylip
```

To simulate an alignment of 1,000,000 base-pairs with tree-ratio 0.75 and defined substitution parameters:
```console
python simulateJC.py -n '(0,1,2,3)' -g 0.75 -p 0.01,0.013,0.02,0.005,0.01,0.006,0.007,0.015 -l 1000000 -s 1 -o alignment_2.phylip
```
Here the edge parameters are in the order a,b,c,d,e,f,g,h (as in the image above).

## simulateK2P.py
As simulateJC.py but simulates according to the K2P substitution model.

## Invariants
The `invariants` directory contains text files describing various phylogenetic invariants. These were computed using the Macaulay2 scripts in the `M2` directory.