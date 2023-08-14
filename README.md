# 4cycle_invariants
Files for inference of 4-leaf 4-cycle phylogenetic networks from aligned DNA sequence data. Simulation scripts for simulating aligned sequence data along a 4-leaf 4-cycle network under the Jukes-Cantor (JC) and Kimura 2-parameter (K2P) substitution models.

## Dependencies
Scripts are written in Python3 and use the libraries sci-kit bio (requires Python v3.8 or later) and mp-math. 
```console
pip install scikit-bio
pip install mpmath
```

Code to calculate phylogenetic invariants is written in Macaulay2, available from https://www.macaulay2.com/.

The 4-leaf, 4-cycle network from which data is simulated is displayed below. Leaves are labelled by the integers 0,1,2, and 3 in some order. 

![4-cycle network](/4cycle_invariants/docs/assets/images/4cycle.pdf)

## evaluate.py

Evaluates the supplied invariants on all 12 possible 4-leaf 4-cycle networks and writes a score for each to standard output.
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


## simulateK2P.py
As simulateJC.py but simulates according to the K2P substitution model.

## Invariants
The `invariants` directory contains text files describing various phylogenetic invariants. These were computed using the Macaulay2 scripts in the `M2` directory.