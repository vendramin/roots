# Computing finite Weyl groupoids

The repository contains the code corresponding to the paper

```
Computing finite Weyl groupoids  
Iván Angiono, Leandro Vendramin
```

# Description

The file ```roots.g``` does...

Some braidings are saved in the following folders: 
```
rank2/
rank3/
rank4/
rank5/
rank6/
rank7/
rank8/
```
There also logfiles are saved.  To obtain the data corresponding, say, to 
tank two braidings, one uses the file ```rank2.g```. Similarly, 
there are scripts for generating data for other ranks. 

## Rank two examples

The file ```rank2.g``` computes positive roots, Cartan roots and the longest word for each rank two braiding included in the folder ```rank2/```. The results are stored in ```.log``` files in the folder ```rank2/```.    

## A concrete example

The file ```rank2/ufo9d.g``` contains the following braiding, corresponding to $\mathfrak{ufo}_{9d}$. 
```
[ [     E(24),  E(24)^19 ],
  [         1,        -1 ] ]

```
In GAP, the command

```
gap> Read("roots.g");
gap> roots("rank2/ufo9d.g");
```
will produce the file ```rank2/ufo9d.log```, containding the following information: 
```
File: rank2/ufo9d.g
The braiding is:
[ [     E(24),  E(24)^19 ],
  [         1,        -1 ] ]
Longest word:
[ 1, 2, 1, 2, 1, 2, 1, 2 ]
Positive roots (8):
[ [ 1, 0 ], [ 5, 1 ], [ 4, 1 ], [ 3, 1 ], [ 5, 2 ], [ 2, 1 ], [ 1, 1 ], [ 0, 1 ] ]
Cartan roots (2):
[ [ 1, 0 ], [ 5, 2 ] ]
```
This means that $s_1s_2s_1s_2s_1s_2s_1s_2$ is the longest word, 
there are eight positive roots: $(1, 0)$, $(5, 1)$, $(4, 1)$, $(3, 1)$, $(5, 2)$, $(2, 1)$, $(1, 1)$, and $(0, 1)$. There are two roots of Cartan type: $(1,0)$ and $(5,2)$. 

