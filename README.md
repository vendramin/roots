# Computing finite Weyl groupoids

The repository contains the [GAP](https://www.gap-system.org) code and the database used in the paper

[Computing finite Weyl groupoids](http://arxiv.org/abs/2505.24555)  
by I. Angiono, L. Vendramin  

[![DOI](https://zenodo.org/badge/827310499.svg)](https://doi.org/10.5281/zenodo.15554003)

Folders
=======
The following folders contain all the exceptional braidings of diagonal type (that are neither of Cartan type nor of super type):

```
rank2/
rank3/ 
rank4/
rank5/
rank6/
rank7/
rank8/
```

The following folders contain some exceptional 
contragradient data of finite-dimensional Lie superalgebras:

```
super_char3/
super_char5/
```

Files
=====

roots.g 
-------
This is the main script. It determines whether a root system of a given
braiding of diagonal type is finite. In this case, the script computes the set
of positive roots. 

### Example:
In this example, we compute the set of positive roots
of the root system of type $\frak{ufo}(7)a$. 

```
gap> Read("roots.g");
gap> roots("rank2/ufo7a.g");
#I  not logging
File: rank2/ufo7a.g
The braiding is:
[ [     -E(4),  -E(12)^7 ],
  [         1,        -1 ] ]
Longest word:
[ 1, 2, 1, 2, 1 ]
Positive roots (5):
[ [ 1, 0 ], [ 3, 1 ], [ 2, 1 ], [ 1, 1 ], [ 0, 1 ] ]
Cartan roots (0):
[  ]
```
Similarly, the script can be used with any of the files in the folder mentioned
above, or with any other braiding matrix (in which case, you need to create a
file containing the corresponding braiding matrix).

relations.g 
-----------
Once the set of positive roots is computed, the script computes the Cartan
roots, the dimension and the defining relations of the corresponding Nichols
algebra.

### Example:
In this example, we compute the defining relations of the Nichols algebra of type $\frak{ufo}(7)a$:
```
gap> Read("relations.g");
gap> relations("rank2/ufo7a.g");
#I  not logging
The braiding is:
[ [     -E(4),  -E(12)^7 ],
  [         1,        -1 ] ]
Dimension: 144
Longest word:
[ 1, 2, 1, 2, 1 ]
Positive roots (5):
[ [ 1, 0 ], [ 3, 1 ], [ 2, 1 ], [ 1, 1 ], [ 0, 1 ] ]
Cartan roots (0):
[  ]
Heights: [ 4, 2, 3, 3, 2 ]
Relations:
x_{1}^{4}
x_{2}^{2}
[x_{1\,1\,2},x_{1\,2}]_c
[x_{1},x_{3\alpha_{1}+2\alpha_{2}}]_c-(0)x_{1\,1\,2}^2
```
lyndon.g
-------
Once the set of positive roots is computed, the script computes
the set of Lyndon words and the hyperwords. 

### Example: 
We now compute the Lyndon words of the Nichols algebra of type $\frak{ufo}(7)a$.
```
gap> Read("lyndon.g");
gap> lyndon("rank2/ufo7a.g");
#I  not logging
File: rank2/ufo7a.g
The braiding is:
[ [     -E(4),  -E(12)^7 ],
  [         1,        -1 ] ]
Longest word:
[ 1, 2, 1, 2, 1 ]
Positive roots (5):
[ [ 1, 0 ], [ 3, 1 ], [ 2, 1 ], [ 1, 1 ], [ 0, 1 ] ]
Cartan roots (0):
[  ]
rec( decompositions := [ 0, 0, [ [ [ 1 ], [ 2 ] ] ], [ [ [ 1 ], [ 1, 2 ] ] ],
      [ [ [ 1 ], [ 1, 1, 2 ] ] ] ], ordering := [ 5, 1, 4, 3, 2 ],
  words := [ [ 2 ], [ 1 ], [ 1, 2 ], [ 1, 1, 2 ], [ 1, 1, 1, 2 ] ] )
```
### Example: 
Finally, we compute the hyperwords of the Nichols algebra of type $\frak{ufo}(7)a$:
```
gap> hyperwords("rank2/ufo7a.g");
#I  not logging
File: rank2/ufo7a.g
The braiding is:
[ [     -E(4),  -E(12)^7 ],
  [         1,        -1 ] ]
Longest word:
[ 1, 2, 1, 2, 1 ]
Positive roots (5):
[ [ 1, 0 ], [ 3, 1 ], [ 2, 1 ], [ 1, 1 ], [ 0, 1 ] ]
Cartan roots (0):
[  ]
x_{\alpha_{2}}=x_{2}
x_{\alpha_{1}}=x_{1}
x_{[ 1, 1 ]}=[x_{[ 1, 0 ]},x_{[ 0, 1 ]}]_c
x_{[ 2, 1 ]}=[x_{[ 1, 0 ]},x_{[ 1, 1 ]}]_c
x_{[ 3, 1 ]}=[x_{[ 1, 0 ]},x_{[ 2, 1 ]}]_c
```
super.g
-------
This script determines whether a root system of a given congragradient Lie
superalgebra is finite. In this case, the script computes the set of positive
roots, the odd roots and the super dimension. 

### Example: 
Here we compute the superdimension corresponding to the Lie superalgebra $\frak{g}(2,3)$:
```
gap> Read("super.g");
gap> super("super_char3/g23.g");
The parity is: [ -1, -1, -1 ]
The normalized matrix is:
 . 1 .
 1 2 1
 . 1 .
Longest word:
[ 1, 2, 1, 2, 3, 2, 1, 3, 2, 1 ]
Positive roots (10):
[ [ 1, 0, 0 ], [ 1, 1, 0 ], [ 1, 2, 0 ], [ 0, 1, 0 ], [ 1, 3, 1 ], [ 1, 2, 1 ], [ 0, 2, 1 ], [ 1, 1, 1 ],
  [ 0, 1, 1 ], [ 0, 0, 1 ] ]
Odd positive roots (7):
[ [ 1, 0, 0 ], [ 1, 2, 0 ], [ 0, 1, 0 ], [ 1, 3, 1 ], [ 0, 2, 1 ], [ 1, 1, 1 ], [ 0, 0, 1 ] ]
Odd nondegenerate (1):
[ [ 0, 1, 0 ] ]
The set Nabla+ is:
[ [ 0, 0, 1 ], [ 0, 1, 0 ], [ 0, 1, 1 ], [ 0, 2, 0 ], [ 0, 2, 1 ], [ 1, 0, 0 ], [ 1, 1, 0 ], [ 1, 1, 1 ],
  [ 1, 2, 0 ], [ 1, 2, 1 ], [ 1, 3, 1 ] ]
superdimension of the contragradient Lie superalgebra g(B,p): (12|14)
```
### Example:
We now compute the superdimension corresponding to the Lie superalgebra $\frak{br}(2,5)$: 
```
gap> Read("super.g");
gap> super("super_char5/brown25.g");
The parity is: [ 1, -1 ]
The normalized matrix is:
 2 2
 1 .
Longest word:
[ 1, 2, 1, 2, 1, 2, 1, 2 ]
Positive roots (8):
[ [ 1, 0 ], [ 3, 1 ], [ 2, 1 ], [ 5, 3 ], [ 3, 2 ], [ 4, 3 ], [ 1, 1 ], [ 0, 1 ] ]
Odd positive roots (6):
[ [ 3, 1 ], [ 2, 1 ], [ 5, 3 ], [ 4, 3 ], [ 1, 1 ], [ 0, 1 ] ]
Odd nondegenerate (2):
[ [ 2, 1 ], [ 1, 1 ] ]
The set Nabla+ is:
[ [ 0, 1 ], [ 1, 0 ], [ 1, 1 ], [ 2, 1 ], [ 2, 2 ], [ 3, 1 ], [ 3, 2 ], [ 4, 2 ], [ 4, 3 ], [ 5, 3 ] ]
superdimension of the contragradient Lie superalgebra g(B,p): (10|12)
```

