# SIL

Sapo Input Language (SIL) is a language designed to provide an input model for Sapo.
It is divided in

- [Header](#header)
- [Symbol definitions](#symdef)
- [Matrices defnition](#matdef)
- [Footer](#footer)

## <a name="header">Header
It consists in

### problem definition (required)

```C++
problem: reachability | synthesis;
```

### number of iterations (required)

```C++
iterations: num;
```

### maximum number of parameter splits (optional)

```C++
max_parameter_splits: num;
```

When parameter synthesis returns an empty set, the result may be due to over-approximation. 
Since the approximation also directly depends on the size of the parameter set, in order to decrease it, 
the initial parameter set can be split in subsets and the computation can be repeated for each them.
The `max_parameter_splits` option declares the maximum number of splits that Sapo can perform 
searching for a non-empty set satisfying the specification.

## <a name="symdef">Symbol definitions

In this section, we define variables with dynamics, parameters, constants and the specification (if required).
Each symbol must be defined before it is used.

### Variables
Variables can be defined with their initial value bounds. Multiple variables can be defined in the same statement. In this case, all variables get the same bounds.

```C++
var v1 in [a,b];
var v2, v3 in [c,d];
```
Otherwise, just the name is provided. Multiple variables can be defined in the same command.

```C++
var v1, v2;
```
Variable names follow the C++ rule for identifiers.

### Parameters

They follow the syntax of variable definition.

```C++
param p1, p2 in [0.01, 0.02];
```
```C++
param p1;
```

### Constants
They define a map between identifiers and numerical values.
The usual rules hold for constant names, and the corresponding expression must be numeric only.

```C++
const c1 = 1.2;
const c2 = 34 * 0.78;
const c3 = c1 + 2*c2;
```

### Definitions
They are used to name subexpressions wich occur in different places. They can contain variables, parameters, constants and other definitions.

```C++
define ex = v + c*p/2;
```

### Dynamics
They represent the evolution of variables with resect to time.

```C++
dynamic(v1) = v1 + p * (v2 - c + v3);
```
where

- `v1` is a name of an already defined variable
- all symbols occurring in the expression must be already defined
- the expression must be linear with regard to parameters
- the expression must be polynomial with respect to the variables

Each variable must have a corresponding dynamic defined.

### Specification
If problem is `synthesis`, a STL specification must be provided. If the problem is `reachability`, this statement is ignored.

+ An atomic formula is of the form
	
	```C++
	ex1 ~ ex2
	```
	where

	- `ex1` and `ex2` are expression involving numbers, variables and constants
	- `~` is one of `<`, `<=`, `>`, `>=` or `=` 

+ Negation is written
	
	```C++
	! f
	```

+ Conjunction is written
	
	```C++
	f1 && f2
	```

+ Disjunction is written
	
	```C++
	f1 || f2
	```

+ Until is written
	
	```C++
	f1 U[a,b] f2
	```
	
	This formula holds at time `t` if there exists a time point `t'` in `[t+a, t+b]` such that `f2` hlds at time `t'` and `f1` hold in all time points from `t` to `t'`.

+ Eventually is written
	
	```C++
	F[a,b] f
	```
	
	This means that there exists a point in `[t+a, t+b]` where `f` holds.

+ Globally is written
	
	```C++
	G[a,b] f
	```
	
	This means that at all time points in '[t+a, t+b]` the formula `f` holds

Finally, the specification is provided as

```C++
spec: STLformula;
```

### Assumption
An assumption is a way to restrict valid states to those satisfying a given formula.
That is, during reachability steps, the states violating the assumptions are discarded.

```C++
assume x + y >= 0;
```

This must be an atomic formula which expression is linear in the variables and contains no parameters.
Multiple assumptions are treated like conjunctions.

Assumptions are currently supported only for `reachability`, not for `synthesis`.

## <a name="matdef">Matrices definition
In this section, we provide the definition of the polytopes used to approximate reachable sets and parameters.

### Directions of initial set
A valid direction is any expression which is linear w.r.t. variables and does not contain parameters.

```C++
direction v1 + 2*v2 - v3/2 in [a, b];
```
It is possible to define named directions

```C++
direction sum: v1 + v2 + v3 in [0,10];
```

It is also possible to define fixed directions

```C++
direction x + y - z = 0;
// this is equivalent to 
// direction x + y - z in [0, 0];
```

To the directions defined in this way we must add the implicit directions defined while declaing variables with bound. For example, the definition
`var x in [a, b];`
is equivalent to 

```C++
var x;
...
direction default_x: x in [a, b];
```

Notice that the implicit direction defined during variable declaration is always named `direction_<var_name>`.

### Template for parallelotpe bundle
In order to represent polytopes, Sapo uses parallelotope bundles. It is possible to define this bundle by explicitly giving a template matrix.

```C++
template = {
	{num_1, ..., num_n},
	{id_1, num_2, ..., id_n},
	.
	.
	.
	{num_1, ..., num_n}
}
```
Each row represent a parallelotope by referencing its directions.
This can be done giving the name, if provided in the definition, or by their number. In fact, each direction has an associated number depending on the order of definition ('0' the first, `1` the second and so on).

The user must pay attention in decribing bounded parallelotopes.

If the template is not provided, then it is automatically generated.
If the template is partial, in that it does not contain all defined directions, a new template is computed that contains all parallelotopes explicitly defined.

### Directions for parameter set
The directions for parameter work in the same way as those for variables

```C++
parameter_direction p1 - p2 in [a, b];

parameter_direction d1: p1 + 4*p3 in [0,1];
```

Again, parameter directions can be implicitly defined in parameter declaration.
For example, `param p in [-1, 1];` is equivalent to

```C++
param p;
...
parameter_direction default_p: p in [-1, 1];
```

Polytopes are not supported for parameter sets, so there can be at most as many directions as the number of parameters (also counting the implicitly defined ones)


## <a name="footer">Footer
In this section, we can provide some options to tune the behavior of Sapo.
All them are optional, so this section can be empty.

We can define
- how we compute the image of a polytope (optional)
	
	```option transformation AFO (default) | OFO;```
	In general, `AFO` gives more accurate results but `OFO` is faster.

- whether or not to perform an optional decomposition phase after computing set images.
	
	```option decomposition;```

- set a weigth `alpha` in `[0,1]` used in decomposition, to define the importance given to polytope volume over maximum side length in the cost function
	
	``` option sapo_alpha num;```

## Comments
SIL understands C/C++ like comments, both sigle and multi-line

```C++
...
// single line comment
...
/*
multi-line
comment
*/
...
``` 

# Examples

### SIR
Here is a SIL input file for a SIR epidemic model

```C++
problem: reachability;

iterations: 10;
max_parameter_splits: 0;

var s in [0.2, 0.3];
var i in [0.001, 0.1];
var r in [0.7, 0.8];

param beta in [0.055, 0.1];
param mu in [0.00001, 0.001];
param gamma in [0.0027, 0.0055];
param alpha in [0.05, 0.07];

dynamic(s) = s - beta*s*i - mu*s + gamma*r;
dynamic(i) = i + beta*s*i - alpha*i;
dynamic(r) = r + mu*s - gamma*r + alpha*i;

spec: r > 0.4;

assume(s <= 1);
assume(i <= 1);
assume(r <= 1);

assume(s >= 0);
assume(i >= 0);
assume(r >= 0);


assume(s + i + r <= 1);
assume(s + i + r >= 0);

option transformation OFO;
```

### Van der Pol oscillator

```C++
problem: reachability;

iterations: 30;


var x in [0, 0.01];
var y in [1.99, 2];

dynamic(x) = x + (y)*0.02;
dynamic(y) = y + (0.5*(1-x*x)*y - x)*0.02;

direction diff: y - x in [-10, 10];
direction sum: x + y in [-10, 10];

template = {
	{default_x, default_y},
	{default_x, diff},
	{default_x, sum},
	{default_y, diff},
	{default_y, sum},
	{diff, sum}
}

```
