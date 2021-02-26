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

### variable modality (optional)

```C++
variable_mode: boxes (default) | parallelotopes | polytopes;
```

### parameter modality (optional)

```C++
parameter_mode: boxes (default) | parallelotopes;
```

### number of iterations (required)

```C++
iterations: num;
```

## <a name="symdef">Symbol definitions

In this section, we define variables with dynamics, parameters, constants and the specification (if required).
Each symbol must be defined before it is used.

### Variables
If the variable modality is `boxes`, they are defined with the initial values

```C++
var v in [a,b];
```
Otherwise, just the name is provided. Multiple variables can be defined in the same command.

```C++
var v1, v2;
```
Variable names follow the C++ rule for identifiers.

### Parameters

They follow the syntax of variable definition. If parameter modality is `boxes`, they require the initial interval.

```C++
param p1, p2 in [0.01, 0.02];
```

Otherwise, no interval must be provided

```C++
param p1;
```

### Constants
They define a map between identifiers and numerical values.
The usual rules hold for constant names, and the corresponding expression must be numeric only.

```C++
const c1 = 1.2;
const c2 = 34 * 0.78;
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
If problem is `synthesis`, a STL specification must be provided.

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

## <a name="matdef">Matrices definition
In this section, we provide the definition of the polytopes used to approximate reachable sets and parameters.

### Variables
If variable modality is `boxes`, no matrix must be defined.
If variable modality is `parallelotopes`, then we must define the direction matrix.
That is, for each direction we write

```C++
direction <num_1, ..., num_n> in [a, b];
```
We have as many elements in a direction as the number of variables.
The meaning is

```C++
a <= num_1*v_1 + ... + num_n * v_n <= b
```
We must provide as many directions as the number ov variables, paying attention to defined a bounded set.

If variable modality is `polytopes`, then we must define a set of directions as described above.
This time we can provide a number of directions which is higher than the number of variables.
Then, we define the *template matrix* in which each row represents a parallelotope.

```C++
template = {
	{num_1, ..., num_n},
	.
	.
	.
	{num_1, ..., num_n}
}
```
where each row has as many elements as the number of variables, and each value corresponds to a direction.
For example, number `0` is the first direction (in order of definition), `1` is the second and so on.
The template matrix can have any number of rows (at least one).

### Parameters
They work similarly to the variables.
If parameter modality is `boxes`, no matrix is required.
Otherwise, we must provide the directions for parameter parallelotope.

```C++
parameter_direction <num_1, ..., num_m> in [a, b];
```
where there is a number for each parameter, and the meaning is

```C++
a <= num_1 * p_1 + ... + num_n * p_m <= b
```

We recall that parameter modality cannot be `polytopes`.

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
