# Thompsons-group

This repository contains a method for computing the normal form of an element in Thompson's group F. 

Definition: Thompson's group F is the group generated by infinitely many generators x_0, x_1, x_2, ... satisfying the relations<br/> 
x_k^{-1} x_n x_k = x_{n+1} for all k<n. 

Every element g of F can be expressed uniquely in the its normal form, i.e. <br/>
g = x_0^{a_0} ... x_n^{a_n} x_n^{-b_n} ... x_n^{-b_0} , <br/>
where a_0, ..., a_n, b_0, ..., b_n are non-negative integers, exactly one of a_n, b_n is non-zero, and: <br/>
(a_i and b_i not equal 0)    implies   (a_{i+1} not equal 0 or b_{i+1} not equal 0) ,<br/>
for all i.

The normal form can be computed via an algorithm from the article "Thompson's group and
public key cryptography" by Shpilrain and Ushakov, availble at
http://www.sci.ccny.cuny.edu/~shpil/thomcryp.pdf

In Thompsonsgroup.py, this algorithm is implemented using Python. 

## Example usage

We represent words in Thompson's group as a list of 2-lists, each 2-list consisting of an index (in {0,1,2,...}) and a number in {1,-1} (the exponent). So, the word x_1^2 x_3 x_2^{-1} will be represented as
```python 
[ [1,1], [1,1], [1,3], [2,-1] ]
```
To compute the normal form of this element, enter
```python 
normalForm([ [1,1], [1,1], [1,3], [2,-1] ])
```
