The purpose of this repo is to generate all isomorphic homogeneous polynomials of degree n over a finite field F_q. 

Because of the scientific nature of this, this code isn't meant to be generally reused.

The code gets slow very fast when we either increase n or q. We ran this code for degree 6 field 2 in +-3 hours.
And degree 5 field 3 in the same amount of time. Increasing either q or n will be probably be either uncalculable or won't fit in memory (probably both).

This code was written to count the amount of smooth homogeneous polynomials of degree n in the P^2 space over a finite field F_q. The repo can be found [here](https://github.com/Chrisvossetje/smooth_polynomial_counter)

Developed by,  
Jacco Hijmans and Chris Vos  
Assisted by,  
Prof. Dr. C.F. Faber at Utrecht University