
# Ex1 – Polynomial Operations (Introduction to Computer Science, Ariel University)

This project implements a full set of static polynomial-processing utilities
as part of the *Introduction to Computer Science (Ex1)* assignment, Ariel University.

The implementation includes evaluating polynomials, finding roots, computing intersections,
parsing from strings, generating polynomials from sample points, computing arc length, area
between curves, and drawing polynomials using the provided GUI.


The assignment includes the following polynomial utilities:

### Core Operations
- `f(poly, x)` – Evaluate polynomial at x  
- `add(p1, p2)` – Add two polynomials  
- `mul(p1, p2)` – Multiply two polynomials  
- `derivative(p)` – Compute derivative polynomial  
- `equals(p1, p2)` – Compare polynomials up to EPS  

### Numeric Algorithms
- `root_rec(...)` – Recursive bisection root-finding  
- `sameValue(p1, p2, x1, x2, eps)` – Find intersection point of two polynomials  
- `length(p, x1, x2, segments)` – Approximate arc length  
- `area(p1, p2, x1, x2, trapezoids)` – Compute area between polynomials  

### Parsing & Construction
- `poly(p)` – Convert polynomial to human-readable string  
- `getPolynomFromString(s)` – Parse polynomial from string  
- `PolynomFromPoints(xx, yy)` – Build polynomial from 2–3 known points


  The following image shows the expected output of running `Ex1_GUI`,  
including the two polynomials and the shaded area between them
<img width="783" height="832" alt="צילום מסך 2025-11-26 152322" src="https://github.com/user-attachments/assets/6f191d9f-79c9-4ecc-a93e-031c54fa573e" />


Oded Magnaji
ID: 212197479
Ex1 – Introduction to Computer Science
Ariel University, 2026

