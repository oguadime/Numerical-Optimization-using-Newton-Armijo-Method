# Numerical-Optimization-using-Newton-Armijo-Method
This repository contains MATLAB implementations of numerical methods for optimizing scalar functions using the Newton-Armijo method. It includes support for both analytical and numerical differentiation, along with visualization tools for convergence analysis.

# Decsription 
This repository includes MATLAB scripts that implement the Newton-Armijo optimization method. Key highlights include:
- Newton-Armijo Optimization (newtsol_opt.m): Implements Newton's method with an Armijo line search for minimizing scalar functions. Supports both analytical and numerical derivatives.
- Test Function (f1.m): Defines a test function with an added perturbation controlled by a global parameter alpha, along with its analytical first and second derivatives.
- Driver Script: Runs the optimization for different initial values and derivative settings while generating performance plots.

# Key Features
- Newton-Armijo Method: Uses Newton’s method for optimization, includes a backtracking line search (Armijo condition) to maintain stability, Supports finite difference approximations for first and second derivatives.
- Derivative Handling: Analytical derivatives (if available), Numerical derivatives using forward differences.
- Visualization and Performance Tracking: Stores iteration history including function values and gradient norms, Generates convergence plots.

# Usage
Each script is modular and can be executed directly in MATLAB. Key parameters (e.g., tolerances, step sizes) can be adjusted in the driver script.
Example 1:  Running the Optimization with Different Initial Values
```matlab
tola = 1e-6;
tolr = tola;
initv = [0.5, 2, 2];
global alpha;
alpha = 0; 
[x, hist] = newtsol_opt(initv(1), 'f1', 1e-4, tola, tolr, 1, 1);
```
Example 2: Using Finite Difference Derivatives
```matlab
[x, hist] = newtsol_opt(2, 'f1', 1e-8, 1e-6, 1e-6, 1, 1);
```
Example 3: Running with Analytical Derivatives
```matlab
[x, hist] = newtsol_opt(0.5, 'f1', 0, 1e-6, 1e-6, 0, 0);
```
# License
This project is licensed under the MIT License - see the LICENSE file for details.
```
Feel free to adjust any part of this README to better fit your specific needs or preferences.
```



