PlaneCuts

The objective of this project is to apply the Avis-Fukuda reverse pivot algorithm to practical computational geometry problems, for instance computing a polytope from its definition as intersection of half spaces.

The main component is a linear programming simplex solver using the pivot method, with reverse pivot exploration.

Input is a linear programming problem, not necessarily in standard or canonical form, possibly with free variables. See SimplexSolver() constructor for details; see LPSolverTest for examples. After solve(), get_solution() gives the LPP solution, then reverse_solve() outputs the list of valid vertices.

The geometric application is not yet implemented - watch this space!

JMA
