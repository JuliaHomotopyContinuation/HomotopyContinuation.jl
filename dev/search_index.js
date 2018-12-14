var documenterSearchIndex = {"docs": [

{
    "location": "#",
    "page": "Introduction",
    "title": "Introduction",
    "category": "page",
    "text": ""
},

{
    "location": "#Introduction-1",
    "page": "Introduction",
    "title": "Introduction",
    "category": "section",
    "text": "Homotopy.Continuation.jl is a package for solving systems of polynomials equations with only finitely many solutions using numerical homotopy continuation. If this is your first time reading this documentation, we recommend you start with the getting started guide."
},

{
    "location": "#Contents-1",
    "page": "Introduction",
    "title": "Contents",
    "category": "section",
    "text": "Pages = [\"solving.md\", \"systems.md\", \"homotopies.md\", \"predictors-correctors.md\", \"pathtracking.md\", \"reference.md\"]"
},

{
    "location": "solving/#",
    "page": "Solving Polynomial Systems",
    "title": "Solving Polynomial Systems",
    "category": "page",
    "text": ""
},

{
    "location": "solving/#Solving-Polynomial-Systems-1",
    "page": "Solving Polynomial Systems",
    "title": "Solving Polynomial Systems",
    "category": "section",
    "text": "At the heart of the package is the solve function. It takes a bunch of different input combinations and returns an AffineResult or ProjectiveResult depending on the input.The solve function works in 3 stages.It takes the input and constructs a homotopy H(xt) such that H(x1)=G(x) and H(x0)=F(x) as well as start solutions mathcalX where for all x_1  mathcalX we have H(x_1 1)  0. This step highly depends on the concrete input you provide.\nThen all start solutions x(1) = x_1  mathcalX will be tracked to solutions x(t_e) such that H(x(t_e) t_e)  0 using a predictor-corrector scheme where t_e is a value between 0 and 1 (by default 01).\nFrom these intermediate solutions x(t_e) we start the so called endgame. This is an algorithm to predict the value x(0) for each intermediate solution.The reason for step 3 is that the final solutions can be singular which provides significant challenges for our gradient based predictor-corrector methods. For more background also check the FAQ."
},

{
    "location": "solving/#HomotopyContinuation.solve",
    "page": "Solving Polynomial Systems",
    "title": "HomotopyContinuation.solve",
    "category": "function",
    "text": "solve(F; options...)\n\nSolve the system F using a total degree homotopy. F can be\n\nVector{<:MultivariatePolynomials.AbstractPolynomial} (e.g. constructed by @polyvar)\nSystems.AbstractSystem (the system has to represent a homogenous polynomial system.)\n\nExample\n\nAssume we want to solve the system F(xy) = (x^2+y^2+1 2x+3y-1).\n\n@polyvar x y\nsolve([x^2+y^2+1, 2x+3y-1])\n\nIf you polynomial system is already homogenous, but you would like to consider it as an affine system you can do\n\n@polyvar x y z\nsolve([x^2+y^2+z^2, 2x+3y-z], homvar=z)\n\nThis would result in the same result as solve([x^2+y^2+1, 2x+3y-1]).\n\nTo solve F by a custom Systems.AbstractSystem you can do\n\n@polyvar x y z\n# The system `F` has to be homgoenous system\nF = Systems.SPSystem([x^2+y^2+z^2, 2x+3y-z]) # Systems.SPSystem <: Systems.AbstractSystem\n# To solve the original affine system we have to tell that the homogenization variable has index 3\nsolve(F, homvar=3)\n\nor equivalently (in this case) by\n\nsolve([x^2+y^2+z^2, 2x+3y-z], system=Systems.SPSystem)\n\nStart Target Homotopy\n\nsolve(G, F, start_solutions; options...)\n\nSolve the system F by tracking the each provided solution of G (as provided by start_solutions).\n\nExample\n\n@polyvar x y\nG = [x^2-1,y-1]\nF = [x^2+y^2+z^2, 2x+3y-z]\nsolve(G, F, [[1, 1], [-1, 1]])\n\nParameter Homotopy\n\nsolve(F::Vector{<:MultivariatePolynomials.AbstractPolynomial},\n    startsolutions; parameters::Vector{<:MP.AbstractVariable}, p₁, p₀, γ₁=nothing, γ₀=nothing)\n\nSolve the parameter homotopy\n\nH(x t) = F(x (tγ₁p₁+(1-t)γ₀p₀)  (tγ₁+(1-t)γ₀))\n\n, where p₁ and p₀ are a vector of parameter values for F and γ₁ and γ₀ are complex numbers. If γ₁ or γ₀ is nothing, it is assumed that γ₁ and γ₀ are 1. The input parameters specifies the parameter variables of F which should be considered as parameters. Neccessarily, length(parameters) == length(p₁) == length(p₀).\n\nsolve(F::Vector{<:MultivariatePolynomials.AbstractPolynomial},\n        startsolutions; parameters::Vector{<:MP.AbstractVariable},\n        startparameters, targetparameters,\n        startgamma=randn(ComplexF64), targetgamma=randn(ComplexF64))\n\nThis is a non-unicode variant where γ₁=start_parameters, γ₀=target_parameters,     γ₁=start_gamma, γ₀=target_gamma.\n\nExample\n\nWe want to solve a parameter homotopy H(xt) = F(x t1 0+(1-t)2 4) where\n\nF(x a) = (x₁^2-a₁ x₁x₂-a₁+a₂)\n\nand let\'s say we are only intersted in tracking of 11. This can be accomplished as follows\n\n@polyvar x[1:2] a[1:2]\nF = [x[1]^2-a[1], x[1]*x[2]-a[1]+a[2]]\nstartsolutions = [[1, 1]]\nsolve(F, startsolutions, parameters=a, p₁=p₁, p₀=p₀)\n# If you don\'t like unicode this is also possible\nsolve(F, startsolutions, parameters=a, startparameters=p₁, targetparameters=p₀)\n\nAbstract Homotopy\n\nsolve(H::Homotopies.AbstractHomotopy, start_solutions; options...)\n\nSolve the homotopy H by tracking the each solution of H( t) (as provided by start_solutions) from t=1 to t=0. Note that H has to be a homotopy between homogenous polynomial systems. If it should be considered as an affine system indicate which is the index of the homogenization variable, e.g. solve(H, startsolutions, homvar=3) if the third variable is the homogenization variable.\n\nOptions\n\nGeneral options:\n\nsystem::Systems.AbstractSystem: A constructor to assemble a Systems.AbstractSystem. The default is Systems.FPSystem. This constructor is only applied to the input of solve. The constructor is called with system(polynomials, variables) where polynomials is a vector of MultivariatePolynomials.AbstractPolynomials and variables determines the variable ordering.\nhomotopy::Systems.AbstractHomotopy: A constructor to construct a Homotopies.AbstractHomotopy. The default is StraightLineHomotopy. The constructor is called with homotopy(start, target) where start and target are homogenous Systems.AbstractSystems.\nseed::Int: The random seed used during the computations.\nhomvar::Union{Int,MultivariatePolynomials.AbstractVariable}: This considers the homogenous system F as an affine system which was homogenized by homvar. If F is an AbstractSystem homvar is the index (i.e. Int) of the homogenization variable. If F is an AbstractVariables (e.g. created by @polyvar x) homvar is the actual variable used in the system F.\nendgame_start=0.1: The value of t for which the endgame is started.\nreport_progress=true: Whether a progress bar should be printed to STDOUT.\nthreading=true: Enable or disable multi-threading.\n\nPathtracking specific:\n\ncorrector::Correctors.AbstractCorrector: The corrector used during in the predictor-corrector scheme. The default is Correctors.Newton.\ncorrector_maxiters=3: The maximal number of correction steps in a single step.\npredictor::Predictors.AbstractPredictor: The predictor used during in the predictor-corrector scheme. The default is Predictors.RK4.\nrefinement_maxiters=corrector_maxiters: The maximal number of correction steps used to refine the final value.\nrefinement_tol=1e-8: The precision used to refine the final value.\ntol=1e-7: The precision used to track a value.\ninitial_steplength=0.1: The initial step size for the predictor.\nminimal_steplength=1e-14: The minimal step size. If the size of step is below this the path is considered failed.\nmaxiters=1000: The maximal number of steps per path.\n\nEndgame specific options\n\ncauchy_loop_closed_tolerance=1e-3: The tolerance for which is used to determine whether a loop is closed. The distance between endpoints is normalized by the maximal difference between any point in the loop and the starting point.\ncauchy_samples_per_loop=6: The number of samples used to predict an endpoint. A higher number of samples should result in a better approximation. Note that the error should be roughly t^n where t is the current time of the loop and n is cauchy_samples_per_loop.\negtol=1e-10: This is the tolerance necessary to declare the endgame converged.\nmaxnorm=1e5: If our original problem is affine we declare a path at infinity if the infinity norm with respect to the standard patch is larger than maxnorm.\nmaxwindingnumber=15: The maximal windingnumber we try to find using Cauchys integral formula.\nmax_extrapolation_samples=4: During the endgame a Richardson extrapolation is used to improve the accuracy of certain approximations. This is the maximal number of samples used for this.\nminradius=1e-15: A path is declared false if the endgame didn\'t finished until then.\nsampling_factor=0.5: During the endgame we approach 0 by the geometric series h^kR₀ where h is sampling_factor and R₀ the endgame start provided in runendgame.\nmaxiters_per_step=100: The maximal number of steps bewtween two samples.\n\n\n\n\n\n"
},

{
    "location": "solving/#The-*solve*-function-1",
    "page": "Solving Polynomial Systems",
    "title": "The solve function",
    "category": "section",
    "text": "solve"
},

{
    "location": "solving/#HomotopyContinuation.Solving.AffineResult",
    "page": "Solving Polynomial Systems",
    "title": "HomotopyContinuation.Solving.AffineResult",
    "category": "type",
    "text": "AffineResult <: Result\n\nThe result of an (non-homogenous) system of polynomials.\n\n\n\n\n\n"
},

{
    "location": "solving/#HomotopyContinuation.Solving.ProjectiveResult",
    "page": "Solving Polynomial Systems",
    "title": "HomotopyContinuation.Solving.ProjectiveResult",
    "category": "type",
    "text": "ProjectiveResult <: Result\n\nThe result of a homogenous system of polynomials.\n\n\n\n\n\n"
},

{
    "location": "solving/#HomotopyContinuation.Solving.results",
    "page": "Solving Polynomial Systems",
    "title": "HomotopyContinuation.Solving.results",
    "category": "function",
    "text": "results(result; onlyreal=false, realtol=1e-6, onlynonsingular=false, onlysigular=false, singulartol=1e10, onlyfinite=true)\n\nReturn all PathResults for which the given conditions apply.\n\nExample\n\nR = solve(F)\n\n# This gives us all PathResults considered non-singular and real (but still as a complex vector).\nrealsolutions = results(R, onlyreal=true, onlynonsingular=true)\n\n\n\n\n\n"
},

{
    "location": "solving/#HomotopyContinuation.Solving.mapresults",
    "page": "Solving Polynomial Systems",
    "title": "HomotopyContinuation.Solving.mapresults",
    "category": "function",
    "text": "mapresults(f::Function, result; conditions...)\n\nApply the function f to all PathResults for which the given conditions apply. For the possible conditions see results.\n\nExample\n\n# This gives us all solutions considered real (but still as a complex vector).\nrealsolutions = results(solution, R, onlyreal=true)\n\n\n\n\n\n"
},

{
    "location": "solving/#HomotopyContinuation.Solving.solutions",
    "page": "Solving Polynomial Systems",
    "title": "HomotopyContinuation.Solving.solutions",
    "category": "function",
    "text": "solutions(result; conditions...)\n\nReturn all solution (as Vectors) for which the given conditions apply. For the possible conditions see results.\n\nExample\n\njulia> @polyvar x y\njulia> result = solve([(x-2)y, y+x+3]);\njulia> solutions(result)\n[[2.0+0.0im, -5.0+0.0im], [-3.0+0.0im, 0.0+0.0im]]\n\n\n\n\n\n"
},

{
    "location": "solving/#HomotopyContinuation.Solving.realsolutions",
    "page": "Solving Polynomial Systems",
    "title": "HomotopyContinuation.Solving.realsolutions",
    "category": "function",
    "text": "realsolutions(result; tol=1e-6, conditions...)\n\nReturn all real solution (as Vectors of reals) for which the given conditions apply. For the possible conditions see results. Note that onlyreal is always true and realtol is now tol.\n\nExample\n\njulia> @polyvar x y\njulia> result = solve([(x-2)y, y+x+3]);\njulia> realsolutions(result)\n[[2.0, -5.0], [-3.0, 0.0]]\n\n\n\n\n\n"
},

{
    "location": "solving/#HomotopyContinuation.Solving.uniquesolutions",
    "page": "Solving Polynomial Systems",
    "title": "HomotopyContinuation.Solving.uniquesolutions",
    "category": "function",
    "text": "uniquesolutions(R::Result; tol=1e-6, multiplicities=false, conditions...)\n\nReturn all unique solutions. If multiplicities is true, then all unique solutions with their correspnding multiplicities as pairs (s, m) where s is the solution and m the multiplicity are returned. For the possible conditions see results.\n\nExample\n\njulia> @polyvar x;\njulia> uniquesolutions([(x-3)^3*(x+2)], multiplicities=true)\n[([3.0+0.0im], 3), ([-2.0+0.0im], 1)]\njulia> uniquesolutions([(x-3)^3*(x+2)])\n[[3.0+0.0im], [-2.0+0.0im]]\n\n\n\n\n\n"
},

{
    "location": "solving/#HomotopyContinuation.Solving.finite",
    "page": "Solving Polynomial Systems",
    "title": "HomotopyContinuation.Solving.finite",
    "category": "function",
    "text": "finite(result::AffineResults; conditions...)\n\nReturn all PathResults for which the solution is finite. This is just a shorthand for results(R; onlyfinite=true, conditions...). For the possible conditions see results.\n\n\n\n\n\n"
},

{
    "location": "solving/#Base.real-Tuple{Union{Result, Array{#s115,1} where #s115<:PathResult}}",
    "page": "Solving Polynomial Systems",
    "title": "Base.real",
    "category": "method",
    "text": "real(result, tol=1e-6)\n\nGet all results where the solutions are real with the given tolerance tol. See isreal for details regarding the determination of \'realness\'.\n\n\n\n\n\n"
},

{
    "location": "solving/#HomotopyContinuation.Solving.atinfinity",
    "page": "Solving Polynomial Systems",
    "title": "HomotopyContinuation.Solving.atinfinity",
    "category": "function",
    "text": "atinfinity(result::AffineResult)\n\nGet all results where the solutions is at infinity.\n\n\n\n\n\n"
},

{
    "location": "solving/#HomotopyContinuation.Solving.singular",
    "page": "Solving Polynomial Systems",
    "title": "HomotopyContinuation.Solving.singular",
    "category": "function",
    "text": "singular(result::Results; conditions...)\n\nReturn all PathResults for which the solution is singular. This is just a shorthand for results(R; onlysingular=true, conditions...). For the possible conditions see results.\n\n\n\n\n\n"
},

{
    "location": "solving/#HomotopyContinuation.Solving.nonsingular",
    "page": "Solving Polynomial Systems",
    "title": "HomotopyContinuation.Solving.nonsingular",
    "category": "function",
    "text": "nonsingular(result::Results; conditions...)\n\nReturn all PathResults for which the solution is non-singular. This is just a shorthand for results(R; onlynonsingular=true, conditions...). For the possible conditions see results.\n\n\n\n\n\n"
},

{
    "location": "solving/#HomotopyContinuation.Solving.failed",
    "page": "Solving Polynomial Systems",
    "title": "HomotopyContinuation.Solving.failed",
    "category": "function",
    "text": "failed(result)\n\nGet all results where the path tracking failed.\n\n\n\n\n\n"
},

{
    "location": "solving/#HomotopyContinuation.Solving.multiplicities",
    "page": "Solving Polynomial Systems",
    "title": "HomotopyContinuation.Solving.multiplicities",
    "category": "function",
    "text": "multiplicities(vectors, tol, distance)\n\nReturns an array of arrays of integers. Each vector v in vectors contains all indices i,j such that V[i] and V[j] have distance at most tol.\n\n\n\n\n\nmultiplicities(V::Results; tol=1e-6)\n\nReturns a Vector of Vector{PathResult}s grouping the PathResults whose solutions appear with multiplicities greater 1 in \'V\'. Two solutions are regarded as equal, when their pairwise distance is less than \'tol\'.\n\n\n\n\n\n"
},

{
    "location": "solving/#HomotopyContinuation.Solving.seed",
    "page": "Solving Polynomial Systems",
    "title": "HomotopyContinuation.Solving.seed",
    "category": "function",
    "text": "seed(result)\n\nThe random seed used in the computation.\n\n\n\n\n\n"
},

{
    "location": "solving/#HomotopyContinuation.Solving.nresults",
    "page": "Solving Polynomial Systems",
    "title": "HomotopyContinuation.Solving.nresults",
    "category": "function",
    "text": "nresults(result; onlyreal=false, realtol=1e-6, onlynonsingular=false, singulartol=1e10, onlyfinite=true)\n\nThe number of solutions which satisfy the corresponding predicates.\n\nExample\n\nresult = solve(F)\n# Get all non-singular results where all imaginary parts are smaller than 1e-8\nnresults(result, onlyreal=true, realtol=1e-8, onlynonsingular=true)\n\n\n\n\n\n"
},

{
    "location": "solving/#HomotopyContinuation.Solving.nfinite",
    "page": "Solving Polynomial Systems",
    "title": "HomotopyContinuation.Solving.nfinite",
    "category": "function",
    "text": "nfinite(affineresult)\n\nThe number of finite solutions.\n\n\n\n\n\n"
},

{
    "location": "solving/#HomotopyContinuation.Solving.nreal",
    "page": "Solving Polynomial Systems",
    "title": "HomotopyContinuation.Solving.nreal",
    "category": "function",
    "text": "nreal(result; tol=1e-6)\n\nThe number of real solutions where all imaginary parts of each solution are smaller than tol.\n\n\n\n\n\n"
},

{
    "location": "solving/#HomotopyContinuation.Solving.nsingular",
    "page": "Solving Polynomial Systems",
    "title": "HomotopyContinuation.Solving.nsingular",
    "category": "function",
    "text": "nsingular(result; tol=1e10)\n\nThe number of singular solutions. A solution is considered singular if its windingnumber is larger than 1 or the condition number is larger than tol.\n\n\n\n\n\n"
},

{
    "location": "solving/#HomotopyContinuation.Solving.nnonsingular",
    "page": "Solving Polynomial Systems",
    "title": "HomotopyContinuation.Solving.nnonsingular",
    "category": "function",
    "text": "nnonsingular(result; tol=1e-10)\n\nThe number of non-singular solutions.\n\n\n\n\n\n"
},

{
    "location": "solving/#HomotopyContinuation.Solving.natinfinity",
    "page": "Solving Polynomial Systems",
    "title": "HomotopyContinuation.Solving.natinfinity",
    "category": "function",
    "text": "natinfinity(affineresult)\n\nThe number of solutions at infinity.\n\n\n\n\n\n"
},

{
    "location": "solving/#HomotopyContinuation.Solving.nfailed",
    "page": "Solving Polynomial Systems",
    "title": "HomotopyContinuation.Solving.nfailed",
    "category": "function",
    "text": "nafailed(result)\n\nThe number of failed paths.\n\n\n\n\n\n"
},

{
    "location": "solving/#The-result-of-*solve*-1",
    "page": "Solving Polynomial Systems",
    "title": "The result of solve",
    "category": "section",
    "text": "Depending on the input solve returns one of the following typesAffineResult\nProjectiveResultA Result is a wrapper around the results of each single path (PathResult) and it contains some additional informations like the used random seed for the computation.In order to analyze a Result we provide the following helper functionsresults\nmapresults\nsolutions\nrealsolutions\nuniquesolutions\nfinite\nBase.real(::Solving.Results)\natinfinity\nsingular\nnonsingular\nfailed\nmultiplicities\nseed\nIf you are interested in the number of solutions of a certain kind we also provide the following helper functions.nresults\nnfinite\nnreal\nnsingular\nnnonsingular\nnatinfinity\nnfailed"
},

{
    "location": "solving/#HomotopyContinuation.Solving.PathResult",
    "page": "Solving Polynomial Systems",
    "title": "HomotopyContinuation.Solving.PathResult",
    "category": "type",
    "text": "PathResult(startvalue, pathtracker_result, endgamer_result, solver)\n\nA PathResult is the result of the tracking of a path (inclusive endgame). Its fields are\n\nreturncode: One of :success, :at_infinity or any error code from the EndgamerResult\nsolution::Vector{T}: The solution vector. If the algorithm computed in projective space and the solution is at infinity then the projective solution is given. Otherwise an affine solution is given if the startvalue was affine and a projective solution is given if the startvalue was projective.\nresidual::Float64: The value of the infinity norm of H(solution, 0).\nnewton_residual: The value of the 2-norm of J_H(textsolution)^-1H(textsolution 0)\ncondition_number: This is the condition number of the Jacobian at the solution. A high condition number indicates a singularity.\nwindingnumber: The estimated winding number\nangle_to_infinity: The angle to infinity is the angle of the solution to the hyperplane where the homogenizing coordinate is 0.\nreal_solution: Indicates whether the solution is real given the defined tolerance at_infinity_tol (from the solver options).\nstartvalue: The startvalue of the path\niterations: The number of iterations the pathtracker needed.\nendgame_iterations: The number of steps in the geometric series the endgamer did.\nnpredictions: The number of predictions the endgamer did.\npredictions: The predictions of the endgamer.\n\n\n\n\n\n"
},

{
    "location": "solving/#HomotopyContinuation.Solving.solution",
    "page": "Solving Polynomial Systems",
    "title": "HomotopyContinuation.Solving.solution",
    "category": "function",
    "text": "solution(pathresult)\n\nGet the solution of the path.\n\n\n\n\n\n"
},

{
    "location": "solving/#HomotopyContinuation.Solving.residual",
    "page": "Solving Polynomial Systems",
    "title": "HomotopyContinuation.Solving.residual",
    "category": "function",
    "text": "residual(pathresult)\n\nGet the residual of the solution x of the path, i.e., H(x t)_infty.\n\n\n\n\n\n"
},

{
    "location": "solving/#HomotopyContinuation.Solving.startsolution",
    "page": "Solving Polynomial Systems",
    "title": "HomotopyContinuation.Solving.startsolution",
    "category": "function",
    "text": "startsolution(pathresult)\n\nGet the start solution of the solution x of the path.\n\n\n\n\n\n"
},

{
    "location": "solving/#Base.isreal-Tuple{PathResult}",
    "page": "Solving Polynomial Systems",
    "title": "Base.isreal",
    "category": "method",
    "text": "isreal(pathresult; tol=1e-6)\nisreal(pathresult, tol)\n\nDetermine whether infinity norm of the imaginary part of the so\n\n\n\n\n\n"
},

{
    "location": "solving/#LinearAlgebra.issuccess-Tuple{PathResult}",
    "page": "Solving Polynomial Systems",
    "title": "LinearAlgebra.issuccess",
    "category": "method",
    "text": "issuccess(pathresult)\n\nChecks whether the path is successfull.\n\n\n\n\n\n"
},

{
    "location": "solving/#HomotopyContinuation.Solving.isfailed",
    "page": "Solving Polynomial Systems",
    "title": "HomotopyContinuation.Solving.isfailed",
    "category": "function",
    "text": "isfailed(pathresult)\n\nChecks whether the path failed.\n\n\n\n\n\n"
},

{
    "location": "solving/#HomotopyContinuation.Solving.isaffine",
    "page": "Solving Polynomial Systems",
    "title": "HomotopyContinuation.Solving.isaffine",
    "category": "function",
    "text": "isaffine(pathresult)\n\nChecks whether the path result is affine.\n\n\n\n\n\n"
},

{
    "location": "solving/#HomotopyContinuation.Solving.isprojective",
    "page": "Solving Polynomial Systems",
    "title": "HomotopyContinuation.Solving.isprojective",
    "category": "function",
    "text": "isprojective(pathresult)\n\nChecks whether the path result is projective.\n\n\n\n\n\n"
},

{
    "location": "solving/#HomotopyContinuation.Solving.isatinfinity",
    "page": "Solving Polynomial Systems",
    "title": "HomotopyContinuation.Solving.isatinfinity",
    "category": "function",
    "text": "isatinfinity(pathresult)\n\nChecks whether the path goes to infinity.\n\n\n\n\n\n"
},

{
    "location": "solving/#HomotopyContinuation.Solving.issingular",
    "page": "Solving Polynomial Systems",
    "title": "HomotopyContinuation.Solving.issingular",
    "category": "function",
    "text": "issingular(pathresult; tol=1e10)\nissingular(pathresult, tol)\n\nChecks whether the path result is singular. This is true if the winding number > 1 or if the condition number of the Jacobian is larger than tol.\n\n\n\n\n\n"
},

{
    "location": "solving/#HomotopyContinuation.Solving.isnonsingular",
    "page": "Solving Polynomial Systems",
    "title": "HomotopyContinuation.Solving.isnonsingular",
    "category": "function",
    "text": "isnonsingular(pathresult; tol=1e10)\n\nChecks whether the path result is non-singular. This is true if it is not singular.\n\n\n\n\n\n"
},

{
    "location": "solving/#PathResult-1",
    "page": "Solving Polynomial Systems",
    "title": "PathResult",
    "category": "section",
    "text": "For each path we return a PathResult containing the detailed information about the single path.PathResultThe following helper functions are providedsolution\nresidual\nstartsolution\nBase.isreal(::PathResult)\nLinearAlgebra.issuccess(::PathResult)\nisfailed\nisaffine\nisprojective\nisatinfinity\nissingular\nisnonsingular"
},

{
    "location": "solving/#HomotopyContinuation.Monodromy.monodromy_solve",
    "page": "Solving Polynomial Systems",
    "title": "HomotopyContinuation.Monodromy.monodromy_solve",
    "category": "function",
    "text": "monodromy_solve(F, sols, p; parameters=..., options..., pathtrackerkwargs...)\n\nSolve a polynomial system F(x;p) with specified parameters and initial solutions sols by monodromy techniques. This makes loops in the parameter space of F to find new solutions.\n\nOptions\n\ntarget_solutions_count=nothing: The computations are stopped if this number of solutions is reached.\ndone_callback=always_false: A callback to end the computation early. This function takes 2 arguments. The first one is the new solution x and the second one are all current solutions (including x). Return true if the compuation is done.\nmaximal_number_of_iterations_without_progress::Int=10: The maximal number of iterations (i.e. loops generated) without any progress.\ngroup_action=nothing: A function taking one solution and returning other solutions if there is a constructive way to obtain them, e.g. by symmetry.\nstrategy: The strategy used to create loops. By default this will be Triangle with weights if F is a real system.\nshowprogress=true: Enable a progress meter.\ntol::Float64=1e-5: The tolerance with which it is decided whether two solutions are identical. The paths are tracked with tolerance tol*1e-2.\ngroup_actions=GroupActions(group_action): If there is more than one group action you can use this to chain the application of them.\ngroup_action_on_all_nodes=false: By default the group_action(s) are only applied on the solutions with the main parameter p. If this is enabled then it is applied for every parameter q.\nparameter_sampler=independent_normal: A function taking the parameter p and returning a new random parameter q. By default each entry of the parameter vector is drawn independently from the unviraite normal distribution.\ntimeout=float(typemax(Int)): The maximal number of seconds the computation is allowed to run.\nminimal_number_of_solutions: The minimal number of solutions before a stopping heuristic is applied. By default this is half of target_solutions_count if applicable otherwise 2.\n\n\n\n\n\n"
},

{
    "location": "solving/#HomotopyContinuation.Monodromy.GroupActions",
    "page": "Solving Polynomial Systems",
    "title": "HomotopyContinuation.Monodromy.GroupActions",
    "category": "type",
    "text": "GroupActions(actions::Function...)\n\nStore a bunch of group actions (f1, f2, f3, ...). Each action has to return a tuple. The actions are applied in the following sense\n\nf1 is applied on the original solution s\nf2 is applied on s and the results of 1\nf3 is applied on s and the results of 1) and 2)\n\nand so on\n\nExample\n\njulia> f1(s) = (s * s,);\n\njulia> f2(s) = (2s, -s, 5s);\n\njulia> f3(s) = (s + 1,);\n\njulia> GroupActions(f1)(3)\n(9,)\n\njulia> GroupActions(f1,f2)(3)\n(9, 18, -9, 45)\n\njulia> GroupActions(f1,f2, f3)(3)\n(9, 18, -9, 45, 10, 19, -8, 46)\n\n\n\n\n\n"
},

{
    "location": "solving/#Monodromy-Solve-1",
    "page": "Solving Polynomial Systems",
    "title": "Monodromy Solve",
    "category": "section",
    "text": "monodromy_solve\nMonodromy.GroupActions"
},

{
    "location": "solving/#HomotopyContinuation.Monodromy.Triangle",
    "page": "Solving Polynomial Systems",
    "title": "HomotopyContinuation.Monodromy.Triangle",
    "category": "type",
    "text": "Triangle(;weights=true)\n\nA triangle is a loop consisting of the main node and two addtional nodes. If weights is true the edges are equipped with additional random weights. Note that this is usually only necessary for real parameters.\n\n\n\n\n\n"
},

{
    "location": "solving/#HomotopyContinuation.Monodromy.Petal",
    "page": "Solving Polynomial Systems",
    "title": "HomotopyContinuation.Monodromy.Petal",
    "category": "type",
    "text": "Petal()\n\nA petal is a loop consisting of the main node and one other node connected by two edges with different random weights.\n\n\n\n\n\n"
},

{
    "location": "solving/#Strategies-1",
    "page": "Solving Polynomial Systems",
    "title": "Strategies",
    "category": "section",
    "text": "Monodromy.Triangle\nMonodromy.Petal"
},

{
    "location": "systems/#",
    "page": "Systems",
    "title": "Systems",
    "category": "page",
    "text": ""
},

{
    "location": "systems/#Polynomial-systems-1",
    "page": "Systems",
    "title": "Polynomial systems",
    "category": "section",
    "text": "Polynomial systems can be represented in numerous ways in a computer and each representation has certain tradeoffs. For our purposes the most important thing is that it is fast to evaluate the system. Therefore we automatically convert an input given by DynamicPolynomials to another representation more suitable for numerically evaluations. The default is currently FPSystem."
},

{
    "location": "systems/#HomotopyContinuation.Systems.FPSystem",
    "page": "Systems",
    "title": "HomotopyContinuation.Systems.FPSystem",
    "category": "type",
    "text": "FPSystem(polynomials, vars) <: AbstractSystem\n\nCreate a polynomial system using the FixedPolynomials package.\n\n\n\n\n\n"
},

{
    "location": "systems/#HomotopyContinuation.Systems.SPSystem",
    "page": "Systems",
    "title": "HomotopyContinuation.Systems.SPSystem",
    "category": "type",
    "text": "SPSystem(polynomials, vars) <: AbstractSystem\n\nCreate a system using the StaticPolynomials package. Note that StaticPolynomials leverages Julias metaprogramming capabilities to automatically generate functions to evaluate the system and its Jacobian. These generated functions are very fast but at the cost of possibly large compile times. The compile time depends on the size of the support of the polynomial system. If you intend to solve a large system or you need to solve a system with the same support but different coefficients even large compile times can be worthwile. As a general rule of thumb this usually is twice as fast as solving the same system using FPSystem.\n\nExample\n\nYou can use SPSystem as follows with solve\n\n@polyvar x y\nF = [x^2+3y^4-2, 2y^2+3x*y+4]\nsolve(F, system=SPSystem)\n\n\n\n\n\n"
},

{
    "location": "systems/#HomotopyContinuation.Systems.FixedHomotopy",
    "page": "Systems",
    "title": "HomotopyContinuation.Systems.FixedHomotopy",
    "category": "type",
    "text": "FixedHomotopy(H, t) <: AbstractSystem\n\nFix a homotopy H(x,t) at t\n\n\n\n\n\n"
},

{
    "location": "systems/#Default-systems-1",
    "page": "Systems",
    "title": "Default systems",
    "category": "section",
    "text": "We provide the following systems by default.FPSystem\nSPSystem\nSystems.FixedHomotopy"
},

{
    "location": "systems/#Interface-for-custom-systems-1",
    "page": "Systems",
    "title": "Interface for custom systems",
    "category": "section",
    "text": "The great thing is that you are not limited to the systems provided by default. Maybe your polynomial system has a particular structure which you want to use to efficiently evaluate it. For this you can define your own homotopy by defining a struct with super type Systems.AbstractSystem. For this the following interface has to be defined."
},

{
    "location": "systems/#HomotopyContinuation.SystemsBase.AbstractSystem",
    "page": "Systems",
    "title": "HomotopyContinuation.SystemsBase.AbstractSystem",
    "category": "type",
    "text": "AbstractSystem\n\nRepresenting a system of polynomials.\n\n\n\n\n\n"
},

{
    "location": "systems/#HomotopyContinuation.SystemsBase.AbstractSystemCache",
    "page": "Systems",
    "title": "HomotopyContinuation.SystemsBase.AbstractSystemCache",
    "category": "type",
    "text": "AbstractSystemCache\n\nA cache to avoid allocations for the evaluation of an AbstractSystem.\n\n\n\n\n\n"
},

{
    "location": "systems/#HomotopyContinuation.SystemsBase.NullCache",
    "page": "Systems",
    "title": "HomotopyContinuation.SystemsBase.NullCache",
    "category": "type",
    "text": "NullCache\n\nAn empty cache if no cache is necessary.\n\n\n\n\n\n"
},

{
    "location": "systems/#Types-1",
    "page": "Systems",
    "title": "Types",
    "category": "section",
    "text": "Systems.AbstractSystem\nSystems.AbstractSystemCache\nSystems.NullCache"
},

{
    "location": "systems/#HomotopyContinuation.SystemsBase.cache",
    "page": "Systems",
    "title": "HomotopyContinuation.SystemsBase.cache",
    "category": "function",
    "text": "cache(F::AbstractSystem, x)::AbstractSystemCache\n\nCreate a cache for the evaluation (incl. Jacobian) of F with elements of the type of x. The default implementation returns a NullCache.\n\n\n\n\n\n"
},

{
    "location": "systems/#HomotopyContinuation.SystemsBase.evaluate!",
    "page": "Systems",
    "title": "HomotopyContinuation.SystemsBase.evaluate!",
    "category": "function",
    "text": "evaluate!(u, F::AbstractSystem, x , cache::AbstractSystemCache)\n\nEvaluate the system F at x and store the result in u.\n\n\n\n\n\n"
},

{
    "location": "systems/#HomotopyContinuation.SystemsBase.evaluate",
    "page": "Systems",
    "title": "HomotopyContinuation.SystemsBase.evaluate",
    "category": "function",
    "text": "evaluate(F::AbstractSystem, x::AbstractVector, cache=cache(F, x))\n\nEvaluate the system F at x.\n\n\n\n\n\n"
},

{
    "location": "systems/#HomotopyContinuation.SystemsBase.jacobian!",
    "page": "Systems",
    "title": "HomotopyContinuation.SystemsBase.jacobian!",
    "category": "function",
    "text": "jacobian!(u, F::AbstractSystem, x , cache::AbstractSystemCache)\n\nEvaluate the Jacobian of the system F at x and store the result in u.\n\n\n\n\n\n"
},

{
    "location": "systems/#HomotopyContinuation.SystemsBase.jacobian",
    "page": "Systems",
    "title": "HomotopyContinuation.SystemsBase.jacobian",
    "category": "function",
    "text": "jacobian(F::AbstractSystem, x, cache=cache(F, x))\n\nEvaluate the Jacobian of the system F at x.\n\n\n\n\n\n"
},

{
    "location": "systems/#Base.size-Tuple{HomotopyContinuation.SystemsBase.AbstractSystem}",
    "page": "Systems",
    "title": "Base.size",
    "category": "method",
    "text": "Base.size(F::AbstractSystem)\n\nReturns a tuple (m, n) indicating that F is a system of m polynomials m in n variables.\n\n\n\n\n\n"
},

{
    "location": "systems/#Mandatory-1",
    "page": "Systems",
    "title": "Mandatory",
    "category": "section",
    "text": "The following methods are mandatory to implement.Systems.cache\nSystems.evaluate!\nSystems.evaluate\nSystems.jacobian!\nSystems.jacobian\nBase.size(::Systems.AbstractSystem)"
},

{
    "location": "systems/#HomotopyContinuation.SystemsBase.evaluate_and_jacobian!",
    "page": "Systems",
    "title": "HomotopyContinuation.SystemsBase.evaluate_and_jacobian!",
    "category": "function",
    "text": "evaluate_and_jacobian!(u, U, F, x , cache::AbstractSystemCache)\n\nEvaluate the system F and its Jacobian at x and store the results in u (evalution) and U (Jacobian).\n\n\n\n\n\n"
},

{
    "location": "systems/#HomotopyContinuation.SystemsBase.evaluate_and_jacobian",
    "page": "Systems",
    "title": "HomotopyContinuation.SystemsBase.evaluate_and_jacobian",
    "category": "function",
    "text": "evaluate_and_jacobian(F::AbstractSystem, x , cache=cache(F, x))\n\nEvaluate the system F and its Jacobian at x.\n\n\n\n\n\n"
},

{
    "location": "systems/#Optional-1",
    "page": "Systems",
    "title": "Optional",
    "category": "section",
    "text": "The following methods are mandatory to implement. The following are optional to implement but usually you want to define at least Systems.cache.Systems.evaluate_and_jacobian!\nSystems.evaluate_and_jacobian"
},

{
    "location": "homotopies/#",
    "page": "Homotopies",
    "title": "Homotopies",
    "category": "page",
    "text": ""
},

{
    "location": "homotopies/#Homotopies-1",
    "page": "Homotopies",
    "title": "Homotopies",
    "category": "section",
    "text": "A homotopy is a functionH mathbbC^N  mathbbC  mathbbC^n (xt)  H(xt)where H( t) is a polynomial system for all tmathbbC."
},

{
    "location": "homotopies/#HomotopyContinuation.Homotopies.StraightLineHomotopy",
    "page": "Homotopies",
    "title": "HomotopyContinuation.Homotopies.StraightLineHomotopy",
    "category": "type",
    "text": "StraightLineHomotopy(G, F; gamma=exp(i * 2π*rand()))\n\nConstruct the homotopy H(x t) = γtG(x) + (1-t)F(x).\n\n\n\n\n\n"
},

{
    "location": "homotopies/#HomotopyContinuation.Homotopies.FixedPointHomotopy",
    "page": "Homotopies",
    "title": "HomotopyContinuation.Homotopies.FixedPointHomotopy",
    "category": "type",
    "text": "FixedPointHomotopy(F, x₀; gamma=exp(i * 2π*rand()))\n\nConstruct the homotopy H(x t) = (1-t)F(x) + γt(x-x₀).\n\n\n\n\n\n"
},

{
    "location": "homotopies/#HomotopyContinuation.Homotopies.ParameterHomotopy",
    "page": "Homotopies",
    "title": "HomotopyContinuation.Homotopies.ParameterHomotopy",
    "category": "type",
    "text": "ParameterHomotopy(F, parameters;\n    variables=setdiff(MP.variables(F), parameters),\n    p₁=randn(ComplexF64, length(parameters)),\n    p₀=randn(ComplexF64, length(parameters)),\n    γ₁=nothing, γ₀=nothing)\n\nConstruct the homotopy\n\nH(x t) = F(x (tγ₁p₁+(1-t)γ₀p₀)  (tγ₁+(1-t)γ₀))\n\n, where p₁ and p₀ are a vector of parameter values for F and γ₁ and γ₀ are complex numbers. If γ₁ or γ₀ is nothing, it is assumed that γ₁ and γ₀ are 1. The input parameters specifies the parameter variables of F. Neccessarily, length(parameters) == length(p₁) == length(p₀).\n\nNote that p₁ and p₀ are stored as a tuple p of SVectors and γ₁ and γ₀ are stored as a tuple γ or as γ=nothing\n\nParameterHomotopy(F, parameters;\n    variables=setdiff(MP.variables(F), parameters),\n    startparameters=randn(ComplexF64, length(parameters)),\n    targetparameters=randn(ComplexF64, length(parameters)),\n    startgamma=nothing, targetgamma=nothing)\n\nThis is a non-unicode variant where γ₁=startparameters, γ₀=targetparameters, γ₁=startgamma, γ₀=targetgamma.\n\n\n\n\n\n"
},

{
    "location": "homotopies/#HomotopyContinuation.Homotopies.PatchedHomotopy",
    "page": "Homotopies",
    "title": "HomotopyContinuation.Homotopies.PatchedHomotopy",
    "category": "type",
    "text": "PatchedHomotopy(H::AbstractHomotopy, patch, v::AbstractProjectiveVector)\n\nAugment the homotopy H with the given patch v. This results in the system [H(x,t); v ⋅ x - 1]\n\n\n\n\n\n"
},

{
    "location": "homotopies/#Default-homotopies-1",
    "page": "Homotopies",
    "title": "Default homotopies",
    "category": "section",
    "text": "The following homotopies are available by defaultStraightLineHomotopy\nFixedPointHomotopy\nParameterHomotopyWe also provide more specialised homotopies, which are mostly used internally currently but could be useful in conjunction with the PathTracking.PathTracker primitive.Homotopies.PatchedHomotopy"
},

{
    "location": "homotopies/#Interface-for-custom-homotopies-1",
    "page": "Homotopies",
    "title": "Interface for custom homotopies",
    "category": "section",
    "text": "The great thing is that you are not limited to the homotopies provided by default. You can define your own homotopy by defining a struct with super type Homotopies.AbstractHomotopy. For this the following interface has to be defined."
},

{
    "location": "homotopies/#HomotopyContinuation.HomotopiesBase.AbstractHomotopy",
    "page": "Homotopies",
    "title": "HomotopyContinuation.HomotopiesBase.AbstractHomotopy",
    "category": "type",
    "text": "AbstractHomotopy\n\nRepresenting a homotopy.\n\n\n\n\n\n"
},

{
    "location": "homotopies/#HomotopyContinuation.HomotopiesBase.AbstractHomotopyCache",
    "page": "Homotopies",
    "title": "HomotopyContinuation.HomotopiesBase.AbstractHomotopyCache",
    "category": "type",
    "text": "AbstractHomotopyCache\n\nA cache to avoid allocations for the evaluation of an AbstractHomotopy.\n\n\n\n\n\n"
},

{
    "location": "homotopies/#HomotopyContinuation.HomotopiesBase.NullCache",
    "page": "Homotopies",
    "title": "HomotopyContinuation.HomotopiesBase.NullCache",
    "category": "type",
    "text": "NullCache\n\nThe default AbstractHomotopyCache containing nothing.\n\n\n\n\n\n"
},

{
    "location": "homotopies/#Types-1",
    "page": "Homotopies",
    "title": "Types",
    "category": "section",
    "text": "Homotopies.AbstractHomotopy\nHomotopies.AbstractHomotopyCache\nHomotopies.NullCache"
},

{
    "location": "homotopies/#HomotopyContinuation.HomotopiesBase.cache",
    "page": "Homotopies",
    "title": "HomotopyContinuation.HomotopiesBase.cache",
    "category": "function",
    "text": "cache(H::AbstractHomotopy, x, t)::AbstractHomotopyCache\n\nCreate a cache for the evaluation (incl. Jacobian) of F with elements of the type of x. The default implementation returns HomotopiesBase.NullCache.\n\n\n\n\n\n"
},

{
    "location": "homotopies/#HomotopyContinuation.HomotopiesBase.evaluate!",
    "page": "Homotopies",
    "title": "HomotopyContinuation.HomotopiesBase.evaluate!",
    "category": "function",
    "text": "evaluate!(u, H::AbstractHomotopy, x, t, cache::AbstractHomotopyCache)\n\nEvaluate the homotopy H at (x, t) and store the result in u.\n\n\n\n\n\n"
},

{
    "location": "homotopies/#HomotopyContinuation.HomotopiesBase.jacobian!",
    "page": "Homotopies",
    "title": "HomotopyContinuation.HomotopiesBase.jacobian!",
    "category": "function",
    "text": "jacobian!(u, H::AbstractHomotopy, x, t, cache::AbstractHomotopyCache)\n\nEvaluate the Jacobian of the homotopy H at (x, t) and store the result in u.\n\n\n\n\n\n"
},

{
    "location": "homotopies/#HomotopyContinuation.HomotopiesBase.dt!",
    "page": "Homotopies",
    "title": "HomotopyContinuation.HomotopiesBase.dt!",
    "category": "function",
    "text": "dt!(u, H::AbstractHomotopy, x, t, cache::AbstractHomotopyCache)\n\nEvaluate the homotopy H at (x, t) and store the result in u.\n\n\n\n\n\n"
},

{
    "location": "homotopies/#Base.size-Tuple{HomotopyContinuation.HomotopiesBase.AbstractHomotopy}",
    "page": "Homotopies",
    "title": "Base.size",
    "category": "method",
    "text": "Base.size(H::AbstractHomotopy)\n\nReturns a tuple (m, n) indicating that H is a homotopy of m polynomials m in n variables.\n\n\n\n\n\n"
},

{
    "location": "homotopies/#Mandatory-1",
    "page": "Homotopies",
    "title": "Mandatory",
    "category": "section",
    "text": "The following methods are mandatory to implement.Homotopies.cache\nHomotopies.evaluate!\nHomotopies.jacobian!\nHomotopies.dt!\nBase.size(::Homotopies.AbstractHomotopy)"
},

{
    "location": "homotopies/#HomotopyContinuation.HomotopiesBase.evaluate_and_jacobian!",
    "page": "Homotopies",
    "title": "HomotopyContinuation.HomotopiesBase.evaluate_and_jacobian!",
    "category": "function",
    "text": "evaluate_and_jacobian!(u, U, F, x, t, cache::AbstractHomotopyCache)\n\nEvaluate the homotopy H and its Jacobian at (x, t) and store the results in u (evalution) and U (Jacobian).\n\n\n\n\n\n"
},

{
    "location": "homotopies/#HomotopyContinuation.HomotopiesBase.evaluate_and_jacobian",
    "page": "Homotopies",
    "title": "HomotopyContinuation.HomotopiesBase.evaluate_and_jacobian",
    "category": "function",
    "text": "evaluate_and_jacobian(H::AbstractHomotopy, x, t, cache::AbstractHomotopyCache)\n\nEvaluate the homotopy H and its Jacobian at (x, t).\n\n\n\n\n\n"
},

{
    "location": "homotopies/#HomotopyContinuation.HomotopiesBase.jacobian_and_dt!",
    "page": "Homotopies",
    "title": "HomotopyContinuation.HomotopiesBase.jacobian_and_dt!",
    "category": "function",
    "text": "jacobian_and_dt!(U, u, H, x, t, cache::AbstractHomotopyCache)\n\nEvaluate the homotopy H and its derivative w.r.t. t at (x, t) and store the results in u (evalution) and v (∂t).\n\n\n\n\n\n"
},

{
    "location": "homotopies/#HomotopyContinuation.HomotopiesBase.evaluate",
    "page": "Homotopies",
    "title": "HomotopyContinuation.HomotopiesBase.evaluate",
    "category": "function",
    "text": "evaluate(H::AbstractHomotopy, x, t, cache::AbstractHomotopyCache)\n\nEvaluate the homotopy H at (x, t).\n\n\n\n\n\n"
},

{
    "location": "homotopies/#HomotopyContinuation.HomotopiesBase.jacobian",
    "page": "Homotopies",
    "title": "HomotopyContinuation.HomotopiesBase.jacobian",
    "category": "function",
    "text": "jacobian(H::AbstractHomotopy, x, t, cache::AbstractHomotopyCache)\n\nEvaluate the Jacobian of the homotopy H at (x, t).\n\n\n\n\n\n"
},

{
    "location": "homotopies/#HomotopyContinuation.HomotopiesBase.dt",
    "page": "Homotopies",
    "title": "HomotopyContinuation.HomotopiesBase.dt",
    "category": "function",
    "text": "dt(H::AbstractHomotopy, x::AbstractVector, cache::AbstractHomotopyCache)\n\nEvaluate the homotopy H at (x, t).\n\n\n\n\n\n"
},

{
    "location": "homotopies/#HomotopyContinuation.HomotopiesBase.basehomotopy",
    "page": "Homotopies",
    "title": "HomotopyContinuation.HomotopiesBase.basehomotopy",
    "category": "function",
    "text": "basehomotopy(H::AbstractHomotopy)\n\nReturns the \'proper\' homotopy describing the problem. Any wrapper homotopy recursively calls wrappedhomotopy with the wrapped homotopy as argument.\n\n\n\n\n\n"
},

{
    "location": "homotopies/#Optional-1",
    "page": "Homotopies",
    "title": "Optional",
    "category": "section",
    "text": "Homotopies.evaluate_and_jacobian!\nHomotopies.evaluate_and_jacobian\nHomotopies.jacobian_and_dt!\nHomotopies.evaluate\nHomotopies.jacobian\nHomotopies.dt\nHomotopies.basehomotopy"
},

{
    "location": "predictors-correctors/#",
    "page": "Predictors and Correctors",
    "title": "Predictors and Correctors",
    "category": "page",
    "text": ""
},

{
    "location": "predictors-correctors/#Predictors-and-Correctors-1",
    "page": "Predictors and Correctors",
    "title": "Predictors and Correctors",
    "category": "section",
    "text": "We use a predictor-corrector scheme to track paths. These are the predictors and correctors currently available."
},

{
    "location": "predictors-correctors/#HomotopyContinuation.Predictors.Euler",
    "page": "Predictors and Correctors",
    "title": "HomotopyContinuation.Predictors.Euler",
    "category": "type",
    "text": "Euler()\n\nThis uses the explicit Euler method for prediction, also known as the tangent predictor.\n\n\n\n\n\n"
},

{
    "location": "predictors-correctors/#HomotopyContinuation.Predictors.Heun",
    "page": "Predictors and Correctors",
    "title": "HomotopyContinuation.Predictors.Heun",
    "category": "type",
    "text": "Heun()\n\nThe Heun predictor of order 2.\n\n\n\n\n\n"
},

{
    "location": "predictors-correctors/#HomotopyContinuation.Predictors.Ralston",
    "page": "Predictors and Correctors",
    "title": "HomotopyContinuation.Predictors.Ralston",
    "category": "type",
    "text": "Ralston()\n\nThe Ralston predictor of order 2.\n\n\n\n\n\n"
},

{
    "location": "predictors-correctors/#HomotopyContinuation.Predictors.RK3",
    "page": "Predictors and Correctors",
    "title": "HomotopyContinuation.Predictors.RK3",
    "category": "type",
    "text": "RK3()\n\nThe classical Runge-Kutta predictor of order 3.\n\n\n\n\n\n"
},

{
    "location": "predictors-correctors/#HomotopyContinuation.Predictors.RK4",
    "page": "Predictors and Correctors",
    "title": "HomotopyContinuation.Predictors.RK4",
    "category": "type",
    "text": "RK4()\n\nThe classical Runge-Kutta predictor of order 4.\n\n\n\n\n\n"
},

{
    "location": "predictors-correctors/#HomotopyContinuation.Predictors.Pade21",
    "page": "Predictors and Correctors",
    "title": "HomotopyContinuation.Predictors.Pade21",
    "category": "type",
    "text": "Pade21()\n\nThis uses a Padé-approximation of type (2,1) for prediction.\n\n\n\n\n\n"
},

{
    "location": "predictors-correctors/#HomotopyContinuation.Predictors.NullPredictor",
    "page": "Predictors and Correctors",
    "title": "HomotopyContinuation.Predictors.NullPredictor",
    "category": "type",
    "text": "NullPredictor()\n\nA predictor which does no prediction step, i.e., it just returns the input as its prediction.\n\n\n\n\n\n"
},

{
    "location": "predictors-correctors/#Predictors-1",
    "page": "Predictors and Correctors",
    "title": "Predictors",
    "category": "section",
    "text": "The following predictors are currently implemented.Predictors.Euler\nPredictors.Heun\nPredictors.Ralston\nPredictors.RK3\nPredictors.RK4\nPredictors.Pade21\nPredictors.NullPredictor"
},

{
    "location": "predictors-correctors/#HomotopyContinuation.Correctors.Newton",
    "page": "Predictors and Correctors",
    "title": "HomotopyContinuation.Correctors.Newton",
    "category": "type",
    "text": "Newton(;simplified_last_step=true)\n\nAn ordinary Newton\'s method. If simplified_last_step is true, then for the last iteration the previously Jacobian will be used. This uses an LU-factorization for square systems and a QR-factorization for overdetermined.\n\n\n\n\n\n"
},

{
    "location": "predictors-correctors/#Correctors-1",
    "page": "Predictors and Correctors",
    "title": "Correctors",
    "category": "section",
    "text": "Correctors.Newton"
},

{
    "location": "pathtracking/#",
    "page": "Path tracker",
    "title": "Path tracker",
    "category": "page",
    "text": ""
},

{
    "location": "pathtracking/#HomotopyContinuation.PathTracking.pathtracker_startsolutions",
    "page": "Path tracker",
    "title": "HomotopyContinuation.PathTracking.pathtracker_startsolutions",
    "category": "function",
    "text": "pathtracker_startsolutions(args...; kwargs...)\n\nConstruct a PathTracking.PathTracker and startsolutions in the same way solve does it. This also takes the same input arguments as solve. This is convenient if you want to investigate single paths.\n\n\n\n\n\n"
},

{
    "location": "pathtracking/#Path-tracking-1",
    "page": "Path tracker",
    "title": "Path tracking",
    "category": "section",
    "text": "We also export a path tracking primitive to make the core path tracking routine available for other applications. At the heart is a PathTracking.PathTracker object which holds all the state. The easiest way to construct a PathTracker is to use the pathtracker_startsolutions routine.pathtracker_startsolutions"
},

{
    "location": "pathtracking/#HomotopyContinuation.PathTracking.PathTracker",
    "page": "Path tracker",
    "title": "HomotopyContinuation.PathTracking.PathTracker",
    "category": "type",
    "text": " PathTracker(H::Homotopies.AbstractHomotopy, x₁, t₁, t₀; options...)::PathTracker\n\nCreate a PathTracker to track x₁ from t₁ to t₀. The homotopy H needs to be homogenous. Note that a PathTracker is also a (mutable) iterator.\n\nOptions\n\ncorrector::Correctors.AbstractCorrector: The corrector used during in the predictor-corrector scheme. The default is Correctors.Newton.\ncorrector_maxiters=3: The maximal number of correction steps in a single step.\ninitial_steplength=0.1: The step length of the first step.\nmaxiters=10_000: The maximal number of iterations the path tracker has available.\nminimal_steplength=1e-14: The minimal step length.\npredictor::Predictors.AbstractPredictor: The predictor used during in the predictor-corrector scheme. The default is [Predictors.RK4](@ref)()`.\nrefinement_maxiters=corrector_maxiters: The maximal number of correction steps used to refine the final value.\nrefinement_tol=1e-8: The precision used to refine the final value.\ntol=1e-7: The precision used to track a value.\n\n\n\n\n\n"
},

{
    "location": "pathtracking/#HomotopyContinuation.PathTracking.PathTrackerResult",
    "page": "Path tracker",
    "title": "HomotopyContinuation.PathTracking.PathTrackerResult",
    "category": "type",
    "text": " PathTrackerResult(tracker)\n\nContaining the result of a tracked path. The fields are\n\nsuccessfull::Bool Indicating whether tracking was successfull.\nreturncode::PathTracking.Status.t If the tracking was successfull then it is PathTracking.Status.success.\nx::V The result.\nt::Float64 The t when the path tracker stopped.\n\n\n\n\n\n"
},

{
    "location": "pathtracking/#Types-1",
    "page": "Path tracker",
    "title": "Types",
    "category": "section",
    "text": "PathTracking.PathTracker\nPathTracking.PathTrackerResult"
},

{
    "location": "pathtracking/#HomotopyContinuation.PathTracking.track!",
    "page": "Path tracker",
    "title": "HomotopyContinuation.PathTracking.track!",
    "category": "function",
    "text": " track!(tracker, x₁, t₁, t₀; setup_patch=true, checkstartvalue=true, compute_ẋ=true)\n\nTrack a value x₁ from t₁ to t₀ using the given PathTracker tracker. Returns one of the enum values of PathTracking.Status.t indicating the status. If the tracking was successfull it is PathTracking.Status.success. If setup_patch is true then AffinePatches.setup! is called at the beginning of the tracking.\n\ntrack!(x₀, tracker, x₁, t₁, t₀; options...)\n\nAdditionally also stores the result in x₀ if the tracking was successfull.\n\n\n\n\n\n"
},

{
    "location": "pathtracking/#HomotopyContinuation.PathTracking.track",
    "page": "Path tracker",
    "title": "HomotopyContinuation.PathTracking.track",
    "category": "function",
    "text": "track(tracker, x₁, t₁, t₀; options...)::PathTrackerResult\n\nTrack a value x₁ from t₁ to t₀ using the given PathTracker tracker. This returns a PathTrackerResult. This modifies tracker. See track! for the possible options.\n\n\n\n\n\n"
},

{
    "location": "pathtracking/#HomotopyContinuation.PathTracking.setup!",
    "page": "Path tracker",
    "title": "HomotopyContinuation.PathTracking.setup!",
    "category": "function",
    "text": "setup!(pathtracker, x₁, t₁, t₀, setup_patch=true, checkstartvalue=true, compute_ẋ=true)\n\nSetup pathtracker to track x₁ from t₁ to t₀. Use this if you want to use the pathtracker as an iterator.\n\n\n\n\n\n"
},

{
    "location": "pathtracking/#HomotopyContinuation.PathTracking.currx",
    "page": "Path tracker",
    "title": "HomotopyContinuation.PathTracking.currx",
    "category": "function",
    "text": "currx(tracker::PathTracker)\n\nReturn the current value of x.\n\n\n\n\n\n"
},

{
    "location": "pathtracking/#HomotopyContinuation.PathTracking.currt",
    "page": "Path tracker",
    "title": "HomotopyContinuation.PathTracking.currt",
    "category": "function",
    "text": " currt(tracker::PathTracker)\n\nCurrent t.\n\n\n\n\n\n"
},

{
    "location": "pathtracking/#HomotopyContinuation.PathTracking.currΔt",
    "page": "Path tracker",
    "title": "HomotopyContinuation.PathTracking.currΔt",
    "category": "function",
    "text": " currΔt(tracker::PathTracker)\n\nCurrent steplength Δt.\n\n\n\n\n\n"
},

{
    "location": "pathtracking/#HomotopyContinuation.PathTracking.curriters",
    "page": "Path tracker",
    "title": "HomotopyContinuation.PathTracking.curriters",
    "category": "function",
    "text": " curriters(tracker::PathTracker)\n\nCurrent number of iterations.\n\n\n\n\n\n"
},

{
    "location": "pathtracking/#HomotopyContinuation.PathTracking.currstatus",
    "page": "Path tracker",
    "title": "HomotopyContinuation.PathTracking.currstatus",
    "category": "function",
    "text": " currstatus(tracker::PathTracker)\n\nCurrent status.\n\n\n\n\n\n"
},

{
    "location": "pathtracking/#HomotopyContinuation.PathTracking.tol",
    "page": "Path tracker",
    "title": "HomotopyContinuation.PathTracking.tol",
    "category": "function",
    "text": " tol(tracker::PathTracker)\n\nCurrent tolerance.\n\n\n\n\n\n"
},

{
    "location": "pathtracking/#HomotopyContinuation.PathTracking.corrector_maxiters",
    "page": "Path tracker",
    "title": "HomotopyContinuation.PathTracking.corrector_maxiters",
    "category": "function",
    "text": " corrector_maxiters(tracker::PathTracker)\n\nCurrent correction maxiters.\n\n\n\n\n\n"
},

{
    "location": "pathtracking/#HomotopyContinuation.PathTracking.refinement_tol",
    "page": "Path tracker",
    "title": "HomotopyContinuation.PathTracking.refinement_tol",
    "category": "function",
    "text": " refinement_tol(tracker::PathTracker)\n\nCurrent refinement tolerance.\n\n\n\n\n\n"
},

{
    "location": "pathtracking/#HomotopyContinuation.PathTracking.refinement_maxiters",
    "page": "Path tracker",
    "title": "HomotopyContinuation.PathTracking.refinement_maxiters",
    "category": "function",
    "text": " refinement_maxiters(tracker::PathTracker)\n\nCurrent refinement maxiters.\n\n\n\n\n\n"
},

{
    "location": "pathtracking/#HomotopyContinuation.PathTracking.set_tol!",
    "page": "Path tracker",
    "title": "HomotopyContinuation.PathTracking.set_tol!",
    "category": "function",
    "text": " set_tol!(tracker::PathTracker, tol)\n\nSet the current tolerance to tol.\n\n\n\n\n\n"
},

{
    "location": "pathtracking/#HomotopyContinuation.PathTracking.set_corrector_maxiters!",
    "page": "Path tracker",
    "title": "HomotopyContinuation.PathTracking.set_corrector_maxiters!",
    "category": "function",
    "text": " set_corrector_maxiters!(tracker::PathTracker, n)\n\nSet the current correction maxiters to n.\n\n\n\n\n\n"
},

{
    "location": "pathtracking/#HomotopyContinuation.PathTracking.set_refinement_tol!",
    "page": "Path tracker",
    "title": "HomotopyContinuation.PathTracking.set_refinement_tol!",
    "category": "function",
    "text": " set_refinement_maxiters!(tracker::PathTracker, tol)\n\nSet the current refinement tolerance to tol.\n\n\n\n\n\n"
},

{
    "location": "pathtracking/#HomotopyContinuation.PathTracking.set_refinement_maxiters!",
    "page": "Path tracker",
    "title": "HomotopyContinuation.PathTracking.set_refinement_maxiters!",
    "category": "function",
    "text": " set_refinement_maxiters!(tracker::PathTracker, n)\n\nSet the current refinement maxiters to n.\n\n\n\n\n\n"
},

{
    "location": "pathtracking/#Methods-1",
    "page": "Path tracker",
    "title": "Methods",
    "category": "section",
    "text": "To track from a start to an endpoint with the PathTracker we provide the following routines.PathTracking.track!\nPathTracking.track\nPathTracking.setup!To introspect the current state and change settings we provide the following routines.PathTracking.currx\nPathTracking.currt\nPathTracking.currΔt\nPathTracking.curriters\nPathTracking.currstatus\nPathTracking.tol\nPathTracking.corrector_maxiters\nPathTracking.refinement_tol\nPathTracking.refinement_maxiters\nPathTracking.set_tol!\nPathTracking.set_corrector_maxiters!\nPathTracking.set_refinement_tol!\nPathTracking.set_refinement_maxiters!"
},

{
    "location": "reference/#",
    "page": "Reference",
    "title": "Reference",
    "category": "page",
    "text": ""
},

{
    "location": "reference/#Reference-1",
    "page": "Reference",
    "title": "Reference",
    "category": "section",
    "text": ""
},

{
    "location": "reference/#Input-1",
    "page": "Reference",
    "title": "Input",
    "category": "section",
    "text": "We support any polynomials which follow the MultivariatePolynomials interface. By default we export the routines @polyvar, PolyVar, differentiate and variables from the DynamicPolynomials implementation. With these you can simply create variables# Create variables x, y, z\n@polyvar x y z\nf = x^2+y^2+z^2\n\n# You can also create an array of variables\n@polyvar x[1:3] # This creates x1, x2, x3 accessed by x[1], x[2], x[3]\nf = dot(x, x) # = x[1]^2+x[2]^2+x[3]^2\n\n# Also you can create matrices of variables\n# This creates x1_1, x1_2, x2_1, x2_2 accessed by\n# x[1,1], x[1,2], x[2,1], x[2,2]\n@polyvar x[1:2, 1:2]"
},

{
    "location": "reference/#HomotopyContinuation.Utilities.ishomogenous",
    "page": "Reference",
    "title": "HomotopyContinuation.Utilities.ishomogenous",
    "category": "function",
    "text": "ishomogenous(f::MP.AbstractPolynomialLike)\n\nChecks whether f is homogenous.\n\nishomogenous(polys::Vector{MP.AbstractPolynomialLike})\n\nChecks whether each polynomial in polys is homogenous.\n\n\n\n\n\nishomogenous(f::MP.AbstractPolynomialLike, v::Vector{<:MP.AbstractVariable})\n\nChecks whether f is homogenous in the variables v.\n\nishomogenous(polys::Vector{<:MP.AbstractPolynomialLike}, v::Vector{<:MP.AbstractVariable})\n\nChecks whether each polynomial in polys is homogenous in the variables v.\n\n\n\n\n\n"
},

{
    "location": "reference/#HomotopyContinuation.Utilities.uniquevar",
    "page": "Reference",
    "title": "HomotopyContinuation.Utilities.uniquevar",
    "category": "function",
    "text": "uniquevar(f::MP.AbstractPolynomialLike, tag=:x0)\nuniquevar(F::Vector{<:MP.AbstractPolynomialLike}, tag=:x0)\n\nCreates a unique variable.\n\n\n\n\n\n"
},

{
    "location": "reference/#HomotopyContinuation.Utilities.homogenize",
    "page": "Reference",
    "title": "HomotopyContinuation.Utilities.homogenize",
    "category": "function",
    "text": "homogenize(f::MP.AbstractPolynomial, variable=uniquevar(f))\n\nHomogenize the polynomial f by using the given variable variable.\n\nhomogenize(F::Vector{<:MP.AbstractPolynomial}, variable=uniquevar(F))\n\nHomogenize each polynomial in F by using the given variable variable.\n\n\n\n\n\nhomogenize(f::MP.AbstractPolynomial, v::Vector{<:MP.AbstractVariable}, variable=uniquevar(f))\n\nHomogenize the variables v in the polynomial f by using the given variable variable.\n\nhomogenize(F::Vector{<:MP.AbstractPolynomial}, v::Vector{<:MP.AbstractVariable}, variable=uniquevar(F))\n\nHomogenize the variables v in each polynomial in F by using the given variable variable.\n\n\n\n\n\n"
},

{
    "location": "reference/#Utilities-1",
    "page": "Reference",
    "title": "Utilities",
    "category": "section",
    "text": "ishomogenous\nuniquevar\nhomogenize"
},

{
    "location": "reference/#HomotopyContinuation.AffinePatches.OrthogonalPatch",
    "page": "Reference",
    "title": "HomotopyContinuation.AffinePatches.OrthogonalPatch",
    "category": "type",
    "text": "OrthogonalPatch()\n\n\n\n\n\n"
},

{
    "location": "reference/#HomotopyContinuation.AffinePatches.EmbeddingPatch",
    "page": "Reference",
    "title": "HomotopyContinuation.AffinePatches.EmbeddingPatch",
    "category": "type",
    "text": "EmbeddingPatch()\n\nHolds an AbstractProjectiveVector onto its affine patch. With this the effect is basically the same as tracking in affine space.\n\n\n\n\n\n"
},

{
    "location": "reference/#HomotopyContinuation.AffinePatches.RandomPatch",
    "page": "Reference",
    "title": "HomotopyContinuation.AffinePatches.RandomPatch",
    "category": "type",
    "text": "RandomPatch()\n\nA random patch. The vector has norm 1.\n\n\n\n\n\n"
},

{
    "location": "reference/#HomotopyContinuation.AffinePatches.FixedPatch",
    "page": "Reference",
    "title": "HomotopyContinuation.AffinePatches.FixedPatch",
    "category": "type",
    "text": "FixedPatch()\n\n\n\n\n\n"
},

{
    "location": "reference/#AffinePatches-1",
    "page": "Reference",
    "title": "AffinePatches",
    "category": "section",
    "text": "Affine patches are there to augment projective system such that they can be considered as (locally) affine system. By default the following patches are definedAffinePatches.OrthogonalPatch\nAffinePatches.EmbeddingPatch\nAffinePatches.RandomPatch\nAffinePatches.FixedPatch"
},

{
    "location": "reference/#HomotopyContinuation.AffinePatches.AbstractAffinePatch",
    "page": "Reference",
    "title": "HomotopyContinuation.AffinePatches.AbstractAffinePatch",
    "category": "type",
    "text": "AbstractAffinePatch\n\nAn affine patch is a hyperplane defined by vx-1=0.\n\n\n\n\n\n"
},

{
    "location": "reference/#HomotopyContinuation.AffinePatches.state",
    "page": "Reference",
    "title": "HomotopyContinuation.AffinePatches.state",
    "category": "function",
    "text": "state(::AbstractAffinePatch, x)::AbstractAffinePatchState\n\nConstruct the state of the path from x.\n\n\n\n\n\n"
},

{
    "location": "reference/#HomotopyContinuation.AffinePatches.onpatch!",
    "page": "Reference",
    "title": "HomotopyContinuation.AffinePatches.onpatch!",
    "category": "function",
    "text": "onpatch!(x::AbstractVector, ::AbstractAffinePatchState)\n\nScale a vector x such that it is on the affine patch.\n\n\n\n\n\n"
},

{
    "location": "reference/#HomotopyContinuation.AffinePatches.setup!",
    "page": "Reference",
    "title": "HomotopyContinuation.AffinePatches.setup!",
    "category": "function",
    "text": "setup!(::AbstractAffinePatchState, x::AbstractVector)\n\nSetup the affine patch depending on x and modify x if necessary. This is only called once at the beginning of a tracked path.\n\n\n\n\n\n"
},

{
    "location": "reference/#HomotopyContinuation.AffinePatches.changepatch!",
    "page": "Reference",
    "title": "HomotopyContinuation.AffinePatches.changepatch!",
    "category": "function",
    "text": "changepatch!(::AbstractAffinePatch, x::AbstractVector)\n\nThe same as setup! but only called during the path tracking.\n\n\n\n\n\n"
},

{
    "location": "reference/#Interface-1",
    "page": "Reference",
    "title": "Interface",
    "category": "section",
    "text": "Each patch has to follow the following interface.AffinePatches.AbstractAffinePatch\nAffinePatches.state\nAffinePatches.onpatch!\nAffinePatches.setup!\nAffinePatches.changepatch!"
},

]}
