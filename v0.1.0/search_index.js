var documenterSearchIndex = {"docs": [

{
    "location": "index.html#",
    "page": "Introduction",
    "title": "Introduction",
    "category": "page",
    "text": ""
},

{
    "location": "index.html#Introduction-1",
    "page": "Introduction",
    "title": "Introduction",
    "category": "section",
    "text": "HomotopyContinuation.jl is a package for solving square polynomial systems via homotopy continuation.The aim of this project is twofold: establishing a fast numerical polynomial solver in Julia and at the same time providing a highly customizable algorithmic environment for researchers for designing and performing individual experiments.You can simply install this package via the Julia package managerPkg.add(\"HomotopyContinuation\");"
},

{
    "location": "index.html#A-first-example-1",
    "page": "Introduction",
    "title": "A first example",
    "category": "section",
    "text": "HomotopyContinuation.jl aims at having easy-to-understand top-level commands. For instance, suppose we wanted to solve the following systemf= beginbmatrix x^2+y y^2-1endbmatrix  First, we have to define f in Julia. For this purpose HomotopyContinuation.jl provides an interface to MultivariatePolynomials.jl for human-readable and easy constructable input. In the following we will use the DynamicPolynomials.jl implementation of that interface.import DynamicPolynomials: @polyvar # @polyvar is a function for initializing variables.\n\n@polyvar x y # initialize the variables x y\nf = [x^2+y, y^2-1]To solve  f=0 we execute the following command.using HomotopyContinuation # load the module HomotopyContinuation\n\nsolve(f) # solves for f=0(see here for a list of options that can be passed to solve).The last command will return a type HomotopyContinuation.Result{Complex{Float64}} of length 4 (one entry for each solution):julia> ans\n\njulia> HomotopyContinuation.Result{Complex{Float64}}\n# paths → 4\n# successfull paths → 4\n# solutions at infinity → 0\n# singular solutions → 0\n# real solutions → 2\nHomotopyContinuation.PathResult{Complex{Float64}}[4]Let us see what is the information that we get. Four paths were attempted to be solved, four of which were completed successfully. Since we tried to solve an affine system, the algorithm checks whether there are solutions at infinity: in this case there are none. None of the solutions is singular and two of them are real. To access the first solution in the array we writejulia> ans[1]\n\njulia> HomotopyContinuation.PathResult{Complex{Float64}}\nreturncode → :success\nsolution → Complex{Float64}[2]\nsingular → false\nresidual → 1.02e-15…\nnewton_residual → 8.95e-16…\nlog10_condition_number → 0.133…\nwindingnumber → 1\nangle_to_infinity → 0.615…\nreal_solution → true\nstartvalue → Complex{Float64}[2]\niterations → 17\nendgame_iterations → 5\nnpredictions → 2\npredictions → Vector{Complex{Float64}}[2]The returncode tells us that the pathtracking was successfull. What do the other entries of the table tell us? Let us consider the most relevant  (for a complete list of explanations consider this section).solution: the zero that is computed (here it is -1-1).\nsingular: boolean that shows whether the zero is singular.\nresidual: the computed value of f(-1-1).\nangle_to_infinity: the algorithms homogenizes the system f and then computes all solutions in projective space. The angle to infinity is the angle of the solution to the hyperplane where the homogenizing coordinate is 0.\nreal_solution: boolean that shows whether the zero is real.Suppose we were only interested in the real solutions. The command to extract them issolutions(solve(f), only_real=true)(a detailed explanation of the solutions function is here). Indeed, we havejulia> [ans[i].solution for i=1:2]\njulia> Vector{Complex{Float64}}[2]\nComplex{Float64}[2]\n1.00… - 2.66e-15…im\n-1.00… + 1.33e-15…im\nComplex{Float64}[2]\n-1.00… + 2.72e-15…im\n-1.00… + 1.44e-15…imwhich are the two real zeros of f. By assigning the boolean values in the solutions function we can filter the solutions given by solve(f) according to our needs.We solve some more elaborate systems in the example section.JuliaHomotopyContinuation also supports input of type BigFloat."
},

{
    "location": "examples.html#",
    "page": "Examples",
    "title": "Examples",
    "category": "page",
    "text": ""
},

{
    "location": "examples.html#examples-1",
    "page": "Examples",
    "title": "Examples",
    "category": "section",
    "text": ""
},

{
    "location": "examples.html#Computing-the-degree-of-a-variety-1",
    "page": "Examples",
    "title": "Computing the degree of a variety",
    "category": "section",
    "text": "Consider the projective variety in the 2-dimensional complex projective space CP^2.V =  x^2 + y^2 - z^2 = 0 The degree of V is the number of intersection points of V with a generic line.   Let us see what it is. First we initialize the defining equation of V.import DynamicPolynomials: @polyvar\n\n@polyvar x y z\nf = x^2 + y^2 - z^2Let us sample the equation of a random line.L = randn(1,3) * [x; y; z]Now we compute the number of solutions to f=0 L=0.using HomotopyContinuation\nsolve([f; L])We find two distinct solutions and conclude that the degree of V is 2."
},

{
    "location": "examples.html#Using-different-types-of-homotopies-1",
    "page": "Examples",
    "title": "Using different types of homotopies",
    "category": "section",
    "text": "The following example is from Section 7.3 of[The numerical solution of systems of polynomials, Sommese, Wampler].Consider a triangle with sides a,b,c and let θ be the angle opposite of c. The goal is to compute θ from a,b,c. We define sθ := sin θ and cθ := cos θ. The polynomial corresponding system is.import DynamicPolynomials: @polyvar\n\na = 5\nb = 4\nc = 3\n\n@polyvar sθ cθ\nf = [cθ^2 + sθ^2 - 1, (a * cθ - b)^2 + (a * sθ)^2 - c^2]To set up a totaldegree homotopy of type StraightLineHomotopy we have to writeusing HomotopyContinuation\nH, s = totaldegree(StraightLineHomotopy, f)This sets up a homotopy H of the specified type using a random starting system that comes with a vector s of solutions. To solve for f = 0 we executesolve(H, s)If instead we wanted to use GeodesicOnTheSphere as homotopy type, we writeH, s = totaldegree(GeodesicOnTheSphere, f)\nsolve(H, s)The angles are of course only the real solutions of f = 0. We get them by usingsolution(ans, only_real=true)"
},

{
    "location": "examples.html#Using-different-types-of-pathrackers-1",
    "page": "Examples",
    "title": "Using different types of pathrackers",
    "category": "section",
    "text": "The following polynomial system is the example from Section 5.1 from[Decoupled molecules with binding polynomials of bidegree (n,2), Ren, Martini, Torres]It is called a binding polynomial.using HomotopyContinuation\nimport DynamicPolynomials: @polyvar\n\n@polyvar w[1:6]\n\nf = [\n    11*(2*w[1]+3*w[3]+5*w[5])+13*(2*w[2]+3*w[4]+5*w[6]),\n    11*(6*w[1]*w[3]+10*w[1]*w[5]+15*w[3]*w[5])+13*(6*w[2]*w[4]+10*w[2]*w[6]+15*w[4]*w[6]),\n    330*w[1]*w[3]*w[5]+390*w[2]*w[4]*w[6],\n    143*(2*w[1]*w[2]+3*w[3]*w[4]+5*w[5]*w[6]),\n    143*(6*w[1]*w[2]*w[3]*w[4]+10*w[1]*w[2]*w[5]*w[6]+15*w[3]*w[4]*w[5]*w[6]),\n    4290*w[1]*w[2]*w[3]*w[4]*w[5]*w[6]\n    ]Suppose we wanted to solve f(w)=a, wherea=[71, 73, 79, 101, 103, 107]To get an initial solution we compute a random forward solution.w_0 = randn(6)\na_0 = map(p -> p(w => w_0), f)Now we set up the homotopy.H = StraightLineHomotopy(f-a_0, f-a)and compute a backward solution with starting value w_0 bysolve(H, w_0)By default the solve function uses SphericalPredictorCorrector as the pathtracking routing. To use the AffinePredictorCorrector instead we must writesolve(H, w_0, AffinePredictorCorrector())The system f=0 has 72 simple non-real roots. The command    S = solve(f-a);\n    solutions(S, singular = false)however, only returns 62. The reason is that the remaining 10 solutions are ill-conditioned. We find all 72 solutions by    S = solve(f-a, singular_tol=1e8);\n    solutions(S, singular = false)The default of singular_tol in JuliaHomotopyContinuation is 1e4."
},

{
    "location": "examples.html#R-Serial-Link-Robots-1",
    "page": "Examples",
    "title": "6-R Serial-Link Robots",
    "category": "section",
    "text": "The following example is from Section 9.4 of[The numerical solution of systems of polynomials, Sommese, Wampler].Consider a robot that consists of 7 links connected by 6 joints. The first link is fixed on the ground and the last link has a hand. The problem of determining the position of the hand when knowing the arrangement of the joints is called forward problem. The problem of determining any arrangement of joints that realized a fixed position of the hand is called backward problem. Let us denote by z_1z_6 the unit vectors that point in the direction of the joint axes.  They satisfy the following polynomial equationsz_i  z_i = 1\n\nz_i  z_i+1 = cos _i\n\na_1 * z_1  z_2 +  + a_5 * z_5  z_6 + a_6 * z_2 +  + a_9 * z_5 = pfor some (a) and a known p (see the aforementioned reference for a detailed explanation on how these numbers are to be interpreted).In this notation the forward problem consists of computing (a) given the z_i and p. The backward problem consists of computing  z_i that realize some fixed (az_1z_6) (knowing z_1z_6 means that the position where the robot is attached to the ground  and the position where its hand should be are fixed).We now compute first a forward solution (_0 a_0), and then use (_0 a_0) to compute a solution for the backward problem imposed by some random ( a).using HomotopyContinuation\nimport DynamicPolynomials: @polyvar\n\n@polyvar z2[1:3] z3[1:3] z4[1:3] z5[1:3]\nz1 = [1, 0, 0]\nz6 = [1, 0, 0]\np = [1, 1, 0]\nz = [z1, z2, z3, z4, z5, z6]\n\nf = [z[i] ⋅ z[i] for i=2:5]\ng = [z[i] ⋅ z[i+1] for i=1:5]\nh = hcat([[z[i] × z[i+1] for i=1:5]; [z[i] for i=2:5]]...)\n\nα = randexp(5)\na = randexp(9)Let us compute a random forward solution.z_0=rand(3,4); # Compute a random assignment for the variable z\nfor i = 1:4\n    z_0[:,i] = z_0[:,i]./ norm(z_0[:,i]) # normalize the columns of z_0 to norm 1\nendWe want to compute the angles arccos g(z_0).z_0 = vec(z_0) # vectorize z_0, because the evaluate function takes vectors as input\n\n# compute the forward solution of α\nα_0 = map(p -> acos(p([z2; z3; z4; z5] => z_0)), g)\n\n# evaluate h at z_0\nh_0 = map(p -> p([z2; z3; z4; z5] => z_0), h)\na_0 = h_0\\pNow we have forward solutions _0 and a_0. From this we construct the following StraightLineHomotopy.H = StraightLineHomotopy([f-1; g-cos.(α_0); h*a_0-p], [f-1; g-cos.(α); h*a-p])To compute a backward solution with starting value z_0 we finally executesolve(H, z_0)To compute all the backward solutions we may perform a totaldegree homotopy. Although the Bezout number of the system is 1024 the generic number of solutions is 16. We find all 16 solutions byH, s = totaldegree(StraightLineHomotopy, [f-1; g-cos.(α_0); h*a_0-p])\nsolutions(solve(H, s), singular=false)On a MacBook Pro with 2,6 GHz Intel Core i7 and 16 GB RAM memory the above operation takes about 572 seconds. With parallel computing provided by the addprocs() command in Julia it takes about 93 seconds."
},

{
    "location": "Homotopy.html#",
    "page": "Setting up homotopies",
    "title": "Setting up homotopies",
    "category": "page",
    "text": ""
},

{
    "location": "Homotopy.html#Setting-up-homotopies-1",
    "page": "Setting up homotopies",
    "title": "Setting up homotopies",
    "category": "section",
    "text": "Homotopies.jl is a package for constructing (polynomial) homotopies H(xt). For the convient use we export in HomotopyContinuation every function from Homotopies.Each homotopy has the same Interface so that you can switch easily between different homotopy types. Based on this interface there are also some convenient higher level constructs provided; e.g., the construction of a total degree system and its start solutions.Homotopies.jl provides an interface to DynamicPolynomials.jl for human-readable input and output. Most of the examples in this introduction are written with DynamicPolynomials.jl . Internally, Homotopies.jl uses FixedPolynomials.jl for fast evaluation."
},

{
    "location": "Homotopy.html#Example-1",
    "page": "Setting up homotopies",
    "title": "Example",
    "category": "section",
    "text": "As an example we construct a homotopy between the polynomial systemsf= beginbmatrix x + y^3  x^2y-2yendbmatrixquad  \ng= beginbmatrixx^3+2 y^3+2endbmatrixCurrently, there are two types of homotopies implemented:StraightLineHomotopy\nGeodesicOnTheSphereThe code to initialize a StraightLineHomotopy is as follows.using HomotopyContinuation\nimport DynamicPolynomials: @polyvar # @polyvar is a function for initializing variables.\n@polyvar x y # initilize the variables x y\n\nf = [x + y^3, x^2*y-2y]\ng = [x^3+2, y^3+2]\n\nH = StraightLineHomotopy(f, g) # H is now StraightLineHomotopy{Int64}\n\n# to avoid unnecessary conversions one could also have\nH = StraightLineHomotopy{Complex128}([x + y^3, x^2*y-2y], [x^3+2, y^3+2])\n\n# we can now evaluate H\nevaluate(H, rand(Complex128, 2), 0.42)\n# or alternatively\nH(rand(Complex128, 2), 0.42)"
},

{
    "location": "Homotopy.html#Homotopies.StraightLineHomotopy",
    "page": "Setting up homotopies",
    "title": "Homotopies.StraightLineHomotopy",
    "category": "Type",
    "text": "StraightLineHomotopy(start, target)\n\nConstruct the homotopy t * start + (1-t) * target.\n\nstart and target have to match and to be one of the following\n\nVector{<:MP.AbstractPolynomial} where MP is MultivariatePolynomials\nMP.AbstractPolynomial\nVector{<:FP.Polynomial} where FP is FixedPolynomials\n\nStraightLineHomotopy{T}(start, target)\n\nYou can also force a specific coefficient type T.\n\n\n\n"
},

{
    "location": "Homotopy.html#Homotopies.GeodesicOnTheSphere",
    "page": "Setting up homotopies",
    "title": "Homotopies.GeodesicOnTheSphere",
    "category": "Type",
    "text": "GeodesicOnTheSphere(start, target)\n\nHomotopy is the geodesic from g=start/|start| (t=1) to f=target/|target| (t=0):\n\nH(x,t) = (cos(tα) - sin (tα)cos(α)/sin(α)) f + sin(tα) / sin(α) * g\n\nwhere  = cos fg. The constructor automatically homgenizes start and target.\n\nstart and target have to match and to be one of the following\n\nVector{<:MP.AbstractPolynomial} where MP is MultivariatePolynomials\nMP.AbstractPolynomial\nVector{<:FP.Polynomial} where FP is FixedPolynomials\n\nGeodesicOnTheSphere{T}(start, target)\n\nYou can also force a specific coefficient type T.\n\n\n\n"
},

{
    "location": "Homotopy.html#Homotopies-1",
    "page": "Setting up homotopies",
    "title": "Homotopies",
    "category": "section",
    "text": "The following homotopies are implemented. They are subtypes of AbstractPolynomialHomotopyStraightLineHomotopy\nGeodesicOnTheSphere"
},

{
    "location": "Homotopy.html#higherlevelconstructs-1",
    "page": "Setting up homotopies",
    "title": "Higher level constructs",
    "category": "section",
    "text": ""
},

{
    "location": "Homotopy.html#Homotopies.totaldegree",
    "page": "Setting up homotopies",
    "title": "Homotopies.totaldegree",
    "category": "Function",
    "text": "totaldegree(H::Type{AbstractPolynomialHomotopy}, F, [unitroots=false])\n\nConstruct a  total degree homotopy of type H with F and an iterator of its solutions. This is the homotopy with start system\n\nbeginalign*\n    z_1^d_1 - b_1\n    z_1^d_2 - b_2\n    vdots \n    z_n^d_n - b_n\nendalign*\n\nand target system F, where d_i is the degree of the i-th polynomial of F. If unitroots=true then b_i=1 otherwise b_i is a random complex number (with real and imaginary part in the unit interval).\n\nExample\n\nH, startsolutions = totaldegree(StraightLineHomotopy{Complex128}, [x^2+y+1, x^3*y-2])\n\n\n\n"
},

{
    "location": "Homotopy.html#Homotopies.TotalDegreeSolutionIterator",
    "page": "Setting up homotopies",
    "title": "Homotopies.TotalDegreeSolutionIterator",
    "category": "Type",
    "text": "TotalDegreeSolutionIterator(degrees, b)\n\nGiven the Vectors degrees and b TotalDegreeSolutionIterator enumerates all solutions of the system\n\nbeginalign*\n    z_1^d_1 - b_1 = 0 \n    z_1^d_2 - b_2 = 0 \n    vdots \n    z_n^d_n - b_n = 0 \nendalign*\n\nwhere d_i is degrees[i] and b_i is b[i].\n\n\n\n"
},

{
    "location": "Homotopy.html#Homotopies.totaldegree_startsystem",
    "page": "Setting up homotopies",
    "title": "Homotopies.totaldegree_startsystem",
    "category": "Function",
    "text": "totaldegree_startsystem(F::Vector{FP.Polynomial{<:Complex}}, [unit_roots=false])\n\nReturn the system\n\nbeginalign*\n    z_1^d_1 - b_1\n    z_1^d_2 - b_2\n    vdots \n    z_n^d_n - b_n\nendalign*\n\nwhere d_i is the degree of the i-th polynomial of F and an iterator of its solutions. If unitroots=true then b_i=1 otherwise b_i is a random complex number (with real and imaginary part in the unit interval).\n\n\n\n"
},

{
    "location": "Homotopy.html#totaldegree-1",
    "page": "Setting up homotopies",
    "title": "Total degree homotopy",
    "category": "section",
    "text": "totaldegree\nTotalDegreeSolutionIterator\ntotaldegree_startsystem"
},

{
    "location": "Homotopy.html#Homotopies.randomhomotopy",
    "page": "Setting up homotopies",
    "title": "Homotopies.randomhomotopy",
    "category": "Function",
    "text": "randomhomotopy(::Type{AbstractPolynomialHomotopy{T}}, size::Int; kwargs...)\n\nCreate a total degree homotopy where the target system is a randomsystem(T, size, size; kwargs...).\n\nExample\n\njulia> H, solutions = randomhomotopy(StraightLineHomotopy{Complex128}, 2, mindegree=3, maxdegree=6);\njulia> length(H)\n3\njulia> nvariables(H)\n3\n\n\n\n"
},

{
    "location": "Homotopy.html#Homotopies.randomsystem",
    "page": "Setting up homotopies",
    "title": "Homotopies.randomsystem",
    "category": "Function",
    "text": "randomsystem([T=Complex128,] nequations::Int, nvars::Int; mindegree=0, maxdegree=5, rng=Base.Random.GLOBAL_RNG, density=rand() * 0.8 + 0.1)\n\nCreates a random polynomial system of nequations equations with nvars variables (named x_1, ...x_nvars). Each polynomial has a total degree uniformly drawn from mindegree maxdegree. The coefficients are drawn independently from the given rng. With density you can control how many coefficients are non-zero. A value of 1.0 creates a dense polynomial (i.e. every coefficient is non-zero). A value of 0.5 creates a polynomial where only half of all monomials are non zero.\n\nrandomsystem([T=Complex128,] degrees::Vector{Int}, variables::Vector{Symbol}; rng=N(0,1))\n\nCreate a random polynomial system with the given degrees and variables.\n\n\n\n"
},

{
    "location": "Homotopy.html#Random-homotopies-1",
    "page": "Setting up homotopies",
    "title": "Random homotopies",
    "category": "section",
    "text": "randomhomotopy\nrandomsystem"
},

{
    "location": "Homotopy.html#Interface-1",
    "page": "Setting up homotopies",
    "title": "Interface",
    "category": "section",
    "text": ""
},

{
    "location": "Homotopy.html#Homotopies.evaluate",
    "page": "Setting up homotopies",
    "title": "Homotopies.evaluate",
    "category": "Function",
    "text": "evaluate(H::AbstractPolynomialHomotopy, x, t)\n\nEvaluate the homotopy H at x to time t, i.e. H(xt).\n\nevaluate(H::AbstractPolynomialHomotopy, x, t, cfg::PolynomialHomotopyConfig)\n\nEvaluate the homotopy H at x to time t using the precompuated values in cfg. Note that this is significantly faster than evaluate(H, x, t).\n\n\n\n"
},

{
    "location": "Homotopy.html#Homotopies.evaluate!",
    "page": "Setting up homotopies",
    "title": "Homotopies.evaluate!",
    "category": "Function",
    "text": "evaluate!(u::Vector, H::AbstractPolynomialHomotopy, x, t)\n\nEvaluate the homotopy H at x to time t, i.e. H(xt), and store the result in u.\n\nevaluate!(u::AbstractVector, H::AbstractPolynomialHomotopy, x, t, cfg::PolynomialHomotopyConfig)\n\nEvaluate the homotopy H at x to time t using the precompuated values in cfg and store the result in u.\n\n\n\n"
},

{
    "location": "Homotopy.html#Evaluation-1",
    "page": "Setting up homotopies",
    "title": "Evaluation",
    "category": "section",
    "text": "evaluate\nevaluate!"
},

{
    "location": "Homotopy.html#Homotopies.jacobian",
    "page": "Setting up homotopies",
    "title": "Homotopies.jacobian",
    "category": "Function",
    "text": "jacobian(H::AbstractPolynomialHomotopy, x, t, cfg::PolynomialHomotopyConfig)\n\nCompute the jacobian of H at x and t using the precomputed values in cfg. The jacobian is constructed w.r.t. x, i.e. it doesn't contain the partial derivatives w.r.t. t.\n\n\n\n"
},

{
    "location": "Homotopy.html#Homotopies.jacobian!",
    "page": "Setting up homotopies",
    "title": "Homotopies.jacobian!",
    "category": "Function",
    "text": "jacobian!(u, H::AbstractPolynomialHomotopy, x, t, cfg::PolynomialHomotopyConfig)\n\nCompute the jacobian of H at x and t using the precomputed values in cfg and store the result in u.\n\njacobian!(r::JacobianDiffResult, H::AbstractPolynomialHomotopy, x, t, cfg::PolynomialHomotopyConfig)\n\nCompute H(x t) and the jacobian of H at x and t at once using the precomputated values in cfg and store thre result in r. This is faster than computing both values separetely.\n\nExample\n\ncfg = PolynomialHomotopyConfig(H)\nr = JacobianDiffResult(cfg)\njacobian!(r, H, x, t, cfg)\n\nvalue(r) == H(x, t)\njacobian(r) == jacobian(H, x, t, cfg)\n\n\n\n"
},

{
    "location": "Homotopy.html#Homotopies.dt",
    "page": "Setting up homotopies",
    "title": "Homotopies.dt",
    "category": "Function",
    "text": "dt(H::AbstractPolynomialHomotopy, x, t, cfg::PolynomialHomotopyConfig)\n\nCompute the derivative of H w.r.t. t at x and t using the precomputed values in cfg.\n\n\n\n"
},

{
    "location": "Homotopy.html#Homotopies.dt!",
    "page": "Setting up homotopies",
    "title": "Homotopies.dt!",
    "category": "Function",
    "text": "dt!(u, H::AbstractPolynomialHomotopy, x, t, cfg::PolynomialHomotopyConfig)\n\nCompute the derivative of H w.r.t. t at x and t using the precomputed values in cfg and store the result in u.\n\ndt!(r::DtDiffResult, H::AbstractPolynomialHomotopy, x, t, cfg::PolynomialHomotopyConfig)\n\nCompute the derivative of H w.r.t. t at x and t using the precomputed values in cfg and store the result in r. This is faster than computing both values separetely.\n\nExample\n\ncfg = PolynomialHomotopyConfig(H)\nr = DtDiffResult(cfg)\ndt!(r, H, x, t, cfg)\n\nvalue(r) == H(x, t)\ndt(r) == dt(H, x, t, cfg)\n\n\n\n"
},

{
    "location": "Homotopy.html#Differentiation-1",
    "page": "Setting up homotopies",
    "title": "Differentiation",
    "category": "section",
    "text": "jacobian\njacobian!\ndt\ndt!"
},

{
    "location": "Homotopy.html#Homotopies.homogenize",
    "page": "Setting up homotopies",
    "title": "Homotopies.homogenize",
    "category": "Function",
    "text": "homogenize(H::AbstractPolynomialHomotopy)\n\nHomogenize the homotopy H. This adds an additional variable. If H is already homogenized, this is the identity.\n\n\n\n"
},

{
    "location": "Homotopy.html#Homotopies.dehomogenize",
    "page": "Setting up homotopies",
    "title": "Homotopies.dehomogenize",
    "category": "Function",
    "text": "dehomogenize(H::AbstractPolynomialHomotopy)\n\nDehomogenize the homotopy H. This removes the first variable. If H is not homogenized, this is the identity.\n\n\n\n"
},

{
    "location": "Homotopy.html#Homotopies.ishomogenized",
    "page": "Setting up homotopies",
    "title": "Homotopies.ishomogenized",
    "category": "Function",
    "text": "ishomogenized(H::AbstractPolynomialHomotopy)\n\nCheck whether the homotopy H was homogenized.\n\n\n\n"
},

{
    "location": "Homotopy.html#Homotopies.ishomogenous",
    "page": "Setting up homotopies",
    "title": "Homotopies.ishomogenous",
    "category": "Function",
    "text": "ishomogenous(H::AbstractPolynomialHomotopy)\n\nCheck whether the homotopy H is homogenous. This does not imply that H was homogenized.\n\n\n\n"
},

{
    "location": "Homotopy.html#Homogenization-1",
    "page": "Setting up homotopies",
    "title": "Homogenization",
    "category": "section",
    "text": "homogenize\ndehomogenize\nishomogenized\nishomogenous"
},

{
    "location": "Homotopy.html#Homotopies.nvariables",
    "page": "Setting up homotopies",
    "title": "Homotopies.nvariables",
    "category": "Function",
    "text": "nvariables(H::AbstractPolynomialHomotopy)\n\nThe number of variables which H expects as input, i.e. to evaluate H(x,t) x has to be a vector of length nvariables(H).\n\n\n\n"
},

{
    "location": "Homotopy.html#Homotopies.weylnorm",
    "page": "Setting up homotopies",
    "title": "Homotopies.weylnorm",
    "category": "Function",
    "text": "weylnorm(H::AbstractPolynomialHomotopy)\n\nCreates a function with variable t that computes the Weyl norm (or Bombieri norm) of H(xt). See here for details about the Weyl norm.\n\n\n\n"
},

{
    "location": "Homotopy.html#Homotopies.gammatrick!",
    "page": "Setting up homotopies",
    "title": "Homotopies.gammatrick!",
    "category": "Function",
    "text": "gammatrick!(H::AbstractPolynomialHomotopy{Complex} [, seed::Int]])\n\nScale the coefficients of the start system of H with a random complex number picked uniformly from the (complex) unit circle. Use this to make the paths z(t) generic.\n\ngammatrick!(H::AbstractPolynomialHomotopy{Complex}, γ::Complex)\n\nYou can also pass a scaling factor directly.\n\n\n\n"
},

{
    "location": "Homotopy.html#Homotopies.gammatrick",
    "page": "Setting up homotopies",
    "title": "Homotopies.gammatrick",
    "category": "Function",
    "text": "gammatrick(H::AbstractPolynomialHomotopy{Complex} , γ::Number)\n\nScale the coefficients of the start system of H with γ.\n\ngammatrick(H::AbstractPolynomialHomotopy{Complex})\n\nA a random complex number γ is picked uniformly from the (complex) unit circle and then scale the coefficients of the start system of H with γ. This returns the new H and γ.\n\n\n\n"
},

{
    "location": "Homotopy.html#Misc-1",
    "page": "Setting up homotopies",
    "title": "Misc",
    "category": "section",
    "text": "nvariables\nweylnorm\ngammatrick!\ngammatrick"
},

{
    "location": "Homotopy.html#Homotopies.κ",
    "page": "Setting up homotopies",
    "title": "Homotopies.κ",
    "category": "Function",
    "text": "κ(H, z, t, cfg)\n\nComputes the condition number of H at (z, t) (with config cfg). See Condition^[1] for details\n\nProposition 16.10: κ(f,z) := ‖f‖ ‖ Df(z)^† diag(‖ z ‖^{d_i-1}) ‖\n\n[1]: Condition, Bürgisser and Cucker\n\n\n\n"
},

{
    "location": "Homotopy.html#Homotopies.κ_norm",
    "page": "Setting up homotopies",
    "title": "Homotopies.κ_norm",
    "category": "Function",
    "text": "κ_norm(H, z, t, cfg)\n\nComputes the condition number of H at (z, t) (with config cfg). See Condition^[1] for details\n\nEq. (16.11): κ_norm(f,z) := ‖f‖ ‖ Df(z)^† diag(√{d_i}‖ z ‖^{d_i-1}) ‖\n\n[1]: Condition, Bürgisser and Cucker\n\n\n\n"
},

{
    "location": "Homotopy.html#Homotopies.μ_norm",
    "page": "Setting up homotopies",
    "title": "Homotopies.μ_norm",
    "category": "Function",
    "text": "μ_norm(H, z, t, cfg)\n\nComputes the condition number of H at (z, t) (with config cfg). See Condition^[1] for details\n\nDefinition 16.43: μ_norm(f,z) := ‖f‖ ‖ (Df(z)-(T_z))^{-1} diag(√{d_i}‖ z ‖^{d_i-1}) ‖\n\n[1]: Condition, Bürgisser and Cucker\n\n\n\n"
},

{
    "location": "Homotopy.html#Condition-numbers-1",
    "page": "Setting up homotopies",
    "title": "Condition numbers",
    "category": "section",
    "text": "κ\nκ_norm\nμ_norm"
},

{
    "location": "solve.html#",
    "page": "Solving homotopies",
    "title": "Solving homotopies",
    "category": "page",
    "text": ""
},

{
    "location": "solve.html#solve-1",
    "page": "Solving homotopies",
    "title": "Solving homotopies",
    "category": "section",
    "text": ""
},

{
    "location": "solve.html#HomotopyContinuation.solve",
    "page": "Solving homotopies",
    "title": "HomotopyContinuation.solve",
    "category": "Function",
    "text": "    solve(H::AbstractHomotopy, startvalues_s, [algorithm]; kwargs...)\n\nSolve the homotopy H via homotopy continuation with the given startvalues_s and the given algorithm.\n\nkwargs are the keyword arguments for the solver options.\n\nThe default pathtracking algorithm is SphericalPredictorCorrector().\n\nTo specify another pathracking algorithm, e.g.AffinePredictorCorrector(), write\n\n    solve(H::AbstractHomotopy, startvalues_s, AffinePredictorCorrector(); kwargs...)\n\nThe function also takes polynomials as inputs:\n\n    solve(f::Vector{<:MP.AbstractPolynomial{T}})\n\nsolves the polynomial system f via a totaldegree homotopy of type StraightLineHomotopy and the SphericalPredictorCorrector pathtracking routine.\n\nTo specify homotopy and pathtracker, use\n\n    solve(f::Vector{<:MP.AbstractPolynomial{T}}, [homotopy], [algorithm]; kwargs...)\n\nDefault is homotopy = StraightLineHomotopy and algorithm = SphericalPredictorCorrector. For instance,\n\n    solve(f, GeodesicOnTheSphere, AffinePredictorCorrector())\n\nsolves f=0 with a GeodesicOnTheSphere homotopy and the AffinePredictorCorrector pathtracking routine.\n\n\n\n"
},

{
    "location": "solve.html#The-solve-function-1",
    "page": "Solving homotopies",
    "title": "The solve function",
    "category": "section",
    "text": "The solve function solves homotopies with given starting values.solve"
},

{
    "location": "solve.html#HomotopyContinuation.Solver",
    "page": "Solving homotopies",
    "title": "HomotopyContinuation.Solver",
    "category": "Type",
    "text": "Solver(homotopy, pathtracking_algorithm=SphericalPredictorCorrector(), endgame=CauchyEndgame(); kwargs...)\n\nCreate a mutable Solver struct. This contains a Pathtracker and an Endgamer, everything you need to solve the given homotopy. Solver supports the following options:\n\nendgame_start=0.1: Where the endgame starts\nabstol=1e-12: The desired accuracy of the final roots\nat_infinity_tol=1e-10: An point is at infinity if the maginitude of the homogenous variable\n\nis less than at_infinity_tol.\n\nsingular_tol=1e4: If the winding number is 1 but the condition number is larger than\n\nsingular_tol then the root is declared as singular.\n\nrefinement_maxiters=100: The maximal number of newton iterations to achieve abstol.\nverbose=false: Print additional warnings / informations\napply_gammatrick=true: This modifies the start system to make it generic.\ngamma=apply_gammatrick ? exp(im*2π*rand()) : complex(1.0): You can overwrite the default gamma.   This is useful if you want to rerun only some paths.\npathcrossing_tolerance=1e-8: The tolerance for when two paths are considered to be crossed.\npathcrossing_check=true: Enable the pathcrossing check.\nparallel_type=:pmap: Currently there are two modes: :pmap will use pmap for parallelism\n\nand :none will use the standard map. :pmap is by defautl enabled since it works reliable, but if you develop new algorithms you probably want to disable parallelism.\n\nbatch_size=1: The batch_size for pmap if parallel_type is :pmap.\n\nFor instance, to solve the homotopy H with starting values s with no endgame and a singular tolerance of 1e5, write\n\n    solve(H, s, endgame_start=0.0, singular_tol=1e5)\n\nTo solve the polynomial system f with the same options write\n\n    solve(f, endgame_start=0.0, singular_tol=1e5)\n\n\n\n"
},

{
    "location": "solve.html#solveroptions-1",
    "page": "Solving homotopies",
    "title": "Solver options",
    "category": "section",
    "text": "Solver"
},

{
    "location": "solve.html#HomotopyContinuation.Result",
    "page": "Solving homotopies",
    "title": "HomotopyContinuation.Result",
    "category": "Type",
    "text": "Result(pathresults, solver)\n\nA thin wrapper around the PathResults of the Solver instance. Result behaves like an array of PathResults but also contains some additional informations. For example you can obtain the γ which was used for the gammatrick.\n\n\n\n"
},

{
    "location": "solve.html#result-1",
    "page": "Solving homotopies",
    "title": "Result",
    "category": "section",
    "text": "The solve function returns a Result struct:Result"
},

{
    "location": "solve.html#HomotopyContinuation.PathResult",
    "page": "Solving homotopies",
    "title": "HomotopyContinuation.PathResult",
    "category": "Type",
    "text": "PathResult(startvalue, pathtracker_result, endgamer_result, solver)\n\nConstruct a PathResult for a given startvalue. pathtracker_result is the PathtrackerResult until the endgame radius is reached. endgamer_result is the EndgamerResult resulting from the corresponding endgame.\n\nA PathResult contains:\n\nreturncode: One of :success, :at_infinity or any error code from the EndgamerResult\nsolution::Vector{T}: The solution vector. If the algorithm computed in projective space\n\nand the solution is at infinity then the projective solution is given. Otherwise an affine solution is given if the startvalue was affine and a projective solution is given if the startvalue was projective.\n\nresidual::Float64: The value of the infinity norm of H(solution, 0).\nnewton_residual: The value of the 2-norm of J_H(textsolution)^-1H(textsolution 0)\nlog10_condition_number: A high condition number indicates singularty. See Homotopies.κ for details.   The value is the logarithmic condition number (with base 10).\nwindingnumber: The estimated winding number\nangle_to_infinity: The angle to infinity is the angle of the solution to the hyperplane where the homogenizing coordinate is 0.\nreal_solution: Indicates whether the solution is real given the defined tolerance at_infinity_tol (from the solver options).\nstartvalue: The startvalue of the path\niterations: The number of iterations the pathtracker needed.\nendgame_iterations: The number of steps in the geometric series the endgamer did.\nnpredictions: The number of predictions the endgamer did.\npredictions: The predictions of the endgamer.\n\n\n\n"
},

{
    "location": "solve.html#PathResult-1",
    "page": "Solving homotopies",
    "title": "PathResult",
    "category": "section",
    "text": "For each tracked path there is a PathResult:PathResult"
},

{
    "location": "solve.html#HomotopyContinuation.solutions",
    "page": "Solving homotopies",
    "title": "HomotopyContinuation.solutions",
    "category": "Function",
    "text": "solutions(r::Result; success=true, at_infnity=true, only_real=false, singular=true)\n\nFilters the solutions which satisfy the constraints.\n\n\n\n"
},

{
    "location": "solve.html#solutions-1",
    "page": "Solving homotopies",
    "title": "The solutions function",
    "category": "section",
    "text": "The solution function helps to extract information from a Result:solutions"
},

{
    "location": "pathtracker.html#",
    "page": "Pathtracking",
    "title": "Pathtracking",
    "category": "page",
    "text": ""
},

{
    "location": "pathtracker.html#HomotopyContinuation.Pathtracker",
    "page": "Pathtracking",
    "title": "HomotopyContinuation.Pathtracker",
    "category": "Type",
    "text": "Pathtracker(H::AbstractHomotopy{T}, alg, [HT::Type=widen(T)]; kwargs...)\n\nConstruct a Pathtracker object. This contains all informations to track a single path for H with the given pathtracking algorithm alg. The optional type HT is used if the pathracker decides to switch to a high precision mode.\n\nThe following keyword arguments are supported:\n\npath_precision=1e-6: The precision for which a correction step is decleared successfull.\ncorrector_maxiters=3: The maximal number of correction iterations. A higher value as 3 is not recommended.\ninitial_steplength=0.1: The initial steplength a preditor-corrector algorithm uses.\nconsecutive_successfull_steps_until_steplength_increase=3:   The number of consecutive successfull steps until the step length is increased multiplied   with the factor steplength_increase_factor.\nsteplength_increase_factor=2.0\nsteplength_decrease_factor=inv(steplength_increase_factor): If a correction step fails the step length is multiplied   with this factor.\nmaxiters=10_000: The maximum number of iterations.\nvebose=false: Print additional informations / warnings during the computation.\n\n\n\n"
},

{
    "location": "pathtracker.html#Pathtracking-routines-1",
    "page": "Pathtracking",
    "title": "Pathtracking routines",
    "category": "section",
    "text": "Pathtracking is at the core of each homotopy continuation method. It is the routine to track a given homotopy H(x t) from a start value x_1 at time t_1 to a target value x_0 at time t_0.At the heart of the pathtracking routine is the  mutable struct Pathtracker.Pathtracker"
},

{
    "location": "pathtracker.html#Examples-1",
    "page": "Pathtracking",
    "title": "Examples",
    "category": "section",
    "text": "The follwing example demonstrates the usual workflow. You first create a Pathtracker object, then you can track a path from a given start value and finally you create a PathtrackerResult.pathtracker = Pathtracker(H, SphericalPredictorCorrector())\ntrack!(pathtracker, x, 1.0, 0.0)\nresult = PathtrackerResult(pathtracker)You can reuse (and should!) resuse a Pathtracker for multiple pathspathtracker = Pathtracker(H, SphericalPredictorCorrector())\nresults = map(xs) do x\n  track!(pathtracker, x, 1.0, 0.0)\n  PathtrackerResult(pathtracker)\nendPathtracker also supports the iterator interface. This returns the complete Pathtracker object at each iteration. This enables all sort of nice features. For example you could store the actual path the pathtracker takes:pathtracker = Pathtracker(H, SphericalPredictorCorrector())\nsetup_pathtracker!(pathtracker, x, 1.0, 0.0)\npath = []\nfor t in pathtracker\n  push!(path, current_value(t))\nend"
},

{
    "location": "pathtracker.html#HomotopyContinuation.PathtrackerResult",
    "page": "Pathtracking",
    "title": "HomotopyContinuation.PathtrackerResult",
    "category": "Type",
    "text": "PathtrackerResult(pathtracker, extended_analysis=false)\n\nReads the result from the current pathtracker state. A PathtrackerResult contains:\n\nreturncode: One of :max_iterations, :singularity, :invalid_startvalue, :success.\nsolution::Vector{T}: The solution.\nresidual::Float64: The value of the infinity norm of H(solution, 0).\niterations: The number of iterations the pathtracker needed.\nangle_to_infinity: The angle to infinity is the angle of the solution to the hyperplane where the homogenizing coordinate is 0.\n\nIf extended_analysis=true there is also:\n\nnewton_residual: The value of the 2-norm of J_H(textsolution)^-1H(textsolution 0)\ncondition_number: A high condition number indicates singularty. See Homotopies.κ for details.\n\n\n\n"
},

{
    "location": "pathtracker.html#Result-1",
    "page": "Pathtracking",
    "title": "Result",
    "category": "section",
    "text": "See also here.PathtrackerResult"
},

{
    "location": "pathtracker.html#HomotopyContinuation.SphericalPredictorCorrector",
    "page": "Pathtracking",
    "title": "HomotopyContinuation.SphericalPredictorCorrector",
    "category": "Type",
    "text": "SphericalPredictorCorrector\n\nThis algorithm uses as an prediction step an explicit Euler method. For the prediction and correction step the Jacobian is augmented by Hermitian transposed x^* to get a square system. Therefore the prediciton step looks like\n\nx_k+1 = x_k - tbeginbmatrix\n    J_H(x t) \n    x^*\nendbmatrix^-1\nfracHt(x t)\n\nand the correction step looks like\n\nx_k+1 = x_k+1 - beginbmatrix\n    J_H(x t) \n    x^*\nendbmatrix^-1\nH(x t)\n\nAfter each prediciton and correction the algorithm normalizes x again, i.e. x is then a point on the sphere again.\n\nThis algorithm tracks the path in the projective space.\n\n\n\n"
},

{
    "location": "pathtracker.html#HomotopyContinuation.AffinePredictorCorrector",
    "page": "Pathtracking",
    "title": "HomotopyContinuation.AffinePredictorCorrector",
    "category": "Type",
    "text": "AffinePredictorCorrector\n\nThis algorithm uses as an prediction step an explicit Euler method. Therefore the prediciton step looks like\n\nx_k+1 = x_k - tJ_H(x t)^-1fracHt(x t)\n\nand the correction step looks like\n\nx_k+1 = x_k+1 - J_H(x t)^-1H(x t)\n\nThis algorithm tracks the path in the affine space.\n\n\n\n"
},

{
    "location": "pathtracker.html#Algorithms-1",
    "page": "Pathtracking",
    "title": "Algorithms",
    "category": "section",
    "text": "Currently, the following pathtracking routines are implementedSphericalPredictorCorrector\nAffinePredictorCorrector"
},

{
    "location": "pathtracker.html#HomotopyContinuation.track!",
    "page": "Pathtracking",
    "title": "HomotopyContinuation.track!",
    "category": "Function",
    "text": "track!(pathtracker, x0, s_start, s_target)\n\nTrack a startvalue x0 from s_start to s_target using the given pathtracker.\n\ntrack!(pathtracker)\n\nRun the given pathtracker. You can use this in combination with setup_pathtracker!.\n\n\n\n"
},

{
    "location": "pathtracker.html#HomotopyContinuation.setup_pathtracker!",
    "page": "Pathtracking",
    "title": "HomotopyContinuation.setup_pathtracker!",
    "category": "Function",
    "text": "setup_pathtracker!(tracker, x0, s_start, s_end)\n\nReset the given pathtracker tracker and set it up to track x0 form s_start to s_end.\n\n\n\n"
},

{
    "location": "pathtracker.html#HomotopyContinuation.current_value",
    "page": "Pathtracking",
    "title": "HomotopyContinuation.current_value",
    "category": "Function",
    "text": "current_value(pathtracker)\n\nGet the current value of the pathtracker.\n\n\n\n"
},

{
    "location": "pathtracker.html#Reference-1",
    "page": "Pathtracking",
    "title": "Reference",
    "category": "section",
    "text": "track!\nsetup_pathtracker!\ncurrent_value"
},

{
    "location": "endgame.html#",
    "page": "Endgame",
    "title": "Endgame",
    "category": "page",
    "text": ""
},

{
    "location": "endgame.html#HomotopyContinuation.Endgamer",
    "page": "Endgame",
    "title": "HomotopyContinuation.Endgamer",
    "category": "Type",
    "text": "Endgamer(endgame_algorithm, pathtracker, [x, R]; kwargs...)\n\nConstruct an Endgamer object. The Endgamer 'plays' the endgame with the given endgame_algorithm and uses the given pathtracker to move forward. The endgame is start at x to time R (the endgame radius). In each iteration the endgame moves forward and then performs one iteration of the endgame algorithm. In each iteration we could get another prediction and an estimate of the winding number. Convergence is declared if two consecutive predictions are smaller than a defined tolerance (endgame_abstol).\n\nThe following options are available:\n\ngeometric_series_factor=0.5: The Endgamer moves forward using the geometric series   ^kR where  is geometric_series_factor.\nmax_winding_number=16 the maximal winding number we allow. If we get a higher winding number\n\nthe path is declared failed.\n\nendgame_abstol=pathtracker.options.abstol: The tolerance necessary to declare convergence.\n\nEndgamer supports similar to Pathtracker an iterator interface.\n\n\n\n"
},

{
    "location": "endgame.html#Endgame-1",
    "page": "Endgame",
    "title": "Endgame",
    "category": "section",
    "text": "Assume we want to find the all solutions of the polynomial (x-2)^4 with homotopy continuation. Then the pathtracker gets into severe trouble near the end of the path since the derivative is 0 at x=2.How do we solve that problem? The idea is to split the pathtracking into two parts. We first track the path x(t) until the so called endgame zone (starting by default at t=01). Then we switch to the endgame. The idea is to estimate the value of x(0.0) without tracking the path all the way to t=00(since this would fail due to a singular Jacobian).There are two well known endgame strategies. The Power Series Endgame and the Cauchy Endgame. Currently only the Cauchy Endgame is implemented.At the heart of the endgame routine is the  mutable struct Endgamer.Endgamer"
},

{
    "location": "endgame.html#Examples-1",
    "page": "Endgame",
    "title": "Examples",
    "category": "section",
    "text": "The follwing example demonstrates the usual workflow. You first create an Endgamer object, then you can track a path from a given start value and finally you create a EndgamerResult.endgamer = Endgamer(CauchyEndgame(), pathtracker)\nendgamer!(endgamer, x, 0.1)\nresult = EndgamerResult(endgamer)You can reuse (and should!) resuse an Endgamer for multiple pathsendgamer = Endgamer(CauchyEndgame(), pathtracker))\nresults = map(xs) do x\n  endgame!(endgamer, x, 0.1)\n  EndgamerResult(endgamer)\nend"
},

{
    "location": "endgame.html#HomotopyContinuation.EndgamerResult",
    "page": "Endgame",
    "title": "HomotopyContinuation.EndgamerResult",
    "category": "Type",
    "text": "EndgamerResult(endgamer, extended_analysis=false)\n\nReads the result from the current pathtracker state. A EndgamerResult contains:\n\nreturncode: One of :ill_conditioned_zone, :success, :windingnumber_too_high\nsolution::Vector{T}: The solution.\nstartvalue::Vector{T}: The solution.\nresidual::Float64: The value of the infinity norm of H(solution, 0).\niterations: The number of iterations the pathtracker needed.\nnpredictions: The number of predictions\npredictions: All predictions for further introspection (e.g. the path failed)\nangle_to_infinity: The angle to infinity is the angle of the solution to the hyperplane where the homogenizing coordinate is 0.\nwindingnumber: The estimated winding number\n\n\n\n"
},

{
    "location": "endgame.html#Result-1",
    "page": "Endgame",
    "title": "Result",
    "category": "section",
    "text": "EndgamerResult"
},

{
    "location": "endgame.html#HomotopyContinuation.CauchyEndgame",
    "page": "Endgame",
    "title": "HomotopyContinuation.CauchyEndgame",
    "category": "Type",
    "text": "CauchyEndgame(;kwargs...)\n\nThe main idea of the Cauchy Endgame is to use Cauchy's integral formula to predict the solution of the path x(t), i.e. x(0). At each iteration we are at some point (x t). We then track the polygon defined by te^i2kn until we end again at x. Here n is the number of samples we take per loop.\n\nThe following options are available:\n\nsamples_per_loop=8: The number of samples we take at one loop.\nloopclosed_tolerance=1e-5: The tolerance when a loop is considered closed.\nL=0.75 and K=05: These are paramters for heuristics. For more details see \"A Parallel Endgame \" by Bates, Hauenstein and Sommese [1],   page 8 and 9.\n\n[1]: Bates, Daniel J., Jonathan D. Hauenstein, and Andrew J. Sommese. \"A Parallel Endgame.\" Contemp. Math 556 (2011): 25-35.\n\n\n\n"
},

{
    "location": "endgame.html#Algorithms-1",
    "page": "Endgame",
    "title": "Algorithms",
    "category": "section",
    "text": "CauchyEndgame"
},

{
    "location": "endgame.html#HomotopyContinuation.endgame!",
    "page": "Endgame",
    "title": "HomotopyContinuation.endgame!",
    "category": "Function",
    "text": "endgame!(endgamer, x, R)\n\nPlay the endgame for x starting from time R.\n\nendgame!(endgamer)\n\nStart the endgamer. You probably want to setup things in prior with setup_endgamer!.\n\n\n\n"
},

{
    "location": "endgame.html#HomotopyContinuation.setup_endgamer!",
    "page": "Endgame",
    "title": "HomotopyContinuation.setup_endgamer!",
    "category": "Function",
    "text": "setup_endgamer!(endgamer, x, R)\n\nSetup endgamer to play the endgame starting from x at time R.\n\n\n\n"
},

{
    "location": "endgame.html#Reference-1",
    "page": "Endgame",
    "title": "Reference",
    "category": "section",
    "text": "endgame!\nsetup_endgamer!"
},

{
    "location": "set_up_homotopy.html#",
    "page": "How to set up your own homotopy",
    "title": "How to set up your own homotopy",
    "category": "page",
    "text": ""
},

{
    "location": "set_up_homotopy.html#How-to-set-up-your-own-homotopy-1",
    "page": "How to set up your own homotopy",
    "title": "How to set up your own homotopy",
    "category": "section",
    "text": "We shall illustrate how to set up your own homotopy by work out the following example.For a polynomial systems fg we want to define the homotopyH(xt) = t * f( U(t) x ) + (1 - t) * g( U(t) x )where U(t) is a random path in the space of unitary matrices with U(0) = U(1) = I, the identity matrix. I.e.,U(t) = U beginbmatrix\ncos(2 t)  -sin(2 t)  0 cdots  0\nsin(2 t)  cos(2 t)  0 cdots  0\n0  0  1 cdots  0\n0  0  0 cdots  0\n0  0  0 cdots  1\nendbmatrix U^Twith a random unitary matrix U.To start we make a copy of the file straigthline.jl (or of any other appropriate file like geodesic_on_the_sphere.jl) and rename it rotation_and_straightline.jl. Now we have to do two things:Adapt the struct and its constructors.\nAdapt the evaluation functions.\nInclude rotation_and_straightline.jl in Homotopies.jl."
},

{
    "location": "set_up_homotopy.html#Adapt-the-struct-and-its-constructors-1",
    "page": "How to set up your own homotopy",
    "title": "Adapt the struct and its constructors",
    "category": "section",
    "text": "We now assume that the file we copied is straigthline.jl. It is convenient to make a search-and-replace on StraightLineHomotopy and replace it by RotationAndStraightLine. First we adapt the contructor. Note that in the initialization of the struct we sample a random matrix and extract a unitary matrix U from its QR-decomposition. From this we define the function U(t) and save it together with its derivative in the struct.\nmutable struct RotationAndStraightLine{T<:Number} <: AbstractPolynomialHomotopy{T}\n    start::Vector{FP.Polynomial{T}}\n    target::Vector{FP.Polynomial{T}}\n    U::Function\n    U_dot::Function\n\n    function RotationAndStraightLine{T}(start::Vector{FP.Polynomial{T}}, target::Vector{FP.Polynomial{T}}) where {T<:Number}\n        @assert length(start) == length(target) \"Expected the same number of polynomials, but got $(length(start)) and $(length(target))\"\n\n\n        s_nvars = maximum(FP.nvariables.(start))\n        @assert all(s_nvars .== FP.nvariables.(start)) \"Not all polynomials of the start system have $(s_nvars) variables.\"\n\n        t_nvars = maximum(FP.nvariables.(target))\n        @assert all(t_nvars .== FP.nvariables.(target)) \"Not all polynomials of the target system have $(t_nvars) variables.\"\n\n        @assert s_nvars == t_nvars \"Expected start and target system to have the same number of variables, but got $(s_nvars) and $(t_nvars).\"\n\n        U = qrfact(randn(s_nvars,s_nvars) + im * randn(s_nvars,s_nvars))[:Q]\n\n        function U_fct(t)\n                (cos(2 * pi * t) - 1) .* U[:,1] * U[:,1]' - sin(2 * pi * t) .* U[:,2] * U[:,1]' + sin(2 * pi * t) .* U[:,1] * U[:,2]' + (cos(2 * pi * t) - 1) .* U[:,2] * U[:,2]' + eye(U)\n        end\n\n        function U_dot(t)\n                2 * pi .* (-sin(2 * pi * t) .* U[:,1] * U[:,1]' - cos(2 * pi * t) .* U[:,2] * U[:,1]' + cos(2 * pi * t) .* U[:,1] * U[:,2]' - sin(2 * pi * t) .* U[:,2] * U[:,2]')\n        end\n\n        new(start, target, U_fct, U_dot)\n    end\n\n    function RotationAndStraightLine{T}(start, target) where {T<:Number}\n        s, t = construct(T, start, target)\n        RotationAndStraightLine{T}(s, t)\n    end\nend\n\n\nfunction RotationAndStraightLine(start, target)\n    T, s, t = construct(start, target)\n    RotationAndStraightLine{T}(s, t)\nendThe conversion functions are adapted easily with copy-and-paste.#\n# SHOW\n#\nfunction Base.deepcopy(H::RotationAndStraightLine)\n    RotationAndStraightLine(deepcopy(H.start), deepcopy(H.target))\nend\n#\n# PROMOTION AND CONVERSION\n#ß\nBase.promote_rule(::Type{RotationAndStraightLine{T}}, ::Type{RotationAndStraightLine{S}}) where {S<:Number,T<:Number} = RotationAndStraightLine{promote_type(T,S)}\nBase.promote_rule(::Type{RotationAndStraightLine}, ::Type{S}) where {S<:Number} = RotationAndStraightLine{S}\nBase.promote_rule(::Type{RotationAndStraightLine{T}}, ::Type{S}) where {S<:Number,T<:Number} = RotationAndStraightLine{promote_type(T,S)}\nBase.convert(::Type{RotationAndStraightLine{T}}, H::RotationAndStraightLine) where {T} = RotationAndStraightLine{T}(H.start, H.target)"
},

{
    "location": "set_up_homotopy.html#Adapt-the-evaluation-functions.-1",
    "page": "How to set up your own homotopy",
    "title": "Adapt the evaluation functions.",
    "category": "section",
    "text": "The essential part of the homotopy struct are the evaluation functions. Here is where we define the orthogonal rotation.The function to be edited are evaluate, jacobian, dt and weylnorm. For fast evaluation there is a function evaluate_start_target that evaluates start and target system efficiently.The function that evaluates the homotopy at x at time t is#\n# EVALUATION + DIFFERENTATION\n#\nfunction evaluate!(u::AbstractVector, H::RotationAndStraightLine{T}, x::Vector, t::Number) where T\n    y = H.U(t) * x\n    for i = 1:length(H.target)\n        f = H.target[i]\n        g = H.start[i]\n        u[i] = (one(T) - t) * FP.evaluate(f, y) + t * FP.evaluate(g, y)\n    end\n    u\nend\n(H::RotationAndStraightLine)(x,t) = evaluate(H,x,t)\n\nfunction evaluate!(u::AbstractVector{T}, H::RotationAndStraightLine, x::Vector, t::Number, cfg::PolynomialHomotopyConfig, precomputed=false) where {T<:Number}\n    y = H.U(t) * x\n    evaluate_start_target!(cfg, H, y, precomputed)\n    u .= (one(t) - t) .* value_target(cfg) .+ t .* value_start(cfg)\nendThe derivative of the homotopy with respect to x isfunction jacobian!(u::AbstractMatrix, H::RotationAndStraightLine{T}, x::AbstractVector, t, cfg::PolynomialHomotopyConfig, precomputed=false) where {T<:Number}\n    U = H.U(t)\n    y = U * x\n    jacobian_start_target!(cfg, H, y, precomputed)\n    u .= ((one(t) - t) .* jacobian_target(cfg) .+ t .* jacobian_start(cfg)) * U\nend\n\nfunction jacobian!(r::JacobianDiffResult, H::RotationAndStraightLine{T}, x::AbstractVector, t, cfg::PolynomialHomotopyConfig, precomputed=false) where {T<:Number}\n    U = H.U(t)\n    y = U * x\n    evaluate_and_jacobian_start_target!(cfg, H, y)\n\n    r.value .= (one(t) - t) .* value_target(cfg) .+ t .* value_start(cfg)\n    r.jacobian .= ((one(t) - t) .* jacobian_target(cfg) .+ t .* jacobian_start(cfg)) * U\n    r\nendThe derivative of the homotopy with respect to t isfunction dt!(u, H::RotationAndStraightLine{T}, x::AbstractVector, t, cfg::PolynomialHomotopyConfig, precomputed=false) where {T<:Number}\n    y = H.U(t) * x\n    evaluate_and_jacobian_start_target!(cfg, H, y)\n\n    u .= value_start(cfg) .- value_target(cfg) .+ ((one(t) - t) .* jacobian_target(cfg) .+ t .* jacobian_start(cfg)) * H.U_dot(t) * x\nend\n\nfunction dt!(r::DtDiffResult, H::RotationAndStraightLine{T}, x::AbstractVector, t, cfg::PolynomialHomotopyConfig, precomputed=false) where {T<:Number}\n    y = H.U(t) * x\n    evaluate_and_jacobian_start_target!(cfg, H, y)\n    r.value .= (one(T) - t) .* value_target(cfg) .+ t .* value_start(cfg)\n    r.dt .= value_start(cfg) .- value_target(cfg) .+ ((one(t) - t) .* jacobian_target(cfg) .+ t .* jacobian_start(cfg)) * H.U_dot(t) * x\n    r\nendFinally, we adapt the function to compute the Weylnorm. Note that precomposing with unitary matrices preserves the Weyl inner product.function weylnorm(H::RotationAndStraightLine{T})  where {T<:Number}\n    f = FP.homogenize.(H.start)\n    g = FP.homogenize.(H.target)\n    λ_1 = FP.weyldot(f,f)\n    λ_2 = FP.weyldot(f,g)\n    λ_3 = FP.weyldot(g,g)\n\n    function (t)\n        sqrt(abs2(one(T) - t) * λ_1 + 2 * real((one(T) - t) * conj(t) * λ_2) + abs2(t) * λ_3)\n    end\nend"
},

{
    "location": "set_up_homotopy.html#Include-rotation_and_straightline.jl-in-Homotopies.jl.-1",
    "page": "How to set up your own homotopy",
    "title": "Include rotation_and_straightline.jl in Homotopies.jl.",
    "category": "section",
    "text": "To enable julia to recognize our new homotopy, we have to include the following line in the Homotopies.jl fileinclude(\"homotopies/rotation_and_straightline.jl\")Now we are ready to use RotationAndStraightLine as homotopy type:import DynamicPolynomials: @polyvar\nusing HomotopyContinuation\n\n@polyvar x y\n\nf = [x^2 - x*y]\nH = RotationAndStraightLine(f,f)\n\nsolve(H,[0.0, 1.0 + im * 0.0])gives another solution of f. The technique of making loops in the space of polynomials to track zeros to other zeros is called monodromy.."
},

{
    "location": "set_up_homotopy.html#The-complete-code-1",
    "page": "How to set up your own homotopy",
    "title": "The complete code",
    "category": "section",
    "text": "After having completed all of the above tasks, we have the following rotation_and_straightline.jl file:export RotationAndStraightLine\n\n\"\"\"\n    RotationAndStraightLine(start, target)\n\nConstruct the homotopy `t * start( U(t) x ) + (1-t) * target( U(t) x)`,\n\nwhere `U(t)` is a path in the space of orthogonal matrices with `U(0)=U(1)=I`, the identity matrix.\n\n`start` and `target` have to match and to be one of the following\n* `Vector{<:MP.AbstractPolynomial}` where `MP` is [`MultivariatePolynomials`](https://github.com/blegat/MultivariatePolynomials.jl)\n* `MP.AbstractPolynomial`\n* `Vector{<:FP.Polynomial}` where `FP` is [`FixedPolynomials`](https://github.com/saschatimme/FixedPolynomials.jl)\n\n\n    RotationAndStraightLine{T}(start, target)\n\nYou can also force a specific coefficient type `T`.\n\"\"\"\nmutable struct RotationAndStraightLine{T<:Number} <: AbstractPolynomialHomotopy{T}\n    start::Vector{FP.Polynomial{T}}\n    target::Vector{FP.Polynomial{T}}\n    U::Function\n    U_dot::Function\n\n    function RotationAndStraightLine{T}(start::Vector{FP.Polynomial{T}}, target::Vector{FP.Polynomial{T}}) where {T<:Number}\n        @assert length(start) == length(target) \"Expected the same number of polynomials, but got $(length(start)) and $(length(target))\"\n\n\n        s_nvars = maximum(FP.nvariables.(start))\n        @assert all(s_nvars .== FP.nvariables.(start)) \"Not all polynomials of the start system have $(s_nvars) variables.\"\n\n        t_nvars = maximum(FP.nvariables.(target))\n        @assert all(t_nvars .== FP.nvariables.(target)) \"Not all polynomials of the target system have $(t_nvars) variables.\"\n\n        @assert s_nvars == t_nvars \"Expected start and target system to have the same number of variables, but got $(s_nvars) and $(t_nvars).\"\n\n        U = qrfact(randn(s_nvars,s_nvars) + im * randn(s_nvars,s_nvars))[:Q]\n\n        function U_fct(t)\n                (cos(2 * pi * t) - 1) .* U[:,1] * U[:,1]' - sin(2 * pi * t) .* U[:,2] * U[:,1]' + sin(2 * pi * t) .* U[:,1] * U[:,2]' + (cos(2 * pi * t) - 1) .* U[:,2] * U[:,2]' + eye(U)\n        end\n\n        function U_dot(t)\n                2 * pi .* (-sin(2 * pi * t) .* U[:,1] * U[:,1]' - cos(2 * pi * t) .* U[:,2] * U[:,1]' + cos(2 * pi * t) .* U[:,1] * U[:,2]' - sin(2 * pi * t) .* U[:,2] * U[:,2]')\n        end\n\n        new(start, target, U_fct, U_dot)\n    end\n\n    function RotationAndStraightLine{T}(start, target) where {T<:Number}\n        s, t = construct(T, start, target)\n        RotationAndStraightLine{T}(s, t)\n    end\nend\n\n\nfunction RotationAndStraightLine(start, target)\n    T, s, t = construct(start, target)\n    RotationAndStraightLine{T}(s, t)\nend\n\n\nconst RotationAndStraightLine{T} = RotationAndStraightLine{T}\n\n#\n# SHOW\n#\nfunction Base.deepcopy(H::RotationAndStraightLine)\n    RotationAndStraightLine(deepcopy(H.start), deepcopy(H.target))\nend\n\n#\n# PROMOTION AND CONVERSION\n#ß\nBase.promote_rule(::Type{RotationAndStraightLine{T}}, ::Type{RotationAndStraightLine{S}}) where {S<:Number,T<:Number} = RotationAndStraightLine{promote_type(T,S)}\nBase.promote_rule(::Type{RotationAndStraightLine}, ::Type{S}) where {S<:Number} = RotationAndStraightLine{S}\nBase.promote_rule(::Type{RotationAndStraightLine{T}}, ::Type{S}) where {S<:Number,T<:Number} = RotationAndStraightLine{promote_type(T,S)}\nBase.convert(::Type{RotationAndStraightLine{T}}, H::RotationAndStraightLine) where {T} = RotationAndStraightLine{T}(H.start, H.target)\n\n\n#\n# EVALUATION + DIFFERENTATION\n#\nfunction evaluate!(u::AbstractVector, H::RotationAndStraightLine{T}, x::Vector, t::Number) where T\n    y = H.U(t) * x\n    for i = 1:length(H.target)\n        f = H.target[i]\n        g = H.start[i]\n        u[i] = (one(T) - t) * FP.evaluate(f, y) + t * FP.evaluate(g, y)\n    end\n    u\nend\n(H::RotationAndStraightLine)(x,t) = evaluate(H,x,t)\n\n\nfunction evaluate!(u::AbstractVector{T}, H::RotationAndStraightLine, x::Vector, t::Number, cfg::PolynomialHomotopyConfig, precomputed=false) where {T<:Number}\n    y = H.U(t) * x\n    evaluate_start_target!(cfg, H, y, precomputed)\n    u .= (one(t) - t) .* value_target(cfg) .+ t .* value_start(cfg)\nend\n\nfunction jacobian!(u::AbstractMatrix, H::RotationAndStraightLine{T}, x::AbstractVector, t, cfg::PolynomialHomotopyConfig, precomputed=false) where {T<:Number}\n    U = H.U(t)\n    y = U * x\n    jacobian_start_target!(cfg, H, y, precomputed)\n    u .= ((one(t) - t) .* jacobian_target(cfg) .+ t .* jacobian_start(cfg)) * U\nend\n\nfunction jacobian!(r::JacobianDiffResult, H::RotationAndStraightLine{T}, x::AbstractVector, t, cfg::PolynomialHomotopyConfig, precomputed=false) where {T<:Number}\n    U = H.U(t)\n    y = U * x\n    evaluate_and_jacobian_start_target!(cfg, H, y)\n\n    r.value .= (one(t) - t) .* value_target(cfg) .+ t .* value_start(cfg)\n    r.jacobian .= ((one(t) - t) .* jacobian_target(cfg) .+ t .* jacobian_start(cfg)) * U\n    r\nend\n\nfunction dt!(u, H::RotationAndStraightLine{T}, x::AbstractVector, t, cfg::PolynomialHomotopyConfig, precomputed=false) where {T<:Number}\n    y = H.U(t) * x\n    evaluate_and_jacobian_start_target!(cfg, H, y)\n\n    u .= value_start(cfg) .- value_target(cfg) .+ ((one(t) - t) .* jacobian_target(cfg) .+ t .* jacobian_start(cfg)) * H.U_dot(t) * x\nend\n\nfunction dt!(r::DtDiffResult, H::RotationAndStraightLine{T}, x::AbstractVector, t, cfg::PolynomialHomotopyConfig, precomputed=false) where {T<:Number}\n    y = H.U(t) * x\n    evaluate_and_jacobian_start_target!(cfg, H, y)\n    r.value .= (one(T) - t) .* value_target(cfg) .+ t .* value_start(cfg)\n    r.dt .= value_start(cfg) .- value_target(cfg) .+ ((one(t) - t) .* jacobian_target(cfg) .+ t .* jacobian_start(cfg)) * H.U_dot(t) * x\n    r\nend\n\nfunction weylnorm(H::RotationAndStraightLine{T})  where {T<:Number}\n    f = FP.homogenize.(H.start)\n    g = FP.homogenize.(H.target)\n    λ_1 = FP.weyldot(f,f)\n    λ_2 = FP.weyldot(f,g)\n    λ_3 = FP.weyldot(g,g)\n\n    function (t)\n        sqrt(abs2(one(T) - t) * λ_1 + 2 * real((one(T) - t) * conj(t) * λ_2) + abs2(t) * λ_3)\n    end\nend"
},

{
    "location": "set_up_pathtracker.html#",
    "page": "How to set up your own pathtracking algorithm",
    "title": "How to set up your own pathtracking algorithm",
    "category": "page",
    "text": ""
},

{
    "location": "set_up_pathtracker.html#How-to-set-up-your-own-pathtracking-algorithm-1",
    "page": "How to set up your own pathtracking algorithm",
    "title": "How to set up your own pathtracking algorithm",
    "category": "section",
    "text": "We want to illustrate how to setup your own pathtracking algorithm by the already implemented AffinePredictorCorrector.First you have to define a subtype of AbstractPathtrackingAlgorithm. This is the user facing part.struct AffinePredictorCorrector <: AbstractPathtrackingAlgorithm\nendNote that you could also allow the user to set certain options, e.g. maybe you want to give him/her the choice between an explicit Euler method and a Runge-Kutta method.You also have to clarify whether the algorithm will work in the projective or affine space. Here we want to work in affine space.is_projective(::AffinePredictorCorrector) = falseNow you have to define a struct which is subtype of AbstractPathtrackerCache. This is used for the internal dispatch and also serves as an cache to avoid memory allocations. We will need a working matrix and a vector. Thus we define the followingstruct AffineCache{T} <: AbstractPathtrackerCache{T}\n    A::Matrix{T}\n    b::Vector{T}\nendThen you have to define a new method for alg_cache(algorithm, homotopy, x) which will create our AffineCache:function alg_cache(alg::AffinePredictorCorrector, H::AbstractHomotopy, x::AbstractVector{T}) where T\n    n = length(x)\n    A = zeros(T, n, n)\n    b = zeros(T, n)\n    AffineCache(A, b)\nendWe are already half way done! Now comes the interesting part. We have to define two methods. The first one is a correction method. For us this is a simple newton iteration.function correct!(x, # the startvalue\n    t, # current 'time'\n    H, # the homotopy itself\n    cfg, # An AbstractHomotopyConfiguration for efficient evaluation\n    abstol::Float64, # the target accuracy\n    maxiters::Int, # the maximal number of iterations\n    cache::AffineCache{Complex{T}} # our defined Cache\n    ) where T\n    @unpack A, b = cache\n    m = size(A,2)\n    k = 0\n    while true\n        k += 1\n        evaluate!(b, H, x, t, cfg)\n\n        if norm(b, Inf) < abstol\n            return true\n        elseif k > maxiters\n            return false\n        end\n\n        # put jacobian in A\n        jacobian!(A, H, x, t, cfg, true)\n\n        # this computes A x = b and stores the result x in b\n        LU = lufact!(A)\n        # there is a bug in v0.6.0 see patches.jl\n        my_A_ldiv_B!(LU, b)\n        x .= x .- b\n    end\nendThis method will be used for the refinement of the final solutions. But we can also use it for the next method. We now want to define the method which will actually be used during the pathtracking! For this we have to define the method perform_step!(pathtracker, values, cache::AffineCache). The Pathtracker will invoke this function at each iteration.function perform_step!(tracker, values::PathtrackerPrecisionValues{T}, cache::AffineCache{Complex{T}}) where T\n    @unpack s, ds = tracker # s is our current 'time', ds the step length\n    @unpack H, cfg, x, xnext = values\n    @unpack A, b = cache\n\n    m = size(A,2)\n\n    # PREDICT\n    # put jacobian in A\n    jacobian!(A, H, x, s, cfg)\n    # put Hdt in b\n    dt!(b, H, x, s, cfg, true)\n\n    # this computes A x = b and stores the result x in b\n    LU = lufact!(A)\n    # there is a bug in v0.6.0 see patches.jl\n    my_A_ldiv_B!(LU, b)\n\n    xnext .= x .- ds .* b\n\n    # CORRECT\n    @unpack abstol, corrector_maxiters = tracker.options\n    tracker.step_sucessfull = correct!(xnext, s + ds, H, cfg, abstol, corrector_maxiters, cache)\n    nothing\nendWith this in place you are ready to go! Now you can simply solve a system using your own pathracking algorithm, e.g. using solve(F, AffinePredictorCorrector())."
},

]}
