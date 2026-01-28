export ExtrinsicSubspaceHomotopy, IntrinsicSubspaceHomotopy, IntrinsicSubspaceProjectiveHomotopy, set_subspaces!


## Implementation details

# The computation is performed using Stiefel coordinates for the Grassmanian. There are three homotopies moving a linear space V to W:
# * ExtrinsicSubspaceHomotopy
#       for V = {x | Ax - a = 0} and W = {x | Bx - b = 0}
#       moves [F; Ax - a] to [F; Bx - b] by taking a geodesic in the 
#       Grassmannian from A to B and interpolating a and b linearly
# * IntrinsicSubspaceHomotopy
#       for V = {x = Av + a} and W = {x = Bx + b}
#       moves F(Av + a) to F(Bx + b) by taking a geodesic in the 
#       Grassmannian from A to B and interpolating a and b linearly
# * IntrinsicSubspaceProjectiveHomotopy
#       for V = {x = Av + a} and W = {x = Bx + b}
#       moves F(Av + a) to F(Bx + b) by taking a geodesic in the 
#       Grassmannian from colspan[A a] to colspan[B b] 
#       (i.e., it embeds V and W in projective space)


## Data structure for geodesics in the Grassmannian

struct GrassmannianGeodesic
    Q::Matrix{ComplexF64}
    Q_cos::Matrix{ComplexF64}
    Θ::Vector{Float64}
    U::Matrix{ComplexF64}
    γ1::Matrix{ComplexF64}
end

function GrassmannianGeodesic(start, target; 
                             embedded_projective::Bool = false)
    if isa(start, ExtrinsicDescription)
        Q, Θ, U = grassmannian_svd(target, start)
        Q_cos = transpose(target.A) * U
    elseif isa(start, IntrinsicDescription)    
        Q, Θ, U = grassmannian_svd(target, start; embedded_projective = embedded_projective) 
        if !embedded_projective 
            Q_cos = target.X * U
        else
            Q_cos = target.Y * U
        end
    end                  
    
    
    γ1 = similar(Q_cos)
    n, k = size(γ1)
    for j = 1:k
        s, c = sincos(Θ[j])
        for i = 1:n
            γ1[i, j] = Q_cos[i, j] * c + Q[i, j] * s
        end
    end
    G = GrassmannianGeodesic(Q, Q_cos, Θ, U, γ1)
    G
end

const AFFINE_EXTRINSIC_LRU = LRU{
    Tuple{LinearSubspace{ComplexF64},LinearSubspace{ComplexF64}},
    GrassmannianGeodesic,
}(
    maxsize = 128,
)

const AFFINE_INTRINSIC_LRU = LRU{
    Tuple{LinearSubspace{ComplexF64},LinearSubspace{ComplexF64}},
    GrassmannianGeodesic,
}(
    maxsize = 128,
)

const PROJECTIVE_INTRINSIC_LRU = LRU{
    Tuple{LinearSubspace{ComplexF64},LinearSubspace{ComplexF64}},
    GrassmannianGeodesic,
}(
    maxsize = 128,
)


"""
    ExtrinsicSubspaceHomotopy(F::System, V::LinearSubspace, W::LinearSubspace)
    ExtrinsicSubspaceHomotopy(F::AbstractSystem, V::LinearSubspace, W::LinearSubspace)

Creates a homotopy ``H(x,t) = [F(x); A(t)x - a(t)]`` from ``V`` to ``W``.
At ``t=1``, we have ``H(v,1) = [F(x); Ax - a]`` where ``V = {x | Ax = a}`` and ``A`` is a Stiefel matrix.
At ``t=0``, we have ``H(v,0) = [F(x); Bx - b]`` where ``W = {x | Bx = b}`` and ``B`` is a Stiefel matrix.
The matrix part ``A(t)`` follows a geodesic in the Grassmannian in Stiefel coordinates from ``A`` to ``B``, while the offset ``b(t)`` interpolates linearly from ``a`` to ``b``. See also [`LinearSubspace`](@ref).
"""
Base.@kwdef mutable struct ExtrinsicSubspaceHomotopy{S<:AbstractSystem} <: AbstractHomotopy
    system::S

    start::LinearSubspace{ComplexF64}
    target::LinearSubspace{ComplexF64}
    
    # For the matrix part (geodesic interpolation)
    path::GrassmannianGeodesic
    a0::Vector{ComplexF64}
    b0::Vector{ComplexF64}

    # For the offset part (linear interpolation)
    a_minus_b::Union{Nothing, Vector{ComplexF64}}
    offset::Union{Nothing, Vector{ComplexF64}}

    # caches for t
    t_cache::Base.RefValue{ComplexF64}
    offset_t_cache::Base.RefValue{ComplexF64}

    # cache for Jacobian
    J::Matrix{ComplexF64}

    # cache for DF
    ū::Vector{ComplexDF64}

    # for taylor
    v::Vector{ComplexF64}
    L::Vector{ComplexF64}
    L̇::Vector{ComplexF64}
    taylor_t_cache::Base.RefValue{ComplexF64}
    taylor_γ::NTuple{5,Matrix{ComplexF64}}
    tL⁴::TaylorVector{5,ComplexF64}
    tL³::TaylorVector{4,ComplexF64}
    tL²::TaylorVector{3,ComplexF64}

    # c
    gamma::Union{Nothing, ComplexF64}
end



ExtrinsicSubspaceHomotopy(
    F::ModelKit.System,
    start::LinearSubspace,
    target::LinearSubspace;
    compile::Union{Bool,Symbol} = COMPILE_DEFAULT[],
    kwargs...
) = ExtrinsicSubspaceHomotopy(fixed(F; compile = compile), start, target; kwargs...)

function ExtrinsicSubspaceHomotopy(
    system::AbstractSystem,
    start::LinearSubspace,
    target::LinearSubspace;
    kwargs...
)
    ExtrinsicSubspaceHomotopy(
        system,
        convert(LinearSubspace{ComplexF64}, start),
        convert(LinearSubspace{ComplexF64}, target);
        kwargs...
    )
end

function ExtrinsicSubspaceHomotopy(
    system::AbstractSystem,
    start::LinearSubspace{ComplexF64},
    target::LinearSubspace{ComplexF64};
    gamma::Union{Nothing, ComplexF64} = cis(2 * pi * rand())
)

    if !isnothing(gamma)
        gamma = gamma / abs(gamma) # make sure |gamma| = 1
        # multiply with random complex number to get generic paths
        start = LinearSubspace(gamma .* extrinsic(start).A, gamma .* extrinsic(start).b)
    end

    # Create geodesic path for the matrix part
    path = get!(AFFINE_EXTRINSIC_LRU, (start, target)) do
        GrassmannianGeodesic(extrinsic(start), extrinsic(target))
    end
    Q = path.Q
    # get correct coordinates for A and b in the Stiefel homotopy
    Uʰ_start = transpose(path.γ1) * start.extrinsic.A'
    Uʰ_target = transpose(path.Q_cos) * target.extrinsic.A'
    
    # Prepare offset data for linear interpolation
    a0 = Uʰ_start * extrinsic(start).b
    b0 = Uʰ_target * extrinsic(target).b
    a_minus_b = a0 - b0
    J = zeros(ComplexF64, size(system))
    ū = zeros(ComplexDF64, length(a0))
    
    k = length(a0)
    tL⁴ = TaylorVector{5}(ComplexF64, k)

    ExtrinsicSubspaceHomotopy(
        system = system,
        start = start,
        target = target,
        path = path,
        a0 = a0,
        b0 = b0,
        a_minus_b = a_minus_b,
        offset = copy(b0),
        t_cache = Ref(complex(NaN, NaN)),
        offset_t_cache = Ref(complex(NaN, NaN)),
        J = J,
        ū = ū,
        v = zeros(ComplexF64, size(system, 2)),
        L = zeros(ComplexF64, k),
        L̇ = zeros(ComplexF64, k),
        taylor_t_cache = Ref(complex(NaN, NaN)),
        taylor_γ = tuple((similar(Q) for i = 0:4)...),
        tL⁴ = tL⁴,
        tL³ = TaylorVector{4}(tL⁴),
        tL² = TaylorVector{3}(tL⁴),
        gamma = gamma
    )

end

Base.size(H::ExtrinsicSubspaceHomotopy) = (first(size(H.system)) + size(H.path.γ1, 2), last(size(H.system)))


"""
    IntrinsicSubspaceHomotopy(F::System, V::LinearSubspace, W::LinearSubspace)
    IntrinsicSubspaceHomotopy(F::AbstractSystem, V::LinearSubspace, W::LinearSubspace)

Creates a homotopy ``H(v,t) = F(A(t)v + a(t))`` from ``V`` to ``W``.
At ``t=1``, we have ``H(v,1) = F(Av + a)`` where ``V = {x = Av + a}`` and ``A`` is a Stiefel matrix.
At ``t=0``, we have ``H(v,0) = F(Bv + b)`` where ``W = {x = Bv + b}`` and ``B`` is a Stiefel matrix.
The matrix part ``A(t)`` follows a geodesic in the Grassmannian in Stiefel coordinates from ``A`` to ``B``, while the offset ``b(t)`` interpolates linearly from ``a`` to ``b``. See also [`LinearSubspace`](@ref).
"""
Base.@kwdef mutable struct IntrinsicSubspaceHomotopy{S<:AbstractSystem} <: AbstractHomotopy
    system::S

    start::LinearSubspace{ComplexF64}
    target::LinearSubspace{ComplexF64}
    
    # For the matrix part (geodesic interpolation)
    path::GrassmannianGeodesic

    # For the offset part (linear interpolation)
    a_minus_b::Union{Nothing, Vector{ComplexF64}}
    offset::Union{Nothing, Vector{ComplexF64}}

    # cache for t 
    t_cache::Base.RefValue{ComplexF64}
    offset_t_cache::Base.RefValue{ComplexF64}

    # caches for Jacobian, point x and its derivative ẋ
    J::Matrix{ComplexF64}
    x::Vector{ComplexF64}
    ẋ::Vector{ComplexF64}
    x_high::Vector{ComplexDF64}

    # for taylor
    taylor_t_cache::Base.RefValue{ComplexF64}
    taylor_γ::NTuple{5,Matrix{ComplexF64}}
    v::Vector{ComplexF64}
    tx⁴::TaylorVector{5,ComplexF64}
    tx³::TaylorVector{4,ComplexF64}
    tx²::TaylorVector{3,ComplexF64}

    # c 
    gamma::Union{Nothing, ComplexF64}
end


IntrinsicSubspaceHomotopy(
    F::ModelKit.System,
    start::LinearSubspace,
    target::LinearSubspace;
    compile::Union{Bool,Symbol} = COMPILE_DEFAULT[],
    kwargs...
) = IntrinsicSubspaceHomotopy(fixed(F; compile = compile), start, target; kwargs...)

function IntrinsicSubspaceHomotopy(
    system::AbstractSystem,
    start::LinearSubspace,
    target::LinearSubspace;
    kwargs...
)
    IntrinsicSubspaceHomotopy(
        system,
        convert(LinearSubspace{ComplexF64}, start),
        convert(LinearSubspace{ComplexF64}, target);
        kwargs...        
    )
end

function IntrinsicSubspaceHomotopy(
    system::AbstractSystem,
    start::LinearSubspace{ComplexF64},
    target::LinearSubspace{ComplexF64};
    gamma::Union{Nothing, ComplexF64} = cis(2 * pi * rand())
)

    if !isnothing(gamma)
        gamma = gamma / abs(gamma) # make sure |gamma| = 1
        # multiply with random complex number to get generic paths
        start = LinearSubspace(gamma .* extrinsic(start).A, gamma .* extrinsic(start).b)
    end

    # Create geodesic path for the matrix part
    path = get!(AFFINE_INTRINSIC_LRU, (start, target)) do
        GrassmannianGeodesic(intrinsic(start), intrinsic(target))
    end
    Q = path.Q
    
    # Prepare offset data for linear interpolation
    a = intrinsic(start).b
    b = intrinsic(target).b
    a_minus_b = a - b
    offset = copy(b)
    J = zeros(ComplexF64, size(system))
    
    tx⁴ = TaylorVector{5}(ComplexF64, size(Q, 1))

    IntrinsicSubspaceHomotopy(
        system = system,
        start = start,
        target = target,
        path = path,
        a_minus_b = a_minus_b,
        offset = offset,
        t_cache = Ref(complex(NaN, NaN)),
        offset_t_cache = Ref(complex(NaN, NaN)),
        J = J,
        x = zeros(ComplexF64, size(Q, 1)),
        ẋ = zeros(ComplexF64, size(Q, 1)),
        x_high = zeros(ComplexDF64, size(Q, 1)),
        taylor_t_cache = Ref(complex(NaN, NaN)),
        taylor_γ = tuple((similar(Q) for i = 0:4)...),
        v = zeros(ComplexF64, size(Q, 2)),
        tx⁴ = tx⁴,
        tx³ = TaylorVector{4}(tx⁴),
        tx² = TaylorVector{3}(tx⁴),
        gamma = gamma
    )
end

Base.size(H::IntrinsicSubspaceHomotopy) = (size(H.system)[1], dim(H.start))


"""
    IntrinsicSubspaceProjectiveHomotopy(F::System, V::LinearSubspace, W::LinearSubspace)
    IntrinsicSubspaceProjectiveHomotopy(F::AbstractSystem, V::LinearSubspace, W::LinearSubspace)

Creates a homotopy ``H(x,t) = F(γ(t)x)`` where ``γ(t)`` is a family of linear subspaces
such that ``γ(1) = V`` and ``γ(0) = W``.
Here ``γ(t)`` is the geodesic between ``V`` and ``W`` in the Grassmanian after embedding ``V`` and ``W`` in projective space via ``x -> [x; 1]``.
See also [`LinearSubspace`](@ref) and [`geodesic`](@ref) and the references therein.
"""
Base.@kwdef mutable struct IntrinsicSubspaceProjectiveHomotopy{S<:AbstractSystem} <: AbstractHomotopy
    system::S

    start::LinearSubspace{ComplexF64}
    target::LinearSubspace{ComplexF64}
    
    # For the matrix part (geodesic interpolation)
    path::GrassmannianGeodesic

    t_cache::Base.RefValue{ComplexF64}
    offset_t_cache::Base.RefValue{ComplexF64}

    J::Matrix{ComplexF64}
    x::Vector{ComplexF64}
    ẋ::Vector{ComplexF64}
    x_high::Vector{ComplexDF64}

    # For AD
    taylor_t_cache::Base.RefValue{ComplexF64}
    taylor_γ::NTuple{5,Matrix{ComplexF64}}
    v::Vector{ComplexF64}
    tx⁴::TaylorVector{5,ComplexF64}
    tx³::TaylorVector{4,ComplexF64}
    tx²::TaylorVector{3,ComplexF64}
    tx¹::TaylorVector{2,ComplexF64}

    # c 
    gamma::Union{Nothing, ComplexF64}
end


IntrinsicSubspaceProjectiveHomotopy(
    F::ModelKit.System,
    start::LinearSubspace,
    target::LinearSubspace;
    compile::Union{Bool,Symbol} = COMPILE_DEFAULT[],
    kwargs...
) = IntrinsicSubspaceProjectiveHomotopy(fixed(F; compile = compile), start, target; kwargs...)

function IntrinsicSubspaceProjectiveHomotopy(
    system::AbstractSystem,
    start::LinearSubspace,
    target::LinearSubspace;
    kwargs...
)
    IntrinsicSubspaceProjectiveHomotopy(
        system,
        convert(LinearSubspace{ComplexF64}, start),
        convert(LinearSubspace{ComplexF64}, target);
        kwargs...
    )
end

function IntrinsicSubspaceProjectiveHomotopy(
    system::AbstractSystem,
    start::LinearSubspace{ComplexF64},
    target::LinearSubspace{ComplexF64};
    gamma::Union{Nothing, ComplexF64} = cis(2 * pi * rand())
)

    if !isnothing(gamma)
        gamma = gamma / abs(gamma) # make sure |gamma| = 1
        # multiply with random complex number to get generic paths
        start = LinearSubspace(gamma .* extrinsic(start).A, gamma .* extrinsic(start).b)
    end

    # Create geodesic path for the matrix part
    path = get!(PROJECTIVE_INTRINSIC_LRU, (start, target)) do
        GrassmannianGeodesic(intrinsic(start), intrinsic(target); embedded_projective = true)
    end
    Q = path.Q
 
    J = zeros(ComplexF64, size(system) .+ (1, 1))
    tx⁴ = TaylorVector{5}(ComplexF64, size(Q, 1))

    IntrinsicSubspaceProjectiveHomotopy(
        system = system,
        start = start,
        target = target,
        path = path,
        t_cache = Ref(complex(NaN, NaN)),
        offset_t_cache = Ref(complex(NaN, NaN)),
        J = J,
        x = zeros(ComplexF64, size(Q, 1)),
        ẋ = zeros(ComplexF64, size(Q, 1)),
        x_high = zeros(ComplexDF64, size(Q, 1)),
        taylor_t_cache = Ref(complex(NaN, NaN)),
        taylor_γ = tuple((similar(Q) for i = 0:4)...),
        v = zeros(ComplexF64, size(Q, 2)),
        tx⁴ = tx⁴,
        tx³ = TaylorVector{4}(tx⁴),
        tx² = TaylorVector{3}(tx⁴),
        tx¹ = TaylorVector{2}(tx⁴),
        gamma = gamma
    )
end

Base.size(H::IntrinsicSubspaceProjectiveHomotopy) = (size(H.system)[1] + 1, dim(H.start) + 1)



const SubspaceHomotopy = Union{ExtrinsicSubspaceHomotopy, IntrinsicSubspaceHomotopy, IntrinsicSubspaceProjectiveHomotopy}

"""
    set_subspaces!(H::SubspaceHomotopy, start::LinearSubspace, target::LinearSubspace)

Update the homotopy `H` to track from the linear subspace `start` to `target`, where `H`. 
"""
function set_subspaces!(
    H::SubspaceHomotopy,
    start::LinearSubspace,
    target::LinearSubspace,
)
    H.start = start
    H.target = target

    if isa(H, ExtrinsicSubspaceHomotopy)
        LRU = AFFINE_EXTRINSIC_LRU
    elseif isa(H, IntrinsicSubspaceHomotopy)
        LRU = AFFINE_INTRINSIC_LRU
    elseif isa(H, IntrinsicSubspaceProjectiveHomotopy)
        LRU = PROJECTIVE_INTRINSIC_LRU
    end

    if isa(H, ExtrinsicSubspaceHomotopy)
        H.path = get!(LRU, (start, target)) do
                GrassmannianGeodesic(extrinsic(start), intrinsic(target))
            end
    else
        H.path = get!(LRU, (start, target)) do
                GrassmannianGeodesic(intrinsic(start), intrinsic(target))
            end
    end

    # Update the offsets
    if isa(H, ExtrinsicSubspaceHomotopy) 
        a = H.a0
        b = H.b0
        H.a_minus_b .= a .- b
        H.offset .= b
    elseif isa(H, IntrinsicSubspaceHomotopy)
        a = intrinsic(start).b
        b = intrinsic(target).b
        H.a_minus_b .= a .- b
        H.offset .= b
    end

    # set caches
    H.t_cache[] = complex(0.0)
    H.offset_t_cache[] = complex(0.0)
    H.taylor_t_cache[] = complex(0.0)
    
    H
end

start_parameters!(H::SubspaceHomotopy, p::LinearSubspace) =
    set_subspaces!(H, convert(LinearSubspace{ComplexF64}, p), H.target)
target_parameters!(H::SubspaceHomotopy, q::LinearSubspace) =
    set_subspaces!(H, H.start, convert(LinearSubspace{ComplexF64}, q))
parameters!(H::SubspaceHomotopy, p::LinearSubspace, q::LinearSubspace) =
    set_subspaces!(
        H,
        convert(LinearSubspace{ComplexF64}, p),
        convert(LinearSubspace{ComplexF64}, q),
    )

# Helper function for computing the geodesic γ(t) at parameter t
function γ!(H::SubspaceHomotopy, t::Number)
    H.t_cache[] != t || return first(H.taylor_γ)
    if isreal(t)
        _γ!(H, real(t))
    else
        _γ!(H, t)
    end
    H.t_cache[] = t

    first(H.taylor_γ)
end

@inline function _γ!(H::SubspaceHomotopy, t::Number)
    @unpack Q, Q_cos, Θ = H.path
    γ = first(H.taylor_γ)
    n, k = size(γ)
    @inbounds for j = 1:k
        Θⱼ = Θ[j]
        s, c = sincos(t * Θⱼ)
        for i = 1:n
            γ[i, j] = Q_cos[i, j] * c + Q[i, j] * s
        end
    end
    γ
end

γ̇!(H::SubspaceHomotopy, t::Number) = isreal(t) ? _γ̇!(H, real(t)) : _γ̇!(H, t)

@inline function _γ̇!(H::SubspaceHomotopy, t::Number)
    @unpack Q, Q_cos, Θ = H.path
    _, γ̇ = H.taylor_γ
    n, k = size(γ̇)
    @inbounds for j = 1:k
        Θⱼ = Θ[j]
        s, c = sincos(t * Θⱼ)
        ċ = -s * Θⱼ
        ṡ = c * Θⱼ
        for i = 1:n
            γ̇[i, j] = Q_cos[i, j] * ċ + Q[i, j] * ṡ
        end
    end
    γ̇
end


# Compute the offset at parameter t: offset(t) = t*a + (1-t)*b = t(a-b) + b
@inline function _compute_offset!(H::SubspaceHomotopy, t::Number)
    if isa(H, ExtrinsicSubspaceHomotopy)
        H.offset .= H.b0
    elseif isa(H, IntrinsicSubspaceHomotopy)
        H.offset .= intrinsic(H.target).b
    end
    LA.axpy!(t, H.a_minus_b, H.offset)
    H.offset_t_cache[] = t
end

@inline function offset_at_t!(H::SubspaceHomotopy, t::Number)
    if H.offset_t_cache[] != t
        _compute_offset!(H, t)
    end
    H.offset
end


# Helper functions to compute taylor series
function _taylor_γ!(H::SubspaceHomotopy, t::Number)
    @unpack path, taylor_γ = H
    @unpack Q, Q_cos, Θ, U = path

    γ, γ¹, γ², γ³, γ⁴ = taylor_γ
    n, k = size(γ)
    @inbounds for j = 1:k
        Θⱼ = Θ[j]
        s, c = sincos(t * Θⱼ)
        c¹ = -s * Θⱼ
        s¹ = c * Θⱼ
        Θⱼ_2 = 0.5 * Θⱼ^2
        c² = -c * Θⱼ_2
        s² = -s * Θⱼ_2
        Θⱼ_3 = Θⱼ_2 * Θⱼ / 3
        c³ = s * Θⱼ_3
        s³ = -c * Θⱼ_3
        Θⱼ_4 = 0.25 * Θⱼ_3 * Θⱼ
        c⁴ = c * Θⱼ_4
        s⁴ = s * Θⱼ_4
        for i = 1:n
            γ[i, j] = Q_cos[i, j] * c + Q[i, j] * s
            γ¹[i, j] = Q_cos[i, j] * c¹ + Q[i, j] * s¹
            γ²[i, j] = Q_cos[i, j] * c² + Q[i, j] * s²
            γ³[i, j] = Q_cos[i, j] * c³ + Q[i, j] * s³
            γ⁴[i, j] = Q_cos[i, j] * c⁴ + Q[i, j] * s⁴
        end
    end
end

function taylor_γ!(H::SubspaceHomotopy, t::Number)
    H.taylor_t_cache[] != t || return H.taylor_γ

    if isreal(t)
        _taylor_γ!(H, real(t))
    else
        _taylor_γ!(H, t)
    end
    H.taylor_t_cache[] = t

    H.taylor_γ
end


############################
## ExtrinsicSubspaceHomotopy
############################ 

function ModelKit.evaluate!(u, H::ExtrinsicSubspaceHomotopy, x::AbstractVector, t)
    γ = γ!(H, t)
    offset_at_t!(H, t)
    n, k = size(γ) # linear equation is: transpose(γ) x - offset

    evaluate!(u, H.system, x)
    m = first(size(H.system))
    for i = 1:k
        u[m+i] = -H.offset[i]
    end
    for j = 1:n
        xj = x[j]
        for i = 1:k
            u[m+i] = muladd(γ[j, i], xj, u[m+i])
        end
    end
    u
end
function ModelKit.evaluate!(u, H::ExtrinsicSubspaceHomotopy, x::Vector{ComplexDF64}, t)
    γ = γ!(H, t)
    offset_at_t!(H, t)
    n, k = size(γ) # linear equation is: tranpose(γ) x - offset

    evaluate!(u, H.system, x)
    m = first(size(H.system))
    for i = 1:k
        H.ū[i] = -H.offset[i]
    end
    for j = 1:n
        xj = x[j]
        for i = 1:k
            H.ū[i] = muladd(γ[j, i], xj, H.ū[i])
        end
    end
    for i = 1:k
        u[m+i] = H.ū[i]
    end

    u
end

function ModelKit.evaluate_and_jacobian!(u, U, H::ExtrinsicSubspaceHomotopy, x::AbstractVector, t)
    γ = γ!(H, t)
    offset_at_t!(H, t)
    n, k = size(γ) # linear equation is: tranpose(γ) x - offset
    
    evaluate_and_jacobian!(u, U, H.system, x)
    m = first(size(H.system))
    for i = 1:k
        u[m+i] = -H.offset[i]
    end
    for j = 1:n
        xj = x[j]
        for i = 1:k
            u[m+i] = muladd(γ[j, i], xj, u[m+i])
        end
    end
    for j = 1:n, i = 1:k
        U[m+i, j] = γ[j, i]
    end

    nothing
end


function ModelKit.taylor!(
    u,
    ::Val{1},
    H::ExtrinsicSubspaceHomotopy,
    v,
    t,
    incr::Bool = true,
)
    # d/dt H(x,t) 
    # = d/dt [F(x); L(t)], where L(t) = A(t) x - b(t)
    # = [0; L̇(t)]  (since x is independent of t)

    γ̇ = γ̇!(H, t)

    m = first(size(H.system))
    k = size(γ̇, 2) # recall that A(t) = transpose(γ(t)), so A(t) is kxn

    H.v .= first.(v)

    LA.mul!(H.L̇, transpose(γ̇), H.v)
    offset_at_t!(H, t)
    LA.axpy!(-1, H.a_minus_b, H.L̇)  # -H.a_minus_b = ḃ = derivative of offset with respect to t 

    # F part is zero
    for i = 1:m
        u[i] = zero(eltype(u))
    end
    # rest is H.L̇
    for i = 1:k
        u[m+i] = H.L̇[i]
    end
    
    u
end

function ModelKit.taylor!(
    u,
    ::Val{2},
    H::ExtrinsicSubspaceHomotopy,
    tv::TaylorVector,
    t,
    incr::Bool = false,
)
    γ, γ¹, γ² = taylor_γ!(H, t)
    L, L¹, L² = vectors(H.tL²)
    v, v¹ = vectors(tv)

    if incr
        L .= H.L
        L¹ .= H.L̇
    else
        LA.mul!(L, transpose(γ), v)
        LA.mul!(L¹, transpose(γ¹), v)
        offset_at_t!(H, t)
        LA.axpy!(-1, H.offset, L)
        LA.axpy!(-1, H.a_minus_b, L¹) 
    end

    H.v .= v¹
    LA.mul!(H.L, transpose(γ), v¹)
    L¹ .+= H.L

    LA.mul!(H.L, transpose(γ¹), H.v)
    H.v .= v
    LA.mul!(H.L, transpose(γ²), H.v, true, true)
    L² .= H.L

    m = first(size(H.system))
    k = size(γ, 2) # recall that A(t) = transpose(γ(t)), so A(t) is kxn

    # F part is zero
    for i = 1:m
        u[i] = zero(eltype(u))
    end
    # rest is L²
    for i = 1:k
        u[m+i] = L²[i]
    end
    
    u
end

function ModelKit.taylor!(
    u,
    ::Val{3},
    H::ExtrinsicSubspaceHomotopy,
    tv::TaylorVector,
    t,
    incr::Bool = false,
)
    γ, γ¹, γ², γ³ = taylor_γ!(H, t)
    L, L¹, L², L³ = vectors(H.tL³)
    v, v¹, v² = vectors(tv)

    if !incr
        LA.mul!(L, transpose(γ), v)
        LA.mul!(L¹, transpose(γ¹), v)
        offset_at_t!(H, t)
        LA.axpy!(-1, H.offset, L)
        LA.axpy!(-1, H.a_minus_b, L¹)
        LA.mul!(L¹, transpose(γ), v¹, true, true)
        LA.mul!(L², transpose(γ²), v)
        LA.mul!(L², transpose(γ¹), v¹, true, true)
    end

    H.v .= v²
    LA.mul!(H.L, transpose(γ), H.v)
    L² .+= H.L

    LA.mul!(H.L, transpose(γ¹), H.v)
    H.v .= v
    LA.mul!(H.L, transpose(γ³), H.v, true, true)
    H.v .= v¹
    LA.mul!(H.L, transpose(γ²), H.v, true, true)
    L³ .= H.L
    
    m = first(size(H.system))
    k = size(γ, 2) # recall that A(t) = transpose(γ(t)), so A(t) is kxn

    # F part is zero
    for i = 1:m
        u[i] = zero(eltype(u))
    end
    # rest is L³
    for i = 1:k
        u[m+i] = L³[i]
    end

    u
end

function ModelKit.taylor!(
    u,
    ::Val{4},
    H::ExtrinsicSubspaceHomotopy,
    tv::TaylorVector,
    t,
    incr::Bool = false,
)
    γ, γ¹, γ², γ³, γ⁴ = taylor_γ!(H, t)
    L, L¹, L², L³, L⁴ = vectors(H.tL⁴)
    v, v¹, v², v³ = vectors(tv)

    if !incr
        LA.mul!(L, transpose(γ), v)
        LA.mul!(L¹, transpose(γ¹), v)
        offset_at_t!(H, t)
        LA.axpy!(-1, H.offset, L)
        LA.axpy!(-1, H.a_minus_b, L¹)
        LA.mul!(L¹, transpose(γ), v¹, true, true)
        LA.mul!(L², transpose(γ²), v)
        LA.mul!(L², transpose(γ¹), v¹, true, true)
        LA.mul!(L², transpose(γ), v², true, true)

        LA.mul!(L³, transpose(γ³), v)
        LA.mul!(L³, transpose(γ²), v¹, true, true)
        LA.mul!(L³, transpose(γ¹), v², true, true)
    end

    LA.mul!(L³, transpose(γ), v³, true, true)

    LA.mul!(L⁴, transpose(γ¹), v³)
    LA.mul!(L⁴, transpose(γ²), v², true, true)
    LA.mul!(L⁴, transpose(γ³), v¹, true, true)
    LA.mul!(L⁴, transpose(γ⁴), v, true, true)
    
    m = first(size(H.system))
    k = size(γ, 2) # recall that A(t) = transpose(γ(t)), so A(t) is kxn

    # F part is zero
    for i = 1:m
        u[i] = zero(eltype(u))
    end
    # rest is L⁴
    for i = 1:k
        u[m+i] = L⁴[i]
    end

    u
end


############################
## IntrinsicSubspaceHomotopy
############################ 

function set_solution!(u::Vector, H::IntrinsicSubspaceHomotopy, x::AbstractVector, t)
    (length(x) == length(H.x)) ||
    throw(ArgumentError("Cannot set solution. Expected extrinsic coordinates."))

    set_solution!(H.x, H.system, x)
    offset_at_t!(H, t)
    H.x .= H.x .- H.offset

    if isone(t)
        LA.mul!(u, H.path.γ1', H.x)
    elseif iszero(t)
        LA.mul!(u, H.path.Q_cos', H.x)
    else
        LA.mul!(u, γ!(H, t)', H.x)
    end
end

function get_solution(H::IntrinsicSubspaceHomotopy, u::AbstractVector, t)
    offset_at_t!(H, t)
    if isone(t)
        out = (H.path.γ1 * u) + H.offset
    elseif iszero(t)
        out = (H.path.Q_cos * u) + H.offset
    else
        γ = γ!(H, t)
        out = (γ * u) + H.offset
    end
    
    out
end

function ModelKit.evaluate!(u, H::IntrinsicSubspaceHomotopy, v::Vector{ComplexDF64}, t)
    γ = γ!(H, t)
    LA.mul!(H.x_high, γ, v)

    offset_at_t!(H, t)
    LA.axpy!(1, H.offset, H.x_high)

    evaluate!(u, H.system, H.x_high)
   
    u
end

function ModelKit.evaluate!(u, H::IntrinsicSubspaceHomotopy, v::AbstractVector, t)
    γ = γ!(H, t)
    LA.mul!(H.x, γ, v)

    offset_at_t!(H, t)
    LA.axpy!(1, H.offset, H.x)

    evaluate!(u, H.system, H.x)
   
    u
end

function ModelKit.evaluate_and_jacobian!(
    u,
    U,
    H::IntrinsicSubspaceHomotopy,
    v::AbstractVector,
    t,
)
    γ = γ!(H, t)
    LA.mul!(H.x, γ, v)

    offset_at_t!(H, t)
    LA.axpy!(1, H.offset, H.x)  
    
    evaluate_and_jacobian!(u, H.J, H.system, H.x)
    LA.mul!(U, H.J, γ)

    nothing
end

function ModelKit.taylor!(
    u,
    ::Val{1},
    H::IntrinsicSubspaceHomotopy,
    v,
    t,
    incr::Bool = true,
)
    # d/dt H(v,t) = d/dt F(A(t)v + b(t))
    #             = J_F(A(t)v + b(t)) * (Ȧ(t)v + ḃ(t))

    γ = γ!(H, t)
    γ̇ = γ̇!(H, t)

    H.v .= first.(v)
    LA.mul!(H.x, γ, H.v)
    LA.mul!(H.ẋ, γ̇, H.v)

    offset_at_t!(H, t)
    LA.axpy!(1, H.offset, H.x)
    LA.axpy!(1, H.a_minus_b, H.ẋ)  # H.a_minus_b = ḃ = derivative of offset with respect to t 
    
    evaluate_and_jacobian!(u, H.J, H.system, H.x)
    LA.mul!(u, H.J, H.ẋ)

    u
end

function ModelKit.taylor!(
    u,
    ::Val{2},
    H::IntrinsicSubspaceHomotopy,
    tv::TaylorVector,
    t,
    incr::Bool = false,
)
    γ, γ¹, γ² = taylor_γ!(H, t)
    x, x¹, x² = vectors(H.tx²)
    v, v¹ = vectors(tv)

    if incr
        x .= H.x
        x¹ .= H.ẋ
    else
        LA.mul!(x, γ, v)
        LA.mul!(x¹, γ¹, v)
        offset_at_t!(H, t)
        LA.axpy!(1, H.offset, x)
        LA.axpy!(1, H.a_minus_b, x¹) 
    end

    H.v .= v¹
    LA.mul!(H.x, γ, v¹)
    x¹ .+= H.x

    LA.mul!(H.x, γ¹, H.v)
    H.v .= v
    LA.mul!(H.x, γ², H.v, true, true)
    x² .= H.x

    taylor!(u, Val(2), H.system, H.tx²)

    u
end

function ModelKit.taylor!(
    u,
    ::Val{3},
    H::IntrinsicSubspaceHomotopy,
    tv::TaylorVector,
    t,
    incr::Bool = false,
)
    γ, γ¹, γ², γ³ = taylor_γ!(H, t)
    x, x¹, x², x³ = vectors(H.tx³)
    v, v¹, v² = vectors(tv)

    if !incr
        LA.mul!(x, γ, v)
        LA.mul!(x¹, γ¹, v)
        offset_at_t!(H, t)
        LA.axpy!(1, H.offset, x)
        LA.axpy!(1, H.a_minus_b, x¹)
        LA.mul!(x¹, γ, v¹, true, true)
        LA.mul!(x², γ², v)
        LA.mul!(x², γ¹, v¹, true, true)
    end

    H.v .= v²
    LA.mul!(H.x, γ, H.v)
    x² .+= H.x

    LA.mul!(H.x, γ¹, H.v)
    H.v .= v
    LA.mul!(H.x, γ³, H.v, true, true)
    H.v .= v¹
    LA.mul!(H.x, γ², H.v, true, true)
    x³ .= H.x

    taylor!(u, Val(3), H.system, H.tx³)

    u
end

function ModelKit.taylor!(
    u,
    ::Val{4},
    H::IntrinsicSubspaceHomotopy,
    tv::TaylorVector,
    t,
    incr::Bool = false,
)
    γ, γ¹, γ², γ³, γ⁴ = taylor_γ!(H, t)
    x, x¹, x², x³, x⁴ = vectors(H.tx⁴)
    v, v¹, v², v³ = vectors(tv)

    if !incr
        LA.mul!(x, γ, v)
        LA.mul!(x¹, γ¹, v)
        offset_at_t!(H, t)
        LA.axpy!(1, H.offset, x)
        LA.axpy!(1, H.a_minus_b, x¹)
        LA.mul!(x¹, γ, v¹, true, true)
        LA.mul!(x², γ², v)
        LA.mul!(x², γ¹, v¹, true, true)
        LA.mul!(x², γ, v², true, true)

        LA.mul!(x³, γ³, v)
        LA.mul!(x³, γ², v¹, true, true)
        LA.mul!(x³, γ¹, v², true, true)
    end

    LA.mul!(x³, γ, v³, true, true)

    LA.mul!(x⁴, γ¹, v³)
    LA.mul!(x⁴, γ², v², true, true)
    LA.mul!(x⁴, γ³, v¹, true, true)
    LA.mul!(x⁴, γ⁴, v, true, true)

    taylor!(u, Val(4), H.system, H.tx⁴)

    u
end


######################################
## IntrinsicSubspaceProjectiveHomotopy
######################################


function set_solution!(u::Vector, H::IntrinsicSubspaceProjectiveHomotopy, x::AbstractVector, t)

    (length(x) == length(H.x) - 1) ||
        throw(ArgumentError("Cannot set solution. Expected extrinsic coordinates."))

    set_solution!(view(H.x, 1:length(x)), H.system, x)
    H.x[end] = 1

    if isone(t)
        LA.mul!(u, H.path.γ1', H.x)
    elseif iszero(t)
        LA.mul!(u, H.path.Q_cos', H.x)
    else
        LA.mul!(u, γ!(H, t)', H.x)
    end
end

function get_solution(H::IntrinsicSubspaceProjectiveHomotopy, u::AbstractVector, t)
    
    if isone(t)
        out = (@view H.path.γ1[1:end-1, :]) * u
    elseif iszero(t)
        out = (@view H.path.Q_cos[1:end-1, :]) * u
    else
        γ = γ!(H, t)
        out = (@view γ[1:end-1, :]) * u
    end

    out
end

function ModelKit.evaluate!(u, H::IntrinsicSubspaceProjectiveHomotopy, v::Vector{ComplexDF64}, t)
    γ = γ!(H, t)
    LA.mul!(H.x_high, γ, v)
    evaluate!(u, H.system, H.x_high)

    n = first(size(H.system))
    u[n+1] = H.x_high[end] - 1.0
 
    u
end

function ModelKit.evaluate!(u, H::IntrinsicSubspaceProjectiveHomotopy, v::AbstractVector, t)
    γ = γ!(H, t)
    LA.mul!(H.x, γ, v)

    evaluate!(u, H.system, H.x)
   
    n = first(size(H.system))
    u[n+1] = H.x[end] - 1.0

    u
end

function ModelKit.evaluate_and_jacobian!(
    u,
    U,
    H::IntrinsicSubspaceProjectiveHomotopy,
    v::AbstractVector,
    t,
)
    γ = γ!(H, t)
    LA.mul!(H.x, γ, v)
    
    evaluate_and_jacobian!(u, H.J, H.system, H.x)
    LA.mul!(U, H.J, γ)

    n = first(size(H.system))
    u[n+1] = H.x[end] - 1
    m = length(v)
    for j = 1:m
        U[n+1, j] = γ[end, j]
    end

    nothing
end

function ModelKit.taylor!(
    u,
    ::Val{1},
    H::IntrinsicSubspaceProjectiveHomotopy,
    v,
    t,
    incr::Bool = true,
)
    # d/dt H(v,t) = d/dt F(A(t)v)
    #             = J_F(A(t)) * (Ȧ(t)v)

    γ = γ!(H, t)
    γ̇ = γ̇!(H, t)

    H.v .= first.(v)
    LA.mul!(H.x, γ, H.v)
    LA.mul!(H.ẋ, γ̇, H.v)
    
    evaluate_and_jacobian!(u, H.J, H.system, H.x)
    LA.mul!(u, H.J, H.ẋ)
    M = size(H, 1)
    u[M] = H.ẋ[end]

    u
end

function ModelKit.taylor!(
    u,
    ::Val{2},
    H::IntrinsicSubspaceProjectiveHomotopy,
    tv::TaylorVector,
    t,
    incr::Bool = false,
)
    γ, γ¹, γ² = taylor_γ!(H, t)
    x, x¹, x² = vectors(H.tx²)
    v, v¹ = vectors(tv)

    if incr
        x .= H.x
        x¹ .= H.ẋ
    else
        LA.mul!(x, γ, v)
        LA.mul!(x¹, γ¹, v)
    end

    H.v .= v¹
    LA.mul!(H.x, γ, v¹)
    x¹ .+= H.x

    LA.mul!(H.x, γ¹, H.v)
    H.v .= v
    LA.mul!(H.x, γ², H.v, true, true)
    x² .= H.x

    taylor!(u, Val(2), H.system, H.tx²)

    n = first(size(H.system))
    u[n+1] = x²[end]

    u
end

function ModelKit.taylor!(
    u,
    ::Val{3},
    H::IntrinsicSubspaceProjectiveHomotopy,
    tv::TaylorVector,
    t,
    incr::Bool = false,
)
    γ, γ¹, γ², γ³ = taylor_γ!(H, t)
    x, x¹, x², x³ = vectors(H.tx³)
    v, v¹, v² = vectors(tv)

    if !incr
        LA.mul!(x, γ, v)
        LA.mul!(x¹, γ¹, v)
        LA.mul!(x¹, γ, v¹, true, true)
        LA.mul!(x², γ², v)
        LA.mul!(x², γ¹, v¹, true, true)
    end

    H.v .= v²
    LA.mul!(H.x, γ, H.v)
    x² .+= H.x

    LA.mul!(H.x, γ¹, H.v)
    H.v .= v
    LA.mul!(H.x, γ³, H.v, true, true)
    H.v .= v¹
    LA.mul!(H.x, γ², H.v, true, true)
    x³ .= H.x

    taylor!(u, Val(3), H.system, H.tx³)

    n = first(size(H.system))
    u[n+1] = x³[end]

    u
end

function ModelKit.taylor!(
    u,
    ::Val{4},
    H::IntrinsicSubspaceProjectiveHomotopy,
    tv::TaylorVector,
    t,
    incr::Bool = false,
)
    γ, γ¹, γ², γ³, γ⁴ = taylor_γ!(H, t)
    x, x¹, x², x³, x⁴ = vectors(H.tx⁴)
    v, v¹, v², v³ = vectors(tv)

    if !incr
        LA.mul!(x, γ, v)
        LA.mul!(x¹, γ¹, v)
        LA.mul!(x¹, γ, v¹, true, true)
        LA.mul!(x², γ², v)
        LA.mul!(x², γ¹, v¹, true, true)
        LA.mul!(x², γ, v², true, true)

        LA.mul!(x³, γ³, v)
        LA.mul!(x³, γ², v¹, true, true)
        LA.mul!(x³, γ¹, v², true, true)
    end

    LA.mul!(x³, γ, v³, true, true)

    # Currently H.v .= v³
    LA.mul!(x⁴, γ¹, v³)
    LA.mul!(x⁴, γ², v², true, true)
    LA.mul!(x⁴, γ³, v¹, true, true)
    LA.mul!(x⁴, γ⁴, v, true, true)

    taylor!(u, Val(4), H.system, H.tx⁴)

    n = first(size(H.system))
    u[n+1] = x⁴[end]

    u
end
