export IntrinsicLinearSubspaceHomotopy

"""
    IntrinsicLinearSubspaceHomotopy

The homotopty H(x,t) = [F(x,t); t * (A_1 * x - b_1) + (1 - t) * (A_0 * x - b_0)] where
(A_1, b_1) and (A_0, b_0) are the extrinsic subspaces of the start and target subspaces
"""
mutable struct IntrinsicLinearSubspaceHomotopy{T<:AbstractSystem} <: AbstractHomotopy
    H::SystemHomotopyComposition{T,StraightLineHomotopy{LinearSystem,LinearSystem}}
end

function IntrinsicLinearSubspaceHomotopy(
    F::System,
    start,
    target;
    compile::Union{Bool,Symbol} = COMPILE_DEFAULT[],
)
    IntrinsicLinearSubspaceHomotopy(fixed(F; compile = compile), start, target)
end

function IntrinsicLinearSubspaceHomotopy(
    F::AbstractSystem,
    start::LinearSubspace,
    target::LinearSubspace;
    kwargs...,
)

    L₁ = LinearSystem(start, Intrinsic)
    L₂ = LinearSystem(target, Intrinsic)
    H = compose(F, StraightLineHomotopy(L₁, L₂))
    IntrinsicLinearSubspaceHomotopy(H)
end

Base.size(H::IntrinsicLinearSubspaceHomotopy) = size(H.H)


# start_parameters!(H::IntrinsicLinearSubspaceHomotopy, start) = parameters!(H, start, H.target)
# target_parameters!(H::IntrinsicLinearSubspaceHomotopy, target) = parameters!(H, H.start, target)
# parameters!(H::IntrinsicLinearSubspaceHomotopy, p, q) = set_subspaces!(H, p, q)

# function set_subspaces!(
#     H::IntrinsicLinearSubspaceHomotopy,
#     start::LinearSubspace,
#     target::LinearSubspace,
# )
#     H.start = start
#     H.target = target
#     H.Ȧ .= extrinsic(start).A .- extrinsic(target).A
#     H.ḃ .= extrinsic(start).b .- extrinsic(target).b
#     H.t_cache[] = NaN
#     H
# end


function ModelKit.evaluate!(u, H::IntrinsicLinearSubspaceHomotopy, x, t)
    evaluate!(u, H.H, x, t)
end

function ModelKit.evaluate_and_jacobian!(u, U, H::IntrinsicLinearSubspaceHomotopy, x, t)
    evaluate_and_jacobian!(u, U, H.H, x, t)
    nothing
end

function ModelKit.taylor!(u, v, H::IntrinsicLinearSubspaceHomotopy, tx, tṫ)
    taylor!(u, v, H.H, tx, tṫ)
    u
end
