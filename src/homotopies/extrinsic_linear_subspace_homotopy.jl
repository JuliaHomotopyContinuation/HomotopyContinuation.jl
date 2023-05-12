export ExtrinsicLinearSubspaceHomotopy

"""
    ExtrinsicLinearSubspaceHomotopy

The homotopty H(x,t) = [F(x,t); t * (A_1 * x - b_1) + (1 - t) * (A_0 * x - b_0)] where
(A_1, b_1) and (A_0, b_0) are the extrinsic subspaces of the start and target subspaces
"""
mutable struct ExtrinsicLinearSubspaceHomotopy{T<:AbstractSystem} <: AbstractHomotopy
    H::SystemHomotopyStack{T,StraightLineHomotopy{LinearSystem,LinearSystem}}
end

function ExtrinsicLinearSubspaceHomotopy(
    F::System,
    start,
    target;
    compile::Union{Bool,Symbol} = COMPILE_DEFAULT[],
)
    ExtrinsicLinearSubspaceHomotopy(fixed(F; compile = compile), start, target)
end

function ExtrinsicLinearSubspaceHomotopy(
    F::AbstractSystem,
    start::LinearSubspace,
    target::LinearSubspace;
    kwargs...,
)

    L₁ = LinearSystem(start)
    L₂ = LinearSystem(target)
    H = stack(F, StraightLineHomotopy(L₁, L₂))
end


# start_parameters!(H::ExtrinsicLinearSubspaceHomotopy, start) = parameters!(H, start, H.target)
# target_parameters!(H::ExtrinsicLinearSubspaceHomotopy, target) = parameters!(H, H.start, target)
# parameters!(H::ExtrinsicLinearSubspaceHomotopy, p, q) = set_subspaces!(H, p, q)

# function set_subspaces!(
#     H::ExtrinsicLinearSubspaceHomotopy,
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


function ModelKit.evaluate!(u, H::ExtrinsicLinearSubspaceHomotopy, x, t)
    evaluate!(u, H.H, x, t)
end

function ModelKit.evaluate_and_jacobian!(u, U, H::ExtrinsicLinearSubspaceHomotopy, x, t)
    evaluate_and_jacobian!(u, U, H.H, x, t)
    nothing
end

function ModelKit.taylor!(u, v, H::ExtrinsicLinearSubspaceHomotopy, tx, tṫ)
    taylor!(u, v, H.H, tx, tṫ)
    u
end
