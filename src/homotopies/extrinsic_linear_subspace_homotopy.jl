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
    stack(F, StraightLineHomotopy(L₁, L₂))
end


start_parameters!(H::ExtrinsicLinearSubspaceHomotopy, p::LinearSubspace) =
    set_linear_subspace!(homotopy(H.H).start, p)
target_parameters!(H::ExtrinsicLinearSubspaceHomotopy, p::LinearSubspace) =
    set_linear_subspace!(homotopy(H.H).target, p)
function parameters!(
    H::ExtrinsicLinearSubspaceHomotopy,
    p::LinearSubspace,
    q::LinearSubspace,
)
    set_linear_subspace!(homotopy(H.H).start, p)
    set_linear_subspace!(homotopy(H.H).target, q)
end


ModelKit.evaluate!(u, H::ExtrinsicLinearSubspaceHomotopy, x, t) = evaluate!(u, H.H, x, t)
ModelKit.evaluate_and_jacobian!(u, U, H::ExtrinsicLinearSubspaceHomotopy, x, t) =
    evaluate_and_jacobian!(u, U, H.H, x, t)
ModelKit.taylor!(u, v, H::ExtrinsicLinearSubspaceHomotopy, tx, tṫ) =
    taylor!(u, v, H.H, tx, tṫ)


