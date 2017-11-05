####
# Currently the generic (non-BLAS) A_ldiv_B! method with LU factorization does not update
# the RHS. It is already fixed on master (https://github.com/JuliaLang/julia/pull/22774)
# but not in v0.6.1.
# This are the revelenvt parts of the PR from the fix

_apply_ipiv!(A::Base.LinAlg.LU, B::StridedVecOrMat) = _ipiv!(A, 1 : length(A.ipiv), B)

function _ipiv!(A::Base.LinAlg.LU, order::OrdinalRange, B::StridedVecOrMat)
    for i = order
        if i != A.ipiv[i]
            _swap_rows!(B, i, A.ipiv[i])
        end
    end
    B
end
function _swap_rows!(B::StridedVector, i::Integer, j::Integer)
    B[i], B[j] = B[j], B[i]
    B
end
function _swap_rows!(B::StridedMatrix, i::Integer, j::Integer)
    for col = 1 : size(B, 2)
        B[i,col], B[j,col] = B[j,col], B[i,col]
    end
    B
end

# just use the normal one
@inline function my_A_ldiv_B!(A::Base.LinAlg.LU{<:Any,<:StridedMatrix}, B::StridedVecOrMat{T}) where {T<:Base.LinAlg.BlasComplex}
    A_ldiv_B!(A, B)
end
function my_A_ldiv_B!(A::Base.LinAlg.LU{<:Any,<:StridedMatrix}, B::StridedVecOrMat)
    _apply_ipiv!(A, B)
    A_ldiv_B!(Base.LinAlg.UpperTriangular(A.factors), A_ldiv_B!(Base.LinAlg.UnitLowerTriangular(A.factors), B))
end
