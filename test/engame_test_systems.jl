# Root (-1, 1) multiplicity k root
function eg_system_1(k::Int)
    @var x y
    System([x * y - 1, (x + 1)^k])
end

# multiplicity 4 root
function eg_system_2()
    @var x y
    System([x^2 + y^2, x^2 - y^2])
end

# Irred. decomp contains 4 components
# {x = 0} ∪ { y^2 - x^3 = 0} ∪ {(1,2)} ∪ {(1,-3)}
function eg_system_3()
    @var x y
    System([x * (y^2 - x^3) * (x - 1), x * (y^2 - x^3) * (y - 2) * (3x + y)])
end

# multiplicity 3 root at (0,0)
function eg_system_4()
    @var x y
    System([y^2 - x^3, x + 2y])
end

# 2 solutions with multiplicity 6, projective
function eg_system_5()
    @var x z
    y = 1
    System(
        [
            0.75 * x^4 + 1.5 * x^2 * y^2 - 2.5 * x^2 * z^2 + 0.75 * y^4 - 2.5 * y^2 * z^2 +
            0.75 * z^4
            10 * x^2 * z + 10 * y^2 * z - 6 * z^3
        ],
    )
end

# 1 singular solution with multiplicity 3
function eg_system_6()
    @var x y
    z = 1
    System([
        x^2 + 2 * y^2 + 2 * im * y * z,
        (18 + 3 * im) * x * y + 7 * im * y^2 - (3 - 18 * im) * x * z - 14 * y * z -
        7 * im * z^2,
    ])

end

# wilkinsonn, no singular solution bad badly conditioned
function eg_system_7(d)
    @var x
    System([expand(prod(x - i for i = 1:d))])
end

function eg_system_8(d)
    @var x
    System([(x - 10)^d])
end

# "Winding Number Family d=$d" for d = 2:2:6
function eg_system_9()
    @var x y
    a = [0.257, -0.139, -1.73, -0.199, 1.79, -1.32]
    f1 = (a[1] * x^d + a[2] * y) * (a[3] * x + a[4] * y) + 1
    f2 = (a[1] * x^d + a[2] * y) * (a[5] * x + a[6] * y) + 1
    System([f1, f2])
end


# 693 non-singular solutions
function eg_system_10()
    # Communicated by Mohab Safey El Din
    @var x y z
    System([
        -9091098778555951517 * x^3 * y^4 * z^2 +
        5958442613080401626 * y^2 * z^7 +
        17596733865548170996 * x^2 * z^6 - 17979170986378486474 * x * y * z^6 -
        2382961149475678300 * x^4 * y^3 - 15412758154771986214 * x * y^3 * z^3 + 133,
        -10798198881812549632 * x^6 * y^3 * z - 11318272225454111450 * x * y^9 -
        14291416869306766841 * y^9 * z - 5851790090514210599 * y^2 * z^8 +
        15067068695242799727 * x^2 * y^3 * z^4 +
        7716112995720175148 * x^3 * y * z^3 +
        171,
        13005416239846485183 * x^7 * y^3 + 4144861898662531651 * x^5 * z^4 -
        8026818640767362673 * x^6 - 6882178109031199747 * x^2 * y^4 +
        7240929562177127812 * x^2 * y^3 * z +
        5384944853425480296 * x * y * z^4 +
        88,
    ])
end