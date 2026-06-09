from sympy import Matrix, Poly, groebner, latex, symbols, trace

a, b, c = symbols('a b c')

# image 1
q1 = [(31, 3), (39, 6),(30, 18),(51, 44),(26, 41)]
# image 2
q2 = [(32, 14),(43, 9),(41, 27),(65, 27),(46, 44)]

# set up 5x9 matrix
Q = Matrix([
    [x1*x2, y1*x2, x2, x1*y2, y1*y2, y2, x1, y1, 1]
    for (x1, y1), (x2, y2) in zip(q1, q2)
])

# row reduce to get nullspace
E_vec = Q.nullspace()

# E=aE_1+bE_2+cE_3+dE_4
# note we rescale so d=1
E = a*E_vec[0] + b*E_vec[1] + c*E_vec[2] + E_vec[3]
E = Matrix(E).reshape(3, 3)

# apply the Demazure constraint to get the Ideal
I = 2*E*E.T*E - trace(E*E.T)*E

# 9 max deg-3 polynom
I = list(I)

# Calculate Groebner basis under grevlex order
G = groebner(I, a, b, c, order='grevlex')

# lex order for easy elimination
G = groebner(G, a, b, c, order='lex')

print(Poly(G[-1]).nroots())

