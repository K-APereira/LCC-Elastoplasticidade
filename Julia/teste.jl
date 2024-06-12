using SparseArrays, LinearSolve, LinearSolvePardiso
dim = 3
@time A = sprand(dim,dim,0.7)
b = rand(dim)

prob = LinearProblem(A, b)
# u = solve(prob,UMFPACKFactorization()).u

# for alg in (
#     MKLPardisoFactorize(),
#     MKLPardisoFactorize(),
#     MKLPardisoIterate(),
#     MKLPardisoIterate())

#     @time u = solve(prob, alg).u
# end
A
b
u = solve(prob, MKLPardisoFactorize())
u
A*u
# a = spzeros(3,3)
# a[1,2] += zero(Float64)
# a