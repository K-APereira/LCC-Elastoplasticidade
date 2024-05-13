

function  met1(m)
    n = broadcast(abs,m)
    return maximum(n[3,:])
end

function met2(n)
    mx = abs(maximum(n[3,:]))
    min = abs(minimum(n[3,:]))
    return max(mx,min)
end

@time for i in 1:5
    a = rand(-1.3:0.01:1,100,100)
    println(met1(a))
end