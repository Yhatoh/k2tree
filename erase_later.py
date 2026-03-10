def count_110100(num):
    x = int(num, 2)
    b1 = x
    b2 = x << 1
    b3 = ~(x << 2)
    b4 = x << 3
    b5 = ~(x << 4)
    b6 = ~(x << 5)

    n = (b1 & b2 & b3 & b4 & b5 & b6) & ((1<<64)-1)
    n = n.bit_count()
    return (n,num.count("110100"))

print(count_110100("1111101001111000"))
