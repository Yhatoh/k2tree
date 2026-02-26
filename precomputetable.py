
table = []

for i in range(31):
    table.append([])

for i in range(31):
    for j in range(2 ** 8):
        s_b = list(bin(j))[2:]
        s_b.reverse()
        while(len(s_b) < 8):
            s_b.append("0")
        depth = i + 1
        pos = 8
        for k in range(8):
            if s_b[k] == "1":
                depth += 1
            else:
                depth += -1
            if depth == 0:
                pos = k
                break
        table[i].append(pos)


print("uint8_t end_tree[31][2048] = {")

for i in range(31):
    print("{", end="")
    for j in range(2**8 - 1):
        print(table[i][j], end=",")
    print(table[i][2**8 - 1], end="}")
    print()
print("}")
