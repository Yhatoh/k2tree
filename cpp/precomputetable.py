
table = []
depth_table = []

for j in range(2 ** 16):
    s_b = list(bin(j))[2:]
    s_b.reverse()
    while(len(s_b) < 16):
        s_b.append("0")
    pos = 8
    depth = 0
    min_exc = 18
    for k in range(16):
        if s_b[k] == "1":
            depth += 1
        else:
            depth += -1
        min_exc = min(depth, min_exc)
    table.append(min_exc)
    depth_table.append(depth)


print("int8_t exc_min_micro[65536] = {")

print("{", end="")
for j in range(2**16 - 1):
    print(table[j], end=",")
print(table[2**16 - 1], end="};")
print()

print("int8_t exc_micro[65536] = {")

print("{", end="")
for j in range(2**16 - 1):
    print(depth_table[j], end=",")
print(depth_table[2**16 - 1], end="};")
print()
