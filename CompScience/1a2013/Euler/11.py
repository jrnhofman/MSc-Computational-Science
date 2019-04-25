b = open('11.dat').readlines()
b = [x.split() for x in b]

product_sum = 0
k = 0
# horizontal
for i in range(len(b)):
    for j in range(len(b)-3):
        product_sum = max(int(b[i][j]) * int(b[i][j+1]) * int(b[i][j+2]) * int(b[i][j+3]),product_sum)

# vertical
for i in range(len(b)-3):
    for j in range(len(b)):
        product_sum = max(int(b[i][j]) * int(b[i+1][j]) * int(b[i+2][j]) * int(b[i+3][j]),product_sum)

# diagonal down
for i in range(len(b)-3):
    for j in range(len(b)-3):
        product_sum = max(int(b[i][j]) * int(b[i+1][j+1]) * int(b[i+2][j+2]) * int(b[i+3][j+3]),product_sum)

# diagonal up
for i in range(3,len(b)):
    for j in range(len(b)-3):
        product_sum = max(int(b[i][j]) * int(b[i-1][j+1]) * int(b[i-2][j+2]) * int(b[i-3][j+3]),product_sum)

print 'The maximum product of 4 numbers in the grid is',product_sum
