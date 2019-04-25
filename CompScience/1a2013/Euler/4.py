palindrome = 0

def isPalindrome(product):
    product = str(product)
    numbers = [x for x in product]
    #print numbers
    return (numbers == numbers[::-1])

for i in range(100,1000):
    for j in range(i,1000):
        product = i*j
        if isPalindrome(product):
            if(product > palindrome):
                palindrome = product

print 'The largest palindrome made from the product of 2 3-digit numbers is',palindrome
