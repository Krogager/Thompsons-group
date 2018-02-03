
"""
This goal of this module is to compute the normal form of elements in
Thompsons group F.

REF: this is based on an algorithm from the article "Thompson's group and
public key cryptography" by Shpilrain and Ushakov, availble at
http://www.sci.ccny.cuny.edu/~shpil/thomcryp.pdf

A word in Thompsons group will be represented as a list of 2-lists, each 2-list
consisting of an index (in {0,1,2,...}) and a number in {1,-1} (the exponent).
So, the word x_1^2x_3x_2^{-1} will be represented as
[ [1,1], [1,1], [1,3], [2,-1] ]
"""

def delta(w, k):
    """ w is a word and k is an integer. Returns the corresponding word where
    every generator has its index increased by k. """
    for i in range(len(w)):
        w[i][0] = w[i][0] + k
    return w

def merge_np(n, p, a, b):
    """ Takes as input two seminormal forms n and p, where n only has negative
    exponents and p only has positive exponents. a and b are integers.
    Returns a seminormal form for the element delta_a(n)delta_b(p). """
    if len(n) == 0 or len(p) == 0:
        return delta(n,a) + delta(p,b)
    elif n[-1][0] + a == p[0][0] + b:
        n.pop()
        p.pop(0)
        return merge_np(n, p, a, b)
    elif n[-1][0] + a < p[0][0] + b:
        x = [ [n[-1][0] + a, -1] ]
        n.pop()
        w = merge_np(n, p, a, b+1)
        return w + x
    elif n[-1][0] + a > p[0][0] + b:
        x = [ [p[0][0] + b, 1] ]
        p.pop(0)
        w = merge_np(n, p, a+1, b)
        return x + w

def merge_pp(p, q):
    """ Takes as input two seminormal forms p and q with only positive exponents.
    Returns a seminormal form for the element pq. """
    if len(p) == 0 or len(q) == 0:
        return p + q
    if p[-1][0] <= q[0][0]:
        return p + q
    elif p[-1][0] > q[0][0]:
        x = [ [p[-1][0] + 1, 1] ]
        p.pop()
        y = [ q[0] ]
        q.pop(0)
        return merge_pp(merge_pp(p, y), merge_pp(x, q))


def merge_nn(n, m):
    """ Takes as input two seminormal forms n and m with only negative exponents.
    Returns a seminormal form for the element nm. """
    if len(n) == 0 or len(m) == 0:
        return n + m
    if n[-1][0] >= m[0][0]:
        return n + m
    elif n[-1][0] < m[0][0]:
        x = [ n[-1] ]
        n.pop()
        y = [ [m[0][0] + 1, -1] ]
        m.pop(0)
        return merge_nn(merge_nn(n, y), merge_nn(x, m))

def pos(w):
    """ Takes as imput a seminormal form w and returns the part of w with
    only positive exponents."""
    if len(w) == 0:
        return w
    s = -1
    for i in range(len(w)):
        if w[i][1] == 1:
            s = i
    if s == -1:
        return []
    else:
        return w[:s+1]

def neg(w):
    """ Takes as imput a seminormal form w and returns the part of w with
    only negative exponents."""
    if len(w) == 0:
        return w
    s = -1
    for i in range(len(w)):
        if w[i][1] == 1:
            s = i
    if s == -1:
        return w
    else:
        return w[s+1:]   

def merge(w, v):
    """ Takes as input two seminormal forms w and v. Returns a seminormal form
    of the element wv."""
    a = merge_np(neg(w), pos(v), 0, 0)
    b = merge_pp(pos(w), pos(a))
    c = merge_nn(neg(a), neg(v))
    return b + c

def seminormalForm(w):
    """ Takes as input a word w and returns a seminormal form for w."""
    if len(w) <= 1:
        return w
    else:
        a = w[:1]
        b = w[1:]
        c = seminormalForm(a)
        d = seminormalForm(b)
        return merge(c, d)

def eraseBadPairs(w):
    """ Takes a seminormal form w as input and returns the normal form of w."""
    d1 = 0
    d2 = 0
    d3 = 0
    u1 = []
    u2 = []
    S1 = []
    S2 = []
    w1 = pos(w)
    w2 = neg(w)
        
    while len(w1) > 0 or len(w2) > 0:
        if len(u1) == 0:
            a = -1
        else:
            a = u1[0][0]
        if len(u2) == 0:
            b = -1
        else:
            b = u2[-1][0]
        if len(S1) == 0:
            eps1 = -1
        else:
            eps1 = S1[-1]
        if len(S2) == 0:
            eps2 = -1
        else:
            eps2 = S2[-1]
        
        if len(w1) > 0 and (len(w2) == 0 or w1[-1][0] > w2[0][0]):
            u1 = [w1[-1]] + u1
            w1.pop()
            S1.append(0)
        elif len(w2) > 0 and (len(w1) == 0 or w2[0][0] > w1[-1][0]):
            u2 = u2 + [w2[0]]
            w2.pop(0)
            S2.append(0)
        elif w1[-1][0] == w2[0][0]:
            if (a != -1 and eps1 != -1 and (a - eps1 == w1[-1][0] or a - eps1 == w1[-1][0] + 1)) or (b != -1 and eps2 != -1 and (b - eps2 == w1[-1][0] or b - eps2 == w1[-1][0] + 1)):
                u1 = [w1[-1]] + u1
                u2 = u2 + [w2[0]]
                w1.pop()
                w2.pop(0)
                S1.append(0)
                S2.append(0)
            else:
                w1.pop()
                w2.pop(0)
                if len(S1) > 0:
                    S1[-1] += 1
                if len(S2) > 0:
                    S2[-1] += 1

    while len(u1) > 0:
        d1 = d1 + S1.pop()
        w1 = w1 + [[u1[0][0] - d1, 1]] 
        u1.pop(0)

    while len(u2) > 0:
        d2 = d2 + S2.pop()
        w2 = [[u2[-1][0] - d2, -1]] + w2
        u2.pop()

    return w1 + w2

def normalForm(w):
    """ Takes a word w as input and returns the normal form of w."""
    u = seminormalForm(w)
    return eraseBadPairs(u)




