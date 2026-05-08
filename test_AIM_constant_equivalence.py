from sage.all import *

from aim import AIM


def test_aim_constant_equivalence():
    A = AIM(n=128)

    c0 = A.F.random_element()
    c1 = A.F.random_element()

    x = A.F.from_integer(123)
    ct = A.eval(x, A.L1, A.L2, (c0, c1))

    a = A.F.random_element()
    b = A.F.random_element()

    L1new = A.post_scaleL(A.pre_scaleL(A.L1, a**(-A.d1)), a**(A.dstar))
    L2new = A.post_scaleL(A.pre_scaleL(A.L2, a**(-A.d2)), a**(A.dstar))

    xnew = a*x+b
    ctnew = A.eval(xnew, L1new, L2new, (a*c0+b, a*c1+b))
    assert ctnew == ct * a + b
    assert ctnew != ct  # we are not doing stupid stuff
    print("passed")


if __name__ == '__main__':
    test_aim_constant_equivalence()
