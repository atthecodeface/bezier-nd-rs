#!/usr/bin/env python3


class Quadratic:
    def __init__(self, pts):
        self.pts = pts
        pass

    def point_at(self, t):
        u = 1.0 - t
        u2 = u * u
        ut = u * t
        t2 = t * t
        return self.pts[0] * u2 + self.pts[1] * (2 * ut) + self.pts[2] * t2

    def gradient_at(self, t):
        u = 1.0 - t
        u_m_t = u - t
        return 2 * (self.pts[2] * t + self.pts[1] * (u_m_t) - self.pts[0] * u)

    def dsq_at_t(self, t, p):
        dp = self.point_at(t) - p
        return dp * dp

    def d_dt_dsq_at_t(self, t, p):
        u = 1.0 - t
        u_m_t = u - t
        u2 = u * u
        ut = u * t
        t2 = t * t
        u3 = u2 * u
        t3 = t2 * t

        # u^2p[0] * p'(t) / 2
        a = (
            -u2 * u * self.pts[0] * self.pts[0]
            + u_m_t * u2 * self.pts[1] * self.pts[0]
            + t * u2 * self.pts[2] * self.pts[0]
        )
        # utp[1] * p'(t) / 2
        b = (
            -u * ut * self.pts[0] * self.pts[1]
            + u_m_t * ut * self.pts[1] * self.pts[1]
            + t * ut * self.pts[2] * self.pts[1]
        )
        # t^2[2] * p'(t) / 2
        c = (
            -u * t2 * self.pts[0] * self.pts[2]
            + u_m_t * t2 * self.pts[1] * self.pts[2]
            + t3 * self.pts[2] * self.pts[2]
        )
        # p * p'(t) / 2
        d = -u * self.pts[0] * p + u_m_t * self.pts[1] * p + t * self.pts[2] * p

        # Tested:
        # 2*d         === self.gradient_at(t)*p)
        # 2*a+4*b+2*c === self.point_at(t)*self.gradient_at(t), a, b, c)

        # d/dt (p(t)-p)^2 = d/dt 2.(p(t)-p).p'(t)
        #                 = 2 (u^2.p0 + 2ut.p1 + t^2.p2 - p).p'(t)
        #                 = 2 ( u^2.p0.p'(t) + 2ut.p1.p'(t) + t^2.p2.p'(t) - p.p'(t) )
        # p'(t)=(-2u.p0 + 2(u-t).p1 + 2t.p2)
        # d/dt (p(t)-p)^2 = 2 ( 2a + 4b + 2c - 2d)
        return (a + 2 * b + c - d) * 4

    def d_dt_dsq_at_t_v2(self, t, p):
        u = 1.0 - t
        u_m_t = u - t
        u2 = u * u
        ut = u * t
        t2 = t * t
        u3 = u2 * u
        t3 = t2 * t

        p0 = self.pts[0]
        p1 = self.pts[1]
        p2 = self.pts[2]

        # u^2*u = (1-t)^3                                      = -t^3 +3t^2 - 3t + 1
        # u_m_t*u^2 = (1-2t)(1+t^2-2t) = 1+t^2-2t-2t-2t^3+4t^2 = -2t^3 + 5t^2 -4t + 1
        # t*u^2 = t*(1+t^2-2t)                                 =  t^3 -2t^2 + t
        # u_m_t*u*t = (1-2t)*(1-t)*t                           = 2t^3 -3t^2 + t
        # t^2*u = t^2*(1-t)                                    = -t^3  +t^2
        # u_m_t*t^2 = (1-2t)*t^2                               =-2t^3  +t^2
        # a = -u2*u*self.pts[0]*self.pts[0] + u_m_t*u2*self.pts[1]*self.pts[0] + t*u2*self.pts[2]*self.pts[0]
        # a = (t^3 -3t^2 + 3t - 1).p0.p0 + (-2t^3 + 5t^2 -4t + 1).p1.p0 + (t^3 -2t^2 + t).p2.p0
        # a = (p0.p0 -2p1.p0 + p2.p0).t^3 + (-3p0.p0+5p1.p0-2p2.p0).t^2 + (3p0.p0-4p1.p0+1p2.p0)t + (-1.p0.p0+1.p0.p1)
        a = (
            (p0 * p0 - 2 * p1 * p0 + p2 * p0) * t3
            + (-3 * p0 * p0 + 5 * p1 * p0 - 2 * p2 * p0) * t2
            + (3 * p0 * p0 - 4 * p1 * p0 + 1 * p2 * p0) * t
            + (-1 * p0 * p0 + 1 * p0 * p1)
        )
        a_orig = (
            -u2 * u * self.pts[0] * self.pts[0]
            + u_m_t * u2 * self.pts[1] * self.pts[0]
            + t * u2 * self.pts[2] * self.pts[0]
        )
        # a == a_orig

        # utp[1] * p'(t) / 2
        # b = -u.ut.self.pts[0].self.pts[1] + u_m_t.ut.self.pts[1].self.pts[1] + t.ut.self.pts[2].self.pts[1]
        # b = (-t^3 +2t^2 - t).p0.p1 + (2t^3 -3t^2 + t).p1.p1 + (-t^3  +t^2).p2.p1
        # b = (-p0.p1 + 2.p1.p1 - 1.p2.p1).t^3 + (2.p0.p1-3.p1.p1+1p2.p1).t^2 + (-1.p0.p1+1.p1.p1).t
        b = (
            (-p0 * p1 + 2 * p1 * p1 - p2 * p1) * t3
            + (2 * p0 * p1 - 3 * p1 * p1 + 1 * p2 * p1) * t2
            + (-1 * p0 * p1 + 1 * p1 * p1) * t
        )
        b_orig = (
            -u * ut * self.pts[0] * self.pts[1]
            + u_m_t * ut * self.pts[1] * self.pts[1]
            + t * ut * self.pts[2] * self.pts[1]
        )
        # b == b_orig

        # t^2[2] * p'(t) / 2
        # c = -u*t2*self.pts[0]*self.pts[2] + u_m_t*t2*self.pts[1]*self.pts[2] + t3*self.pts[2]*self.pts[2]
        # c = -(-t^3  +t^2).p0.p2 + (-2t^3  +t^2).p1.p2 + t3.p2.p2
        # c = (1.p0.p2 -2.p1.p2 + 1.p2.p2)*t^3  + (-1.p0.p2 + 1.p1.p2).t^2
        c = (p0 * p2 - 2 * p1 * p2 + p2 * p2) * t3 + (-1 * p0 * p2 + 1 * p1 * p2) * t2
        c_orig = (
            -u * t2 * self.pts[0] * self.pts[2]
            + u_m_t * t2 * self.pts[1] * self.pts[2]
            + t3 * self.pts[2] * self.pts[2]
        )
        # c == c_orig
        # p * p'(t) / 2
        # d = -u*self.pts[0]*p + u_m_t*self.pts[1]*p + t*self.pts[2]*p
        # d = (t-1)*self.pts[0]*p + (1-2t)*self.pts[1]*p + t*self.pts[2]*p
        # d = (p0.p -2.p1.p + p.p2).t + (-1.p0.p+1.p1.p)
        d = (p0 * p - 2 * p1 * p + p2 * p) * t + (-p0 * p + p1 * p)
        d_orig = -u * self.pts[0] * p + u_m_t * self.pts[1] * p + t * self.pts[2] * p
        # d == d_orig

        return (a + 2 * b + c - d) * 4

    def d_dt_dsq_at_t_v3(self, t, p):
        u = 1.0 - t
        u_m_t = u - t
        u2 = u * u
        ut = u * t
        t2 = t * t
        u3 = u2 * u
        t3 = t2 * t

        p0 = self.pts[0]
        p1 = self.pts[1]
        p2 = self.pts[2]

        a = (
            (p0 * p0 - 2 * p1 * p0 + p2 * p0) * t3
            + (-3 * p0 * p0 + 5 * p1 * p0 - 2 * p2 * p0) * t2
            + (3 * p0 * p0 - 4 * p1 * p0 + 1 * p2 * p0) * t
            + (-1 * p0 * p0 + 1 * p0 * p1)
        )
        b = (
            (-p0 * p1 + 2 * p1 * p1 - p2 * p1) * t3
            + (2 * p0 * p1 - 3 * p1 * p1 + 1 * p2 * p1) * t2
            + (-1 * p0 * p1 + 1 * p1 * p1) * t
        )
        c = (p0 * p2 - 2 * p1 * p2 + p2 * p2) * t3 + (-1 * p0 * p2 + 1 * p1 * p2) * t2
        d = (p0 * p - 2 * p1 * p + p2 * p) * t + (-p0 * p + p1 * p)
        dsq_dt_orig = a + 2 * b + c - d

        dsq_dt_t3 = (
            p0 * p0
            - 2 * p1 * p0
            + p2 * p0
            - 2 * p0 * p1
            + 4 * p1 * p1
            - 2 * p2 * p1
            + p0 * p2
            - 2 * p1 * p2
            + p2 * p2
        )
        dsq_dt_t2 = (
            -3 * p0 * p0
            + 5 * p1 * p0
            - 2 * p2 * p0
            + 4 * p0 * p1
            - 6 * p1 * p1
            + 2 * p2 * p1
            - 1 * p0 * p2
            + 1 * p1 * p2
        )
        dsq_dt_t = (
            3 * p0 * p0
            - 4 * p1 * p0
            + 1 * p2 * p0
            - 2 * p0 * p1
            + 2 * p1 * p1
            - p0 * p
            + 2 * p1 * p
            - p2 * p
        )
        dsq_dt_c = -1 * p0 * p0 + 1 * p0 * p1 + p0 * p - p1 * p

        dsq_dt = dsq_dt_c + t * dsq_dt_t + t2 * dsq_dt_t2 + t3 * dsq_dt_t3
        # dsq_dt == dsq_dt_orig

        return dsq_dt * 4

    def d_dt_dsq_at_t_v4(self, t, p):
        u = 1.0 - t
        u_m_t = u - t
        u2 = u * u
        ut = u * t
        t2 = t * t
        u3 = u2 * u
        t3 = t2 * t

        p0 = self.pts[0]
        p1 = self.pts[1]
        p2 = self.pts[2]

        # dsq_dt_t3 = p0*p0 -2*p1*p0 + p2*p0 -2*p0*p1 + 4*p1*p1 - 2*p2*p1 + p0*p2 -2*p1*p2 + p2*p2
        # dsq_dt_t2 = -3*p0*p0+5*p1*p0-2*p2*p0 + 4*p0*p1-6*p1*p1+2*p2*p1 -1*p0*p2+1*p1*p2
        # dsq_dt_t = 3*p0*p0-4*p1*p0+1*p2*p0 -2*p0*p1+2*p1*p1 - p0*p + 2*p1*p - p2*p
        # dsq_dt_c = -1*p0*p0+1*p0*p1 +p0*p - p1*p
        dsq_dt_t3 = (p0 - 2 * p1 + p2) * (p0 - 2 * p1 + p2)
        dsq_dt_t2 = 3 * (p1 - p0) * (p0 - 2 * p1 + p2)
        dsq_dt_t = (p0 - p1) * (p0 - p1) * 2 + (p0 - 2 * p1 + p2) * (p0 - p)
        dsq_dt_c = (p0 - p1) * (p - p0)

        dsq_dt = dsq_dt_c + t * dsq_dt_t + t2 * dsq_dt_t2 + t3 * dsq_dt_t3

        return dsq_dt * 4

    def d_dt_dsq_at_t_v5(self, t, p):
        t2 = t * t
        t3 = t2 * t

        p0 = self.pts[0]
        p1 = self.pts[1]
        p2 = self.pts[2]
        p012 = p0 - 2 * p1 + p2
        p01 = p1 - p0
        p0p = p0 - p

        dsq_dt_t3 = p012 * p012
        dsq_dt_t2 = 3 * p01 * p012
        dsq_dt_t = 2 * p01 * p01 + p012 * p0p
        dsq_dt_c = p01 * p0p

        dsq_dt = dsq_dt_c + t * dsq_dt_t + t2 * dsq_dt_t2 + t3 * dsq_dt_t3

        return dsq_dt * 4

    pass


B = {}


def bernoulli(i, n):
    """
    n=0:   1
    n=1:   1   1
    n=2:   1   2   1
    n=3:   1   3   3   1
    """
    if (n, i) not in B:
        B[(n, i)] = calc_bernoulli(i, n)
        pass
    return B[(n, i)]


def calc_bernoulli(i, n):
    if i < 0:
        return 0
    if i >= n:
        return 0
    if i == 0 or i == n - 1:
        return 1
    return bernoulli(i - 1, n - 1) + bernoulli(i, n - 1)


print("pub const BINOMIALS_U: &[&[usize]] = &[")
for n in range(1, 40):
    result = f"{1 << (n - 1)}"
    for i in range(n):
        result += f", {bernoulli(i, n)}"
        pass
    print(f"&[{result}],")
    pass
print("];")

b = Quadratic([0.5, 2.5, 3.5])
# b = Quadratic([1.0, 2.5, 0.0])
ts = [t * 0.01 for t in range(101)]
p = 1.2
last_t = 0
last_p_t = 0
last_dsq_t = 0
for t in ts:
    p_t = b.point_at(t)
    dp_dt = b.gradient_at(t)
    dsq_t = b.dsq_at_t(t, p)
    d_dt_dsq_t = b.d_dt_dsq_at_t(t, p)
    d_dt_dsq_t_2 = b.d_dt_dsq_at_t_v2(t, p)
    d_dt_dsq_t_3 = b.d_dt_dsq_at_t_v3(t, p)
    d_dt_dsq_t_4 = b.d_dt_dsq_at_t_v4(t, p)
    d_dt_dsq_t_5 = b.d_dt_dsq_at_t_v5(t, p)
    if last_t > 0:
        e_dp_dt = dp_dt - (p_t - last_p_t) / (t - last_t)
        e_dsq_t = dsq_t - (p_t - p) * (p_t - p)
        e_d_dt_dsq = d_dt_dsq_t - (dsq_t - last_dsq_t) / (t - last_t)
        p_t_m_p_times_dp_dt = (p_t - p) * dp_dt * 2
        e_dt = d_dt_dsq_t - d_dt_dsq_t_5
        print(
            f"{t:14.5} {p_t:14.5} {dp_dt:14.5} {dsq_t:14.5} {e_dsq_t:14.5}~=0 {p_t_m_p_times_dp_dt:14.5}~={d_dt_dsq_t:14.5} {e_dp_dt:14.5}~=0 {e_d_dt_dsq:14.5}~=0 {e_dt:14.5}~=0"
        )
        pass
    last_t = t
    last_p_t = p_t
    last_dsq_t = dsq_t
    pass
