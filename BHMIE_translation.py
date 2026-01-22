import math

def bhmie(x: float, refrel: complex, nang: int):
    y = x * refrel
    xstop = x + 4.0 * (x ** (1.0 / 3.0)) + 2.0
    nstop = int(xstop)
    ymod = abs(y)
    nmx = int(max(xstop, ymod)) + 15
    dang = (math.pi / 2.0) / float(nang - 1)
    theta = [j * dang for j in range(nang)]
    amu = [math.cos(t) for t in theta]
    d = [0j] * (nmx + 1)   
    d[nmx] = 0j
    for n in range(1, nmx):  
        rn = float(nmx - n + 1)
        i = nmx - n
        d[i] = (rn / y) - (1.0 / (d[i + 1] + rn / y))
    pio = [0.0] * nang
    pi1 = [1.0] * nang
    nn = 2 * nang - 1
    s1 = [0j] * nn
    s2 = [0j] * nn
    psi0 = math.cos(x)
    psi1 = math.sin(x)
    chi0 = -math.sin(x)
    chi1 =  math.cos(x)
    xi0 = complex(psi0, -chi0)
    xi1 = complex(psi1, -chi1)
    qsca = 0.0
    for n in range(1, nstop + 1):
        dn = float(n)
        rn = dn
        fn = (2.0 * rn + 1.0) / (rn * (rn + 1.0))
        psi = (2.0 * dn - 1.0) * psi1 / x - psi0
        chi = (2.0 * rn - 1.0) * chi1 / x - chi0
        xi = complex(psi, -chi)
        an_num = (d[n] / refrel + rn / x) * psi - psi1
        an_den = (d[n] / refrel + rn / x) * xi  - xi1
        an = an_num / an_den
        bn_num = (refrel * d[n] + rn / x) * psi - psi1
        bn_den = (refrel * d[n] + rn / x) * xi  - xi1
        bn = bn_num / bn_den
        qsca += (2.0 * rn + 1.0) * (abs(an) ** 2 + abs(bn) ** 2)
        pi = [0.0] * nang
        tau = [0.0] * nang
        p = 1.0 if ((n - 1) % 2 == 0) else -1.0
        t = 1.0 if (n % 2 == 0) else -1.0
        for jF in range(1, nang + 1):  
            j = jF - 1
            jjF = 2 * nang - jF 
            jj = jjF - 1
            pi[j] = pi1[j]
            tau[j] = rn * amu[j] * pi[j] - (rn + 1.0) * pio[j]
            s1[jF - 1] += fn * (an * pi[j] + bn * tau[j])
            s2[jF - 1] += fn * (an * tau[j] + bn * pi[j])
            if jF != jjF:
                s1[jj] += fn * (an * pi[j] * p + bn * tau[j] * t)
                s2[jj] += fn * (an * tau[j] * t + bn * pi[j] * p)
        psi0, psi1 = psi1, psi
        chi0, chi1 = chi1, chi
        xi1 = complex(psi1, -chi1)
        if n < nstop:
            rn_next = float(n + 1)
            for j in range(nang):
                pi1[j] = ((2.0 * rn_next - 1.0) / (rn_next - 1.0)) * amu[j] * pi[j] \
                         - (rn_next / (rn_next - 1.0)) * pio[j]
                pio[j] = pi[j]
    qsca = (2.0 / (x * x)) * qsca
    qext = (4.0 / (x * x)) * (s1[0].real)
    qback = (4.0 / (x * x)) * (abs(s1[nn - 1]) ** 2)
    return s1, s2, qext, qsca, qback

def main():
    print("\nSPHERE SCATTERING PROGRAM\n")
    refmed = 1.0
    refre = 1.55
    refim = 0.0
    refrel = complex(refre, refim) / refmed
