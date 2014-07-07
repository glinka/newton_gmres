import numpy as np

def gmres(A, x0, b, tol, kmax):
    n = A.shape[0]
    r = b - np.dot(A, x0)
    V = np.empty((n, kmax+1))
    H = np.zeros((kmax+1, kmax))
    rho = np.linalg.norm(r, 2)
    beta = rho
    b_norm = np.linalg.norm(b, 2)
    V[:,0] = r/rho
    g = np.zeros(kmax + 1)
    c = np.zeros(kmax + 1)
    s = np.zeros(kmax + 1)
    g[0] = rho
    k = 0
    while rho > tol*b_norm and k < kmax:
        k = k + 1
        V[:, k] = np.dot(A, V[:, k-1])
        for i in range(k):
            H[i, k-1] = np.dot(V[:, k], V[:, i])
            V[:, k] = V[:, k] - H[i, k-1]*V[:, i]
        v_norm = np.linalg.norm(V[:, k])
        H[k, k-1] = v_norm
        # should test for orthogonality
        V[:, k] = V[:, k]/v_norm
        for i in range(k-1):
            temp1 = c[i]*H[i, k-1] - s[i]*H[i+1, k-1]
            temp2 = s[i]*H[i, k-1] + c[i]*H[i+1, k-1]
            H[i, k-1] = temp1
            H[i+1, k-1] = temp2
        nu = np.sqrt(H[k-1, k-1]*H[k-1, k-1] + H[k, k-1]*H[k, k-1])
        c[k-1] = H[k-1, k-1]/nu
        s[k-1] = -H[k, k-1]/nu
        H[k-1, k-1] = c[k-1]*H[k-1, k-1] - s[k-1]*H[k, k-1]
        H[k, k-1] = 0
        temp1 = g[k-1]*c[k-1] - s[k-1]*g[k]
        temp2 = g[k-1]*s[k-1] + c[k-1]*g[k]
        g[k-1] = temp1
        g[k] = temp2
        rho = np.abs(g[k])
        # print rho
    y = np.zeros(k)
    for i in range(k-1, -1, -1):
        sum_ = 0.0
        for j in range(i+1, k):
            sum_ = sum_ + y[j]*H[i,j]
        y[i] = (g[i] - sum_)/H[i, i]
    return x0 + np.dot(V[:,:k], y), k, rho, b_norm, tol

def fd(F, x, v, h):
    v_norm = np.linalg.norm(v, 2)
    x_norm = np.linalg.norm(x, 2)
    return v_norm*(F(x + h*x_norm*v/v_norm) - F(x))/(h*x_norm)

def fd_xzero(F, x, v, h):
    v_norm = np.linalg.norm(v, 2)
    return v_norm*(F(h*v/v_norm) - F(x))/h

def fdgmres(F, x, tol, h, kmax, x0):
    float_tol = 1e-9
    fd_ = fd
    if np.linalg.norm(x, 2) < float_tol:
        fd_ = fd_xzero
    r = 0
    if np.linalg.norm(x0, 2) < float_tol:
        r = -F(x)
    else:
        r = -fd_(F, x, x0) - F(x)
    n = x.shape[0]
    V = np.empty((n, kmax+1))
    H = np.zeros((kmax+1, kmax))
    rho = np.linalg.norm(r, 2)
    beta = rho
    b_norm = np.linalg.norm(F(x), 2)
    V[:,0] = r/rho
    g = np.zeros(kmax + 1)
    c = np.zeros(kmax + 1)
    s = np.zeros(kmax + 1)
    g[0] = rho
    k = 0
    while rho > tol*b_norm and k < kmax:
        k = k + 1
        V[:, k] = fd_(F, x, V[:, k-1], h)
        for i in range(k):
            H[i, k-1] = np.dot(V[:, k], V[:, i])
            V[:, k] = V[:, k] - H[i, k-1]*V[:, i]
        v_norm = np.linalg.norm(V[:, k])
        H[k, k-1] = v_norm
        # should test for orthogonality
        V[:, k] = V[:, k]/v_norm
        for i in range(k-1):
            temp1 = c[i]*H[i, k-1] - s[i]*H[i+1, k-1]
            temp2 = s[i]*H[i, k-1] + c[i]*H[i+1, k-1]
            H[i, k-1] = temp1
            H[i+1, k-1] = temp2
        nu = np.sqrt(H[k-1, k-1]*H[k-1, k-1] + H[k, k-1]*H[k, k-1])
        c[k-1] = H[k-1, k-1]/nu
        s[k-1] = -H[k, k-1]/nu
        H[k-1, k-1] = c[k-1]*H[k-1, k-1] - s[k-1]*H[k, k-1]
        H[k, k-1] = 0
        temp1 = g[k-1]*c[k-1] - s[k-1]*g[k]
        temp2 = g[k-1]*s[k-1] + c[k-1]*g[k]
        g[k-1] = temp1
        g[k] = temp2
        rho = np.abs(g[k])
        # print rho
    y = np.zeros(k)
    for i in range(k-1, -1, -1):
        sum_ = 0.0
        for j in range(i+1, k):
            sum_ = sum_ + y[j]*H[i,j]
        y[i] = (g[i] - sum_)/H[i, i]
    return x0 + np.dot(V[:,:k], y)

def newton_fdgmres(F, x0, rel_tol, abs_tol, gmres_tol, gmres_h, itermax):
    n = x0.shape[0]
    sqrtn = np.sqrt(n)
    iter_ = 0
    x = x0
    rho = np.linalg.norm(F(x), 2)/sqrtn
    err = rho
    gmres_x0 = np.zeros(n)
    while rho > rel_tol*err + abs_tol and iter_ < itermax:
        # assume constant forcing of gmres_tol
        iter_ = iter_ + 1
        x = x + fdgmres(F, x, gmres_tol, gmres_h, n, gmres_x0)
        rho = np.linalg.norm(F(x), 2)/sqrtn
    return x

def test_f(x):
    return np.array((x[0]*x[0] - 1, x[1]*x[1] - 4, x[2]*x[2] - 9, x[3]*x[3] - 16, x[4]*x[4] - 25))

if __name__=="__main__":
    tol = 1e-6
    h = tol
    itermax = 40
    print newton_fdgmres(test_f, np.random.rand(5), 1e-6, tol, tol, h, itermax)
    # n = 500
    # A = np.random.rand(n, n)
    # b = np.random.rand(n)
    # print np.dot(np.linalg.inv(A), b)
    # print gmres(A, np.dot(np.linalg.inv(A), b)+.05*np.random.rand(n), b, 1e-6, n)
