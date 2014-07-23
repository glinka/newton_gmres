import numpy as np
import matplotlib.pyplot as plt

def threshold(A):
    n = A.shape[0]
    m = A.shape[1]
    tol = 1e-14
    for i in range(n):
        for j in range(m):
            if np.abs(A[i,j]) < tol:
                A[i,j] = 0
    return A

def qr_tridiag(A):
    R = np.copy(A)
    # assume R is square upper-Hessenberg
    n = R.shape[0]
    if n == 1:
        return np.ones((1,1)), R
    c = lambda diag_elem, subdiag_elem: diag_elem/np.sqrt(diag_elem*diag_elem + subdiag_elem*subdiag_elem)
    s = lambda diag_elem, subdiag_elem: -subdiag_elem/np.sqrt(diag_elem*diag_elem + subdiag_elem*subdiag_elem)
    Q = np.identity(n)
    Q[0,0] = c(R[0,0], R[1,0])
    Q[1,1] = Q[0,0]
    Q[1,0] = -s(R[0,0], R[1,0])
    Q[0,1] = s(R[0,0], R[1,0])
    coeffs = np.empty(2)
    # calc G_1 H
    R = np.dot(Q.T, R)
    # multiply G^T_1 G^T_2 ... G^T_k
    # from left to right
    # [k,k] will mark the upper left
    # entry of the Given's rotation
    # in the matrix G^T_k being multiplied into
    # G^T_1 G^T_2 ... G^T_{k-1}
    for k in range(1, n-2):
        coeffs[0] = c(R[k,k], R[k+1,k])
        coeffs[1] = s(R[k,k], R[k+1,k])
        counter = 0
        copy = np.copy(Q[:,k])
        for j in range(k, k+2):
            for i in range(k+1):
                Q[i,j] = coeffs[counter]*copy[i]
            counter = 1
        Q[k+1, k] = -coeffs[1]
        Q[k+1, k+1] = coeffs[0]
        # calc Q_n R
        # since tridiagonal, iterate from k -> k+1
        for j in range(k, k+3):
            temp1 = R[k, j]
            temp2 = R[k+1, j]
            R[k, j] = coeffs[0]*temp1 - coeffs[1]*temp2
            R[k+1, j] = coeffs[1]*temp1 + coeffs[0]*temp2
    # to remove if/else, treat k = n-1
    # separately, otherwise updating of G_n R
    # would need some check "if k == n-1"
    # notice that the loop is now range(k, k+2)
    # instead of range(k, k+3)
    k = n-2
    coeffs[0] = c(R[k,k], R[k+1,k])
    coeffs[1] = s(R[k,k], R[k+1,k])
    counter = 0
    copy = np.copy(Q[:,k])
    for j in range(k, k+2):
        for i in range(k+1):
            Q[i,j] = coeffs[counter]*copy[i]
        counter = 1
    Q[k+1, k] = -coeffs[1]
    Q[k+1, k+1] = coeffs[0]
    # calc Q_n R
    # since tridiagonal, iterate from k -> k+1
    for j in range(k, k+2):
        temp1 = R[k, j]
        temp2 = R[k+1, j]
        R[k, j] = coeffs[0]*temp1 - coeffs[1]*temp2
        R[k+1, j] = coeffs[1]*temp1 + coeffs[0]*temp2
    return [Q, R]

# ain't nerr use dis

def qr_upper_hessenberg(A):
    R = np.copy(A)
    # assume R is square upper-Hessenberg
    n = R.shape[0]
    if n == 1:
        return np.ones((1,1)), R
    c = lambda diag_elem, subdiag_elem: diag_elem/np.sqrt(diag_elem*diag_elem + subdiag_elem*subdiag_elem)
    s = lambda diag_elem, subdiag_elem: -subdiag_elem/np.sqrt(diag_elem*diag_elem + subdiag_elem*subdiag_elem)
    Q = np.identity(n)
    Q[0,0] = c(R[0,0], R[1,0])
    Q[1,1] = Q[0,0]
    Q[1,0] = -s(R[0,0], R[1,0])
    Q[0,1] = s(R[0,0], R[1,0])
    coeffs = np.empty(2)
    # calc G_1 H
    R = np.dot(Q.T, R)
    # multiply G^T_1 G^T_2 ... G^T_k
    # from left to right
    # [k,k] will mark the upper left
    # entry of the Given's rotation
    # in the matrix G^T_k being multiplied into
    # G^T_1 G^T_2 ... G^T_{k-1}
    for k in range(1, n-1):
        coeffs[0] = c(R[k,k], R[k+1,k])
        coeffs[1] = s(R[k,k], R[k+1,k])
        counter = 0
        copy = np.copy(Q[:,k])
        for j in range(k, k+2):
            for i in range(k+1):
                Q[i,j] = coeffs[counter]*copy[i]
            counter = 1
        Q[k+1, k] = -coeffs[1]
        Q[k+1, k+1] = coeffs[0]
        # calc Q_n R
        for j in range(k, n):
            temp1 = R[k, j]
            temp2 = R[k+1, j]
            R[k, j] = coeffs[0]*temp1 - coeffs[1]*temp2
            R[k+1, j] = coeffs[1]*temp1 + coeffs[0]*temp2
    return [Q, R]

def arnoldi_iter(A, V_init, kmax, H_init=None, f_init=None, zero_tol=1e-14):
    # V_init must be either a single vector with norm 1
    # or a set of k orthonormal basis vectors for the k^{th}
    # Krylov subspace
    n = A.shape[0]
    k_init = V_init.shape[1]

    V = np.empty((n, kmax))
    H = np.zeros((kmax, kmax))
    f = np.zeros(n)

    V[:,:k_init] = np.copy(V_init)
    if k_init == 1:
        f = np.dot(A, V[:, 0])
        H[0,0] = np.dot(V[:,0], f)
        f = f - H[0,0]*V[:,0]
    elif H_init is not None and f_init is not None:
        f = np.copy(f_init)
        H[:k_init, :k_init] = np.copy(H_init)
    else:
        print "incorrect input"
        quit()

    for k in range(k_init, kmax):
        fnorm = np.linalg.norm(f)
        if fnorm < zero_tol:
            print 'initial vector only spans', k, 'dimensional subspace'
            print '|| f || = 0'
            return [V, H, f]
        V[:, k] = f/fnorm
        H[k, k-1] = fnorm
        # k = k + 1
        f = np.dot(A, V[:, k])
        for i in range(k+1):
            H[i, k] = np.dot(V[:, i], f)
            f = f - H[i, k]*V[:, i]
    return [V, H, f]

def qr_implicitly_shifted_tridiag(M, tol=1e-11, zero_tol=1e-13, qr_maxiter=10):

    # Francis shifts ftw

    A = np.copy(M)
    n = A.shape[0]
    thetas = np.empty(n)
    S = np.zeros((n, n))

    if n == 1:
        return [1, M]
    # should we compute the eigenpairs of 2x2 matrices directly?
    # for now, do
    elif n == 2:
        a = A[0,0]
        b = A[0,1]
        c = A[1,1]
        d = a*c - b*b
        thetas[0] = (a+c + np.sqrt(np.power(a+c, 2) - 4*d))/2.0
        thetas[1] = (a+c - np.sqrt(np.power(a+c, 2) - 4*d))/2.0
        S[1,0] = 1
        S[0,0] = (thetas[0] - b - c)*S[1,0]/(a + b - thetas[0])
        S[:,0] = S[:,0]/np.linalg.norm(S[:,0])
        S[1,1] = -1
        S[0,1] = (thetas[1] - b - c)*S[1,1]/(a + b - thetas[1])
        S[:,1] = S[:,1]/np.linalg.norm(S[:,1])
        return [S, thetas]

    I = np.identity(n)
    gersh_rings = np.empty(n)
    Q = np.identity(n)
    err = 1
    iters = 0
    # Lower right matrix of the form
    # | a b |
    # | b c |
    wilk_shift = lambda a, b, c: c - np.sign((a - c)/2.0)*b*b/(np.abs((a - c)/2.0) + np.sqrt(np.power((a - c)/2.0, 2) + b*b))
    while err > tol and iters < qr_maxiter:
        # use lower right value of A^{k} as shift
        shift = wilk_shift(A[-2, -2], A[-2, -1], A[-1, -1])
        if np.iscomplex(shift):
            print '***********************'
            print 'francis shifting, bitch'
            print '***********************'
            # use the double, implicit Francis shifts
            shift1 = shift
            shift2 = np.conj(shift)
            realshift = np.real(shift)
            complexshift = np.imag(shift)
            shift_squaredmag = realshift*realshift + complexshift*complexshift
            s = 2*realshift
            t = shift_squaredmag
            l = np.zeros(n)
            l[0] = A[0,0]*A[0,0] + A[0,1]*A[0,1] - s*A[0,0] + t
            l[1] = A[0,1]*(A[0,0] + A[1,1] - s)
            l[2] = A[0,1]*A[2,1]
            reflect_vect = np.copy(l)
            reflect_vect[0] = reflect_vect[0] + np.sign(l[0])*np.linalg.norm(l)
            reflect_vect.shape = (n, 1)
            reflect_mat = np.identity(n) - 2*np.dot(reflect_vect, reflect_vect.T)/np.dot(reflect_vect.T, reflect_vect)
            # francis shifts ftloss
            A = np.dot(reflect_mat.T, np.dot(A, reflect_mat))
            Q = np.dot(reflect_mat, Q)
            reflect_vect = np.empty(3)
            I_small = np.identity(3)
            # chase dat bulj
            for i in range(1, n-3):
                # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! #
                # whole loop could be likely be sped up by   #
                # manually multiplying matrices, given their #
                # tridiagonal structure                      #
                # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! #
                reflect_vect[:] = A[i+1:i+4, i]
                reflect_vect[0] = reflect_vect[0] + np.sign(reflect_vect[0])*np.linalg.norm(reflect_vect)
                reflect_mat = np.copy(I)
                reflect_mat[i+1:i+4, i+1:i+4] = I_small - 2*np.dot(reflect_vect, reflect_vect.T)/np.dot(reflect_vect.T, reflect_vect)
                A = np.dot(reflect_mat.T, np.dot(A, reflect_mat))
                Q = np.dot(reflect_mat, Q)
            # last reflection vector will be 2x1
            reflect_vect = A[-2:, -3]
            reflect_vect[0] = reflect_vect[0] + np.sign(reflect_vect[0])*np.linalg.norm(reflect_vect)
            reflect_mat = np.copy(I)
            reflect_mat[-2:, -2:] = np.identity(2) - 2*np.dot(reflect_vect, reflect_vect.T)/np.dot(reflect_vect.T, reflect_vect)
            A = np.dot(reflect_mat.T, np.dot(A, reflect_mat))
            Q = np.dot(reflect_mat, Q)

        else:
            # use single Wilkinson
            Qi, Ri = qr_tridiag(A - shift*I)
            Q = np.dot(Q, Qi)
            A = np.dot(Ri, Qi) + shift*I
        for i in range(n-1):
            if np.abs(A[i, i+1]) < zero_tol:
                # print 'splitting at n =', n, 'after', iters, 'iterations'
                A[i, i+1] = 0
                A[i+1, i] = 0
                S_top, thetas_top = qr_implicitly_shifted_tridiag(A[:i+1, :i+1], tol=tol, zero_tol=zero_tol, qr_maxiter=qr_maxiter)
                S_bottom, thetas_bottom = qr_implicitly_shifted_tridiag(A[i+1:, i+1:], tol=tol, zero_tol=zero_tol, qr_maxiter=qr_maxiter)
                thetas[:i+1] = thetas_top
                thetas[i+1:] = thetas_bottom
                S[:i+1, :i+1] = S_top
                S[i+1:, i+1:] = S_bottom
                S = np.dot(Q, S)
                # print threshold(np.dot(S.T, S)), S.shape
                return [S, thetas]
            # error is measured by Gershgorin's rings
            gersh_rings[0] = np.abs(A[0, 1])
            gersh_rings[n-1] = np.abs(A[n-1, n-2])
            for j in range(1, n-1):
                gersh_rings[j] = np.abs(A[j, j-1]) + np.abs(A[j, j+1])
        err = np.max(gersh_rings)
        iters = iters + 1
    if iters == qr_maxiter:
        print 'qr failed to converge with n =', n
    S = Q
    thetas = np.diagonal(A)
    return [S, thetas]
    


def qr_explicitly_shifted_tridiag(M, tol=1e-11, zero_tol=1e-13, qr_maxiter=10):
    A = np.copy(M)
    n = A.shape[0]
    I = np.identity(n)
    thetas = np.empty(n)
    gersh_rings = np.empty(n)
    S = np.zeros((n, n))
    # compute directly the eigenpairs of 2x2 matrices
    if n == 2:
        a = A[0,0]
        b = A[0,1]
        c = A[1,1]
        d = a*c - b*b
        thetas[0] = (a+c + np.sqrt(np.power(a+c, 2) - 4*d))/2.0
        thetas[1] = (a+c - np.sqrt(np.power(a+c, 2) - 4*d))/2.0
        S[1,0] = 1
        S[0,0] = (thetas[0] - b - c)*S[1,0]/(a + b - thetas[0])
        S[:,0] = S[:,0]/np.linalg.norm(S[:,0])
        S[1,1] = -1
        S[0,1] = (thetas[1] - b - c)*S[1,1]/(a + b - thetas[1])
        S[:,1] = S[:,1]/np.linalg.norm(S[:,1])
        return [S, thetas]
    Q = np.identity(n)
    err = 1
    iters = 0
    while err > tol and iters < qr_maxiter:
        # use lower right value of A^{k} as shift
        shift = A[-1, -1]
        Qi, Ri = qr_tridiag(A - shift*I)
        Q = np.dot(Q, Qi)
        A = np.dot(Ri, Qi) + shift*I
        for i in range(n-1):
            if np.abs(A[i, i+1]) < zero_tol:
                # print 'splitting at n =', n, 'after', iters, 'iterations'
                A[i, i+1] = 0
                A[i+1, i] = 0
                S_top, thetas_top = qr_explicitly_shifted_tridiag(A[:i+1, :i+1], tol=tol, zero_tol=zero_tol, qr_maxiter=qr_maxiter)
                S_bottom, thetas_bottom = qr_explicitly_shifted_tridiag(A[i+1:, i+1:], tol=tol, zero_tol=zero_tol, qr_maxiter=qr_maxiter)
                thetas[:i+1] = thetas_top
                thetas[i+1:] = thetas_bottom
                S[:i+1, :i+1] = S_top
                S[i+1:, i+1:] = S_bottom
                S = np.dot(Q, S)
                # print threshold(np.dot(S.T, S)), S.shape
                return [S, thetas]
            # error is measured by Gershgorin's rings
            #!!!!!!!!!!!!!!!!!!!! HOW SHOULD YOU ACTUALLY MEASURE ERROR? !!!!!!!!!!!!!!!!!!!!
            gersh_rings[0] = np.abs(A[0, 1])
            gersh_rings[n-1] = np.abs(A[n-1, n-2])
            for j in range(1, n-1):
                gersh_rings[j] = np.abs(A[j, j-1]) + np.abs(A[j, j+1])
        err = np.max(gersh_rings)
        iters = iters + 1
    if iters == qr_maxiter:
        print 'qr failed to converge with n =', n
    S = Q
    thetas = np.diagonal(A)
    return [S, thetas]


def implicitly_restarted_arnoldi_symmetric(A, v, k, p, tol=1e-14, iram_maxiter=10, **kwargs):
    m = k+p
    V, H, f = arnoldi_iter(A, v/np.linalg.norm(v), m)

    beta = H[k, k-1]
    e_m = np.zeros(m)
    e_m[-1] = 1
    S = np.empty((m,m))
    thetas = np.empty(m)
    mu = np.empty(p)
    sorted_indices = None

    iters = 0
    err = 1
    while err > tol and iters < iram_maxiter:
        V, H, f = arnoldi_iter(A, V, k+p, H, f)
        S, thetas = qr_implicitly_shifted_tridiag(H, **kwargs)
        sorted_indices = np.argsort(np.abs(thetas))
        errs = np.empty(k)
        fnorm = np.linalg.norm(f)
        for i in range(1, k+1):
            errs[i-1] = fnorm*np.abs(np.dot(e_m, S[:, sorted_indices[-i]]))
        err = np.max(errs)
        mu[:] = thetas[sorted_indices[:p]]
        Q = np.identity(m)
        I = np.identity(m)
        if err > tol:
            for i in range(p):
                Qi, Ri = qr_tridiag(H - mu[i]*I)
                H = np.dot(Qi.T, np.dot(H, Qi))
                Q = np.dot(Q, Qi)
            beta = H[k, k-1]
            sigma = Q[m-1, k-1]
            f = beta*V[:, k] + sigma*f
            V = np.dot(V, Q[:,:k])
            H = H[:k, :k]
            iters = iters + 1
    if iters == iram_maxiter:
        print 'iram failed to converge'
        quit()
    print 'iram exit err:', err
    return [S[:,sorted_indices[-k:]], thetas[sorted_indices[-k:]], V]

if __name__=="__main__":
    # tests 
    n = 500
    A = np.empty((n,n))
    for i in range(n):
        A[i,i:] = np.random.randn(n-i)
        A[i:, i] = A[i, i:]
    v = np.random.rand(n)
    v = v/np.linalg.norm(v)
    v.shape = (n, 1)
    # V, H, f = arnoldi_iter(A, v, k)
    # f.shape = (n, 1)
    # e_k = np.zeros((k, 1))
    # e_k[k-1] = 1
    # print threshold(np.dot(A, V) - (np.dot(V, H) + np.dot(f, e_k.T)))
    # print threshold(np.dot(V.T, V))
    # # print threshold(np.dot(A, V) - (np.dot(V, H) + np.dot(f, e_k.T)))
    # V, H, f = arnoldi_iter(A, V[:,:5], k, H[:5,:5], H[5,4]*V[:,5])
    # f.shape = (n, 1)
    # e_k = np.zeros((k, 1))
    # e_k[k-1] = 1
    # # print H
    # print threshold(np.dot(A, V) - (np.dot(V, H) + np.dot(f, e_k.T)))
    # print threshold(np.dot(V.T, V))
    k = 4
    p = 8
    # kwargs = {'maxiter': 25}
    S, thetas, V = implicitly_restarted_arnoldi_symmetric(A, v, k, p, iram_maxiter=100, qr_maxiter=2000)
    eigs = np.linalg.eigvalsh(A)
    si = np.argsort(np.abs(eigs))
    print 'squared eigenvalue error:', np.linalg.norm(eigs[si[-k:]] - thetas)
    print thetas
    print A.shape, V.shape, S.shape
    print 'squared first eigenvector error:', np.linalg.norm(np.dot(A, np.dot(V, S[:,-1])) - np.dot(V, S[:,-1])*thetas[-1])
