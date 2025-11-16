      SUBROUTINE PENTA(N,E,A,D,C,F,B,X)

!   RESULTS:  matrix has 5 bands, EADCF, with D being the main diagonal,
!   E and A are the lower diagonals, and C and F are the upper diagonals.
!
!     E is defined for rows i = 3:N, but is defined as E(1) to E(N-2)
!     A is defined for rows i = 2:N, but is defined as A(1) to A(N-1)
!     D is defined for rows i = 1:N
!     C is defined for rows i = 1:N-1, but the last element isn't used
!     F is defined for rows i = 1:N-2, but the last 2 elements aren't used
!
!   B is the right-hand side
!   X is the solution vector

      implicit none
      integer i, n
      double precision E(N), A(N), D(N), C(N), F(N), B(N), X(N), XMULT

      DO 2 i = 2, N-1
        XMULT = A(i-1+1) / D(i-1)
        D(i)  = D(i)  - XMULT * C(i-1)
        C(i)  = C(i)  - XMULT * F(i-1)
        B(i)  = B(i)  - XMULT * B(i-1)

        XMULT   = E(i-1+2) / D(i-1)
        A(i+1)  = A(i+1)  - XMULT * C(i-1)
        D(i+1)  = D(i+1)  - XMULT * F(i-1)
        B(i+1)  = B(i+1)  - XMULT * B(i-1)
    2 CONTINUE

      XMULT = A(N-1+1) / D(N-1)
      D(N)  = D(N) - XMULT * C(N-1)
      X(N)  = (B(N) - XMULT * B(N-1)) / D(N)
      X(N-1) = (B(N-1) - C(N-1) * X(N)) / D(N-1)

      DO 3 i = N-2, 1, -1
        X(i) = (B(i) - F(i) * X(i+2) - C(i) * X(i+1)) / D(i)
    3 CONTINUE

      RETURN
      END
