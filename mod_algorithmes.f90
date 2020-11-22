!--Module contenant les algorithmes de résolution de systèmes linéaires
Module mod_algorithmes
  Use mod_donnees
  Implicit None

  !Algorithmes
  Contains
    !Algorithme de vérification du résultat
    Subroutine Validation(A,X,b)
      !Déclaration des variables locales
      Real(PR),Dimension(:,:),Intent(In)::A
      Real(PR),Dimension(:),Intent(In)::X, b

      Real(PR),Dimension(:),Allocatable::Diff
      Integer:: i

      !Instructions
      Allocate(Diff(Size(b)))
      Diff = Matmul(A,X)-b
      Do i = 1,Size(b)
        If (abs(Diff(i)) > (10._PR**(-10))) Then !Vérification avec une tolérance de 10**-4
          Print*,"Le résultat du système est incorrect"
          Print*,"AX = ",Matmul(Mat_A,Vect_X)
          Print*,"b = ",Vect_b
          Stop
        End If
      End Do
      Print*,"AX = ",Matmul(Mat_A,Vect_X)
      Print*,"b = ",Vect_b
      Print*,"Le résultat du système est correct"
    End Subroutine Validation

    !Algorithme affichage Matrice
    Subroutine Affichage(A)
      !Déclaration des variables locales
      Real(PR),Dimension(:,:),Intent(In)::A
      Integer::k

      !Instructions
      Do k = 1, Size(A,1)
        Print*,A(k,:)
      End Do
    End Subroutine Affichage

    !Algorithme d'enregistrement des résidus
    Subroutine SaveSol(Liste)
      !Déclaration des variables locales
      Real(PR),Dimension(:),Intent(In)::Liste

      Integer::k
      Character(len=30)::n_string

      !Instructions
      Write(n_string,*)taille_n
      Open(Unit=100, file=trim(adjustl(sys))//'_'//trim(adjustl(meth))//'_'//trim(adjustl(n_string))//'.dat')
      Do k = 1, Size(Liste)
        Write(100,*)k, Liste(k)
      End Do
      Close(100)
    End Subroutine SaveSol

    !Algorithme lecture fichier entrée
    Subroutine Lecture()

    End Subroutine Lecture

    !Algorithme écriture fichier sortie
    Subroutine Ecriture()

    End Subroutine Ecriture

    !----------------------------------------------------------------------!

    Subroutine Init_Test(A,b)
      !Déclaration des variables locales
      Real(PR),Dimension(:,:),Allocatable,Intent(InOut)::A
      Real(PR),Dimension(:),Allocatable,Intent(InOut)::b

      !Instructions
      Allocate(A(3,3)); Allocate(b(3))
      b(1) = 1._PR; b(2) = 2._PR; b(3) = 3._PR

      A(1,1) = 0._PR; A(1,2) = 2._PR; A(1,3) = 2._PR
      A(2,1) = 2._PR; A(2,2) = 0._PR; A(2,3) = -2._PR
      A(3,1) = 2._PR; A(3,2) = -2._PR; A(3,3) = 0._PR

    End Subroutine Init_Test

    Subroutine Init_An(A,b)
      !Déclaration des variables locales
      Real(PR),Dimension(:,:),Allocatable,Intent(InOut)::A
      Real(PR),Dimension(:),Allocatable,Intent(InOut)::b

      Real(PR),Dimension(:,:),Allocatable::In, Bn, BnTBn
      Real(PR),Dimension(:),Allocatable::Somme
      Real(PR)::alpha
      Integer::i,j

      !Instructions
      Allocate(A(taille_n,taille_n)); Allocate(b(taille_n))
      b = 1._PR

      Allocate(In(taille_n,taille_n)); Allocate(Bn(taille_n,taille_n))
      Allocate(BnTBn(taille_n,taille_n))

      In = 0._PR
      Do i = 1, taille_n
        In(i,i) = 1._PR
      End Do

      Call random_seed()
      Do i = 1, taille_n
        Do j = 1, taille_n
          Call random_number(Bn(i,j))
        End Do
      End Do
      BnTBn = Matmul(Transpose(Bn),Bn)

      Allocate(Somme(taille_n))
      Do i = 1, taille_n
        Do j = 1, taille_n
          Somme(i) = Somme(i) + BnTBn(i,j)
        End Do
      End Do
      alpha = maxval(Somme)

      A = In + alpha * BnTBn

    End Subroutine Init_An


    !----------------------------------------------------------------------!

    !Algorithme du Gradient à pas optimal
    Function GPO(A,b) Result(X)
      !Déclaration des variables locales
      Real(PR),Dimension(:,:),Intent(In)::A !Matrice A
      Real(PR),Dimension(:),Intent(In)::b !Vecteur b
      Real(PR),Dimension(:),Allocatable::X !Vecteur X

      Real(PR),Dimension(:),Allocatable::xk, rk, zk, list
      Integer::k
      Real(PR)::alpha,beta

      !Initialisation
      Allocate(xk(Size(A,2)))
      xk = 1._PR !Choix arbitraire de la condition initiale
      Allocate(rk(Size(b)))
      rk = b - Matmul(A,xk) !Erreur initiale
      beta = sqrt(Dot_Product(rk,rk))
      k = 0 !Initialisation du garde-fou
      Print*,"Le résidu à l'itération ",k," vaut : ",beta

      Allocate(list(kmax))
      list(1) = beta

      !Instructions
      Do While ((beta > eps) .and. (k < kmax))
        zk = Matmul(A,rk)
        alpha = Dot_Product(rk,rk) / Dot_Product(zk,rk)
        xk = xk + alpha * rk
        rk = rk - alpha * zk
        beta = sqrt(Dot_Product(rk,rk))
        k = k + 1
        Print*,"Le résidu à l'itération ",k," vaut : ",beta
        list(k) = beta
      End Do

      !Sortie
      Allocate(X(Size(xk)))
      X = xk
      list = Reshape(list,(/k/))
      Call SaveSol(list)

      !Arret du programme si tolérence non atteinte
      If (k > kmax) Then
        Print*,"Tolérence non atteinte : ", sqrt(Dot_Product(rk,rk))
        Print*,"Arret du programme"
        Print*,"!--------------------------------------------------!"
        Stop
      End If

    End Function GPO

    !Algorithme du Gradient Conjugué
    Function GC(A,b) Result(X)
      !Déclaration des variables locales
      Real(PR),Dimension(:,:),Intent(In)::A !Matrice A
      Real(PR),Dimension(:),Intent(In)::b !Vecteur b
      Real(PR),Dimension(:),Allocatable::X !Vecteur X

      Real(PR),Dimension(:),Allocatable::xk, rk, rk1, zk, p, list
      Real(PR)::alpha, beta, gamma
      Integer::k

      !Initialisation
      Allocate(xk(Size(A,2)))
      xk = 1._PR !Choix arbitraire de la condition initiale
      Allocate(rk(Size(b)))
      rk = Matmul(A,xk) - b !Erreur initiale
      p = -rk
      beta = sqrt(Dot_Product(rk,rk))
      k = 0 !Initialisation du garde-fou
      Print*,"Le résidu à l'itération ",k," vaut : ",beta

      Allocate(list(kmax))
      list(1) = beta

      !Instructions
      Do While ((beta > eps) .and. (k < kmax))
        zk = Matmul(A,p)
        alpha = -Dot_Product(rk,p) / Dot_Product(p,zk)
        xk = xk + alpha * p
        rk1 = rk + alpha * zk
        gamma = Dot_Product(rk1,rk1) / Dot_Product(rk,rk)
        p = -rk1 + gamma * p
        rk = rk1
        beta = sqrt(Dot_Product(rk,rk))
        k = k + 1
        Print*,"Le résidu à l'itération ",k," vaut : ",beta
        list(k) = beta
      End Do

      !Sortie
      Allocate(X(Size(xk)))
      X = xk
      list = Reshape(list,(/k/))
      Call SaveSol(list)

      !Arret du programme si tolérence non atteinte
      If (k > kmax) Then
        Print*,"Tolérence non atteinte :", sqrt(Dot_Product(rk,rk))
        Print*,"Arret du programme"
        Print*,"!--------------------------------------------------!"
        Stop
      End If

    End Function GC

    !Algorithme du Résidu Minimum
    Function ResMin(A,b) Result(X)
      !Déclaration des variables locales
      Real(PR),Dimension(:,:),Intent(In)::A !Matrice A
      Real(PR),Dimension(:),Intent(In)::b !Vecteur b
      Real(PR),Dimension(:),Allocatable::X !Vecteur X

      Real(PR),Dimension(:),Allocatable::xk, rk, zk, list
      Real(PR)::alpha,beta
      Integer::k

      !Initialisation
      Allocate(xk(Size(A,2)))
      xk = 1._PR !Choix arbitraire de la condition initiale
      Allocate(rk(Size(b)))
      rk = b - Matmul(A,xk) !Erreur initiale
      k = 0 !Initialisation du garde-fou
      beta = sqrt(Dot_Product(rk,rk))
      Print*,"Le résidu à l'itération ",k," vaut : ",beta
      Allocate(list(kmax))
      list(1) = beta

      !Instructions
      Do While ((beta > eps) .and. (k < kmax))
        zk = Matmul(A,rk)
        alpha = Dot_Product(rk,zk) / Dot_Product(zk,zk)
        xk = xk + alpha * rk
        rk = rk - alpha * zk
        k = k + 1
        beta = sqrt(Dot_Product(rk,rk))
        Print*,"Le résidu à l'itération ",k," vaut : ",beta
        list(k) = beta
      End Do

      !Sortie
      Allocate(X(Size(xk)))
      X = xk
      list = Reshape(list,(/k/))
      Call SaveSol(list)

      !Arret du programme si tolérence non atteinte
      If (k > kmax) Then
        Print*,"Tolérence non atteinte :", sqrt(Dot_Product(rk,rk))
        Print*,"Arret du programme"
        Print*,"!--------------------------------------------------!"
        Stop
      End If
    End Function ResMin

    !Algorithme de décomposition d'un matrice A en Matrice Q orthogonale et R triangulaire supèrieure
    !Méthode de Givens
    Subroutine Decomp_QR(A,Q,R)
      !Déclaration des variables locales
      Real(PR),Dimension(:,:),Intent(In)::A
      Real(PR),Dimension(:,:),Allocatable,Intent(Out)::Q,R

      Real(PR),Dimension(:,:),Allocatable::G
      Real(PR),Dimension(:),Allocatable::x
      Real(PR)::c, s
      Integer::i, j, k, m, n

      !Instructions
      m = Size(A,1); n = Size(A,2)
      Allocate(Q(m,m))
      Q = 0._PR
      Do k = 1, m
        Q(k,k) = 1._PR
      End Do
      Allocate(R(m,n))
      R = A

      Allocate(G(m,m))
      Allocate(x(m))
      Do j = 1, n
        Do i = m, j+1, -1
          x = R(:,j)
          If (sqrt((x(i-1)**2) + (x(i)**2)) > 0) Then
            c = x(i-1) / sqrt((x(i-1)**2) + (x(i)**2))
            s = -x(i) / sqrt((x(i-1)**2) + (x(i)**2))
            G = 0._PR
            Do k = 1, m
              G(k,k) = 1._PR
            End Do
            G(i-1,i-1) = c ; G(i-1,i) = s
            G(i, i-1) = -s ; G(i,i) = c
            Q = Matmul(Q,G)
            R = Matmul(Transpose(G),R)
          End If
        End Do
      End Do

    End Subroutine Decomp_QR

    !Algorithme de résolution de système QR
    Function ResolvQR(A,b) Result(X)
      !Déclaration des variables locales
      Real(PR),Dimension(:,:),Intent(In)::A
      Real(PR),Dimension(:),Intent(In)::b
      Real(PR),Dimension(:),Allocatable::X

      Real(PR)::sum
      Integer::i, k, n

      !Instructions
      n = Size(b)
      Allocate(X(n))

      !Algorithme de remonté seulement puisque R est triangulaire supèrieure
      X(n) = b(n) / A(n,n)
      Do i = n-1, 1, -1
        sum = 0._PR
        Do k = i+1, n
          sum = sum + A(i,k) * X(k)
        End Do
        X(i) = (b(i) - sum) / A(i,i)
      End Do

    End Function ResolvQR

    !Algorithme d'Arnoldi
    Subroutine Arnoldi(A,x0,V,H)
      !Déclaration des variables locales
      Real(PR),Dimension(:,:),Allocatable,Intent(In)::A
      Real(PR),Dimension(:),Allocatable,Intent(In)::x0
      Real(PR),Dimension(:,:),Allocatable,Intent(Out)::V, H

      Integer::m, i, j
      Real(PR),Dimension(:),Allocatable::wj

      !Initialisation
      m = taille_K  !Taille arbitraire de l'espace de Krylov
      Allocate(V(Size(A,2),m+1)) ; V(:,1) = (1 / sqrt(Dot_Product(x0,x0))) * x0
      Allocate(H(m+1,m)); H = 0._PR
      Allocate(wj(Size(A,2)))

      !Instructions
      Do j = 1, m
        wj = Matmul(A,V(:,j))
        Do i = 1, j
          H(i,j) = Dot_Product(wj,V(:,i))
          wj = wj - H(i,j) * V(:,i)
        End Do
        H(j+1,j) = sqrt(Dot_Product(wj,wj))
        If (H(j+1,j) == 0) Then
          Exit
        End If
        V(:,j+1) = (1 / H(j+1,j)) * wj
      End Do

    End Subroutine Arnoldi

    !Algorithme de résolution FOM
    Function FOM(A,b) Result(X)
      !Déclaration des variables locales
      Real(PR),Dimension(:,:),Allocatable,Intent(In)::A !Matrice A
      Real(PR),Dimension(:),Allocatable,Intent(In)::b !Vecteur b
      Real(PR),Dimension(:),Allocatable::X !Vecteur X

      Real(PR),Dimension(:,:),Allocatable::Vm, Hm, Qm, QmT, Rm
      Real(PR),Dimension(:),Allocatable::rk, yk, xk, list
      Real(PR)::beta
      Integer::m, k

      !Initialisation
      Allocate(xk(Size(A,2)))
      xk = 1._PR !Choix arbitraire de la condition initiale
      Allocate(rk(Size(b)))
      rk = b - Matmul(A,xk) !Erreur initiale
      beta = sqrt(Dot_Product(rk,rk))
      m = taille_K
      k = 0
      Print*,"Le résidu à l'itération ",k," vaut : ",beta
      Allocate(list(kmax))
      list(1) = beta

      !Instructions
      Do While ((beta > eps) .and. (k < kmax))
        Call Arnoldi(A,rk,Vm,Hm)
        Call Decomp_QR(Hm(1:m,:),Qm,Rm)
        QmT = Transpose(Qm)
        yk = ResolvQR(Rm,beta * QmT(:,1))
        xk = xk + Matmul(Vm(:,1:m),yk)
        rk = -Hm(m+1,m) * yk(m) * Vm(:,m+1)
        beta = sqrt(Dot_Product(rk,rk))
        k = k+1
        Print*,"Le résidu à l'itération ",k," vaut : ",beta
        list(k) = beta
      End Do

      !Sortie
      Allocate(X(Size(xk)))
      X = xk
      list = Reshape(list,(/k/))
      Call SaveSol(list)

      !Arret du programme si tolérence non atteinte
      If (k > kmax) Then
        Print*,"Tolérence non atteinte :", sqrt(Dot_Product(rk,rk))
        Print*,"Arret du programme"
        Print*,"!--------------------------------------------------!"
        Stop
      End If

    End Function FOM

    !Algorithme de la résolution GMRes
    Function GMRes(A,b) Result(X)
      !Déclaration des variables locales
      Real(PR),Dimension(:,:),Allocatable,Intent(In)::A !Matrice A
      Real(PR),Dimension(:),Allocatable,Intent(In)::b !Vecteur b
      Real(PR),Dimension(:),Allocatable::X !Vecteur X

      Real(PR),Dimension(:,:),Allocatable::Vm, Hm, Qm, QmT, Rm
      Real(PR),Dimension(:),Allocatable::rk, yk, xk, gm, list
      Real(PR)::beta
      Integer::m, k

      !Initialisation
      Allocate(xk(Size(A,2)))
      xk = 1._PR !Choix arbitraire de la condition initiale
      Allocate(rk(Size(b)))
      rk = b - Matmul(A,xk) !Erreur initiale
      beta = sqrt(Dot_Product(rk,rk))
      m = taille_K
      k = 0
      Print*,"Le résidu à l'itération ",k," vaut : ",beta
      Allocate(list(kmax))
      list(1) = beta

      !Instructions
      Do While ((beta > eps) .and. (k < kmax))
        Call Arnoldi(A,rk,Vm,Hm)
        Call Decomp_QR(Hm,Qm,Rm)
        QmT = Transpose(Qm)
        gm = beta * QmT(:,1)
        yk = ResolvQR(Rm(1:m,:),gm(1:m))
        xk = xk + Matmul(Vm(:,1:m),yk)
        rk = gm(m+1) * Matmul(Vm,QmT(:,m+1))
        beta = abs(gm(m+1))
        k = k + 1
        Print*,"Le résidu à l'itération ",k," vaut : ",beta
        list(k) = beta
      End Do

      !Sortie
      Allocate(X(Size(xk)))
      X = xk
      list = Reshape(list,(/k/))
      Call SaveSol(list)

      !Arret du programme si tolérence non atteinte
      If (k > kmax) Then
        Print*,"Tolérence non atteinte :", sqrt(Dot_Product(rk,rk))
        Print*,"Arret du programme"
        Print*,"!--------------------------------------------------!"
        Stop
      End If

    End Function GMRes

End Module mod_algorithmes
