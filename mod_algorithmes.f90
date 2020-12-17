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
      Print*,"Vérification de la solution :"
      Allocate(Diff(Size(b)))
      Diff = Matmul(A,X)-b
      Do i = 1,Size(b)
        If (abs(Diff(i)) > (10._PR**(-5))) Then !Vérification avec une tolérance de 10**-5
          Print*,"Le résultat du système est incorrect"
          Print*,"AX = ",Matmul(Mat_A,Vect_X)
          Print*,""
          Print*,"b = ",Vect_b
          Print*,"!--------------------------------------------------!"
          Print*,"Fin du Programme"
          Print*,"!--------------------------------------------------!"
          Stop
        End If
      End Do
      Print*,"AX = ",Matmul(Mat_A,Vect_X)
      Print*,""
      Print*,"b = ",Vect_b
      Print*,"Le résultat du système est correct"
      Print*,"!--------------------------------------------------!"

      Deallocate(Diff)
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

    !Algorithme d'affichage du produit scalaire des vecteurs de Vm
    Subroutine test_orthogonalite(V)
      Real(PR),Dimension(:,:),Intent(In)::V

      Integer::i, j

      Print*,"Produit scalaire entre les vecteurs :"
      Do i = 1, Size(V,2)-1
        Do j = 1, Size(V,2)-1
          If (i /= j) Then !Test tout les produits scalaires sauf ceux donnant la norme au carré
            Print*, i, j, "Produit scalaire : ", Dot_Product(V(:,i),V(:,j))
          End If
        End Do
      End Do
      Print*,"Fin produit scalaire"
    End Subroutine test_orthogonalite

    !Algorithme d'enregistrement des résidus et du temps de calcul
    Subroutine SaveSol(Liste,time)
      !Déclaration des variables locales
      Real(PR),Dimension(:),Intent(In)::Liste, time

      Integer::k
      Character(len=30)::n_string

      !Instructions
      Write(n_string,*)taille_n
      Open(Unit=100, file=trim(adjustl(sys))//'_'//trim(adjustl(meth))//'_'//trim(adjustl(n_string))//'.dat')
      Do k = 1, Size(Liste)
        Write(100,*)k, Liste(k), time(k)
      End Do
      Close(100)
    End Subroutine SaveSol

    !Algorithme lecture fichier Matrice S3RMT3M3
    Subroutine Lecture(A,b)
      !Déclaration des variables locales
      Real(PR),Dimension(:,:),Allocatable,Intent(InOut)::A
      Real(PR),Dimension(:),Allocatable,Intent(InOut)::b

      Integer::n_line, dim1, dim2, i, j, io
      Character(len=500)::temp

      !Instructions
      n_line = 1
      Open(file = "s3rmt3m3.mtx",Unit = 100)
      Do While (Is_Iostat_end(io) .eqv. .False.)
        If (n_line < 36) Then !Saut des premières lignes contenant des informations sur la matrice
          Read(100,*,iostat = io)temp
          n_line = n_line + 1
        Else If (n_line == 36) Then !Lecture dimension + Initialisation matrice A et vecteur b
          Read(100,*, iostat = io)dim1, dim2, temp
          taille_n = dim1
          Allocate(A(dim1,dim2))
          A = 0._PR
          Allocate(b(dim2))
          b = 1._PR
          n_line = n_line + 1
        Else  !Ecriture des coeffs dans la matrice
          Read(100,*, iostat = io)i, j, A(i,j)
          A(j,i) = A(i,j)
        End If
      End Do
      Close(100)
    End Subroutine Lecture

    !----------------------------------------------------------------------!

    !Initialisation pour test sur Matrice 3x3 SDP
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

    !Initialisation pour Matrice An = In + alpha * BnT Bn
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

    !Initialisation pour Matrice S3RMT3M3
    Subroutine Init_S3RMT3M3(A,b,choice)
      !Déclaration des variables locales
      Real(PR),Dimension(:,:),Allocatable,Intent(InOut)::A
      Real(PR),Dimension(:),Allocatable,Intent(InOut)::b
      Integer,Intent(In)::choice

      Integer::i

      !Instructions
      Call Lecture(A,b)   !initialisation de la matrice A et du vecteur b
      If (choice == 1) Then   !Changement de la matrice si cas Instationnaire
        A = dt * A
        Do i = 1, Size(A,1)
          A(i,i) = 1 + A(i,i)
        End Do
      Else If (choice > 1) Then
        Print*,"Ce choix n'est pas disponible, veuillez recommencer."
        Stop
      End If
    End Subroutine Init_S3RMT3M3

    !----------------------------------------------------------------------!

    !Algorithme du Gradient à pas optimal
    Function GPO(A,b) Result(X)
      !Déclaration des variables locales
      Real(PR),Dimension(:,:),Intent(In)::A !Matrice A
      Real(PR),Dimension(:),Intent(In)::b !Vecteur b
      Real(PR),Dimension(:),Allocatable::X !Vecteur X

      Real(PR),Dimension(:),Allocatable::xk, rk, zk, list, time
      Integer::k
      Real(PR)::alpha,beta
      Real(PR)::t_init, t_iteration

      !Initialisation
      Allocate(time(kmax))
      Allocate(xk(Size(A,2)))
      xk = 1._PR !Choix arbitraire de la condition initiale
      Allocate(rk(Size(b)))
      Call CPU_TIME(t_init)
      rk = b - Matmul(A,xk) !Erreur initiale
      beta = sqrt(Dot_Product(rk,rk))
      Call CPU_time(t_iteration)
      time(1) = t_iteration - t_init
      k = 0 !Initialisation du garde-fou
      If (Aff .eqv. .True.) Then
        Print*,"Le résidu à l'itération ",k," vaut : ",beta, " après : ", time(k+1), " s"
      End If
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
        Call CPU_time(t_iteration)
        time(k) = t_iteration - t_init
        If (Aff .eqv. .True.) Then
          Print*,"Le résidu à l'itération ",k," vaut : ",beta, " après : ", time(k), " s"
        End If
        list(k) = beta
      End Do

      !Sortie
      Allocate(X(Size(xk)))
      X = xk
      list = Reshape(list,(/k/))
      time = Reshape(time,(/k/))
      Call SaveSol(list,time)

      !Arret du programme si tolérence non atteinte
      If (k == kmax) Then
        Print*,"Tolérence non atteinte : ", sqrt(Dot_Product(rk,rk))
        Print*,"Arret du programme"
        Print*,"!--------------------------------------------------!"
        Stop
      End If

      Deallocate(xk); Deallocate(rk); Deallocate(zk); Deallocate(list); Deallocate(time)

    End Function GPO

    !Algorithme du Gradient Conjugué
    Function GC(A,b) Result(X)
      !Déclaration des variables locales
      Real(PR),Dimension(:,:),Intent(In)::A !Matrice A
      Real(PR),Dimension(:),Intent(In)::b !Vecteur b
      Real(PR),Dimension(:),Allocatable::X !Vecteur X

      Real(PR),Dimension(:),Allocatable::xk, rk, rk1, zk, p, list, time
      Real(PR)::alpha, beta, gamma
      Integer::k
      Real(PR)::t_init, t_iteration

      !Initialisation
      Allocate(time(kmax))
      Allocate(xk(Size(A,2)))
      xk = 1._PR !Choix arbitraire de la condition initiale
      Allocate(rk(Size(b)))
      Call CPU_TIME(t_init)
      rk = Matmul(A,xk) - b !Erreur initiale
      p = -rk
      beta = sqrt(Dot_Product(rk,rk))
      Call CPU_TIME(t_iteration)
      time(1) = t_iteration - t_init
      k = 0 !Initialisation du garde-fou
      If (Aff .eqv. .True.) Then
        Print*,"Le résidu à l'itération ",k," vaut : ",beta, " après : ", time(k+1), " s"
      End If

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
        Call CPU_time(t_iteration)
        time(k) = t_iteration - t_init
        If (Aff .eqv. .True.) Then
          Print*,"Le résidu à l'itération ",k," vaut : ",beta, " après : ", time(k), " s"
        End If
        list(k) = beta
      End Do

      !Sortie
      Allocate(X(Size(xk)))
      X = xk
      list = Reshape(list,(/k/))
      time = Reshape(time,(/k/))
      Call SaveSol(list,time)

      !Arret du programme si tolérence non atteinte
      If (k == kmax) Then
        Print*,"Tolérence non atteinte :", sqrt(Dot_Product(rk,rk))
        Print*,"Arret du programme"
        Print*,"!--------------------------------------------------!"
        Stop
      End If

      Deallocate(xk); Deallocate(rk); Deallocate(p); Deallocate(list); Deallocate(time)

    End Function GC

    !Algorithme du Résidu Minimum
    Function ResMin(A,b) Result(X)
      !Déclaration des variables locales
      Real(PR),Dimension(:,:),Intent(In)::A !Matrice A
      Real(PR),Dimension(:),Intent(In)::b !Vecteur b
      Real(PR),Dimension(:),Allocatable::X !Vecteur X

      Real(PR),Dimension(:),Allocatable::xk, rk, zk, list, time
      Real(PR)::alpha,beta
      Integer::k
      Real(PR)::t_init, t_iteration

      !Initialisation
      Allocate(time(kmax))
      Allocate(xk(Size(A,2)))
      xk = 1._PR !Choix arbitraire de la condition initiale
      Allocate(rk(Size(b)))
      Call CPU_TIME(t_init)
      rk = b - Matmul(A,xk) !Erreur initiale
      k = 0 !Initialisation du garde-fou
      beta = sqrt(Dot_Product(rk,rk))
      Call CPU_TIME(t_iteration)
      time(1) = t_iteration - t_init
      If (Aff .eqv. .True.) Then
        Print*,"Le résidu à l'itération ",k," vaut : ",beta, " après : ", time(k+1), " s"
      End If
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
        Call CPU_time(t_iteration)
        time(k) = t_iteration - t_init
        If (Aff .eqv. .True.) Then
          Print*,"Le résidu à l'itération ",k," vaut : ",beta, " après : ", time(k), " s"
        End If
        list(k) = beta
      End Do

      !Sortie
      Allocate(X(Size(xk)))
      X = xk
      list = Reshape(list,(/k/))
      time = Reshape(time,(/k/))
      Call SaveSol(list,time)

      !Arret du programme si tolérence non atteinte
      If (k == kmax) Then
        Print*,"Tolérence non atteinte :", sqrt(Dot_Product(rk,rk))
        Print*,"Arret du programme"
        Print*,"!--------------------------------------------------!"
        Stop
      End If

      Deallocate(xk); Deallocate(rk); Deallocate(zk); Deallocate(list); Deallocate(time)

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

      Deallocate(G); Deallocate(x)

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

      Deallocate(wj)

    End Subroutine Arnoldi

    !Algorithme de résolution FOM
    Function FOM(A,b) Result(X)
      !Déclaration des variables locales
      Real(PR),Dimension(:,:),Allocatable,Intent(In)::A !Matrice A
      Real(PR),Dimension(:),Allocatable,Intent(In)::b !Vecteur b
      Real(PR),Dimension(:),Allocatable::X !Vecteur X

      Real(PR),Dimension(:,:),Allocatable::Vm, Hm, Qm, QmT, Rm
      Real(PR),Dimension(:),Allocatable::rk, yk, xk, list, time
      Real(PR)::beta
      Integer::m, k
      Real(PR)::t_init, t_iteration

      !Initialisation
      Allocate(time(kmax))
      Allocate(xk(Size(A,2)))
      xk = 1._PR !Choix arbitraire de la condition initiale
      Allocate(rk(Size(b)))
      Call CPU_TIME(t_init)
      rk = b - Matmul(A,xk) !Erreur initiale
      beta = sqrt(Dot_Product(rk,rk))
      Call CPU_TIME(t_iteration)
      time(1) = t_iteration - t_init
      m = taille_K
      Print*,"La taille de l'espace de Krylov est : ",m
      k = 0
      If (Aff .eqv. .True.) Then
        Print*,"Le résidu à l'itération ",k," vaut : ",beta, " après : ", time(k+1), " s"
      End If
      Allocate(list(kmax))
      list(1) = beta

      !Instructions
      Do While ((beta > eps) .and. (k < kmax))
        Call Arnoldi(A,rk,Vm,Hm)
        !Call test_orthogonalite(Vm)
        Call Decomp_QR(Hm(1:m,:),Qm,Rm)
        QmT = Transpose(Qm)
        yk = ResolvQR(Rm,beta * QmT(:,1))
        xk = xk + Matmul(Vm(:,1:m),yk)
        rk = -Hm(m+1,m) * yk(m) * Vm(:,m+1)
        beta = sqrt(Dot_Product(rk,rk))
        k = k+1
        Call CPU_time(t_iteration)
        time(k) = t_iteration - t_init
        If (Aff .eqv. .True.) Then
          Print*,"Le résidu à l'itération ",k," vaut : ",beta, " après : ", time(k), " s"
        End If
        list(k) = beta
        Deallocate(Hm); Deallocate(Vm); Deallocate(Qm); Deallocate(QmT); Deallocate(Rm); Deallocate(yk)
      End Do

      !Sortie
      Allocate(X(Size(xk)))
      X = xk
      list = Reshape(list,(/k/))
      time = Reshape(time,(/k/))
      Call SaveSol(list,time)

      !Arret du programme si tolérence non atteinte
      If (k == kmax) Then
        Print*,"Tolérence non atteinte :", sqrt(Dot_Product(rk,rk))
        Print*,"Arret du programme"
        Print*,"!--------------------------------------------------!"
        Stop
      End If

      Deallocate(xk); Deallocate(rk); Deallocate(list)

    End Function FOM

    !Algorithme de la résolution GMRes
    Function GMRes(A,b) Result(X)
      !Déclaration des variables locales
      Real(PR),Dimension(:,:),Allocatable,Intent(In)::A !Matrice A
      Real(PR),Dimension(:),Allocatable,Intent(In)::b !Vecteur b
      Real(PR),Dimension(:),Allocatable::X !Vecteur X

      Real(PR),Dimension(:,:),Allocatable::Vm, Hm, Qm, QmT, Rm
      Real(PR),Dimension(:),Allocatable::rk, yk, xk, gm, list, time
      Real(PR)::beta
      Integer::m, k
      Real(PR)::t_init, t_iteration

      !Initialisation
      Allocate(time(kmax))
      Allocate(xk(Size(A,2)))
      xk = 1._PR !Choix arbitraire de la condition initiale
      Allocate(rk(Size(b)))
      Call CPU_TIME(t_init)
      rk = b - Matmul(A,xk) !Erreur initiale
      beta = sqrt(Dot_Product(rk,rk))
      Call CPU_TIME(t_iteration)
      time(1) = t_iteration - t_init
      m = taille_K
      Print*,"La taille de l'espace de Krylov est : ",m
      k = 0
      If (Aff .eqv. .True.) Then
        Print*,"Le résidu à l'itération ",k," vaut : ",beta, " après : ", time(k+1), " s"
      End If
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
        rk = gm(m+1) * Matmul(Vm,Qm(:,m+1))
        beta = abs(gm(m+1))
        k = k + 1
        Call CPU_time(t_iteration)
        time(k) = t_iteration - t_init
        If (Aff .eqv. .True.) Then
          Print*,"Le résidu à l'itération ",k," vaut : ",beta, " après : ", time(k), " s"
        End If
        list(k) = beta
        Deallocate(Hm); Deallocate(Vm); Deallocate(Qm); Deallocate(QmT); Deallocate(Rm); Deallocate(yk); Deallocate(gm)
      End Do

      !Sortie
      Allocate(X(Size(xk)))
      X = xk
      list = Reshape(list,(/k/))
      time = Reshape(time,(/k/))
      Call SaveSol(list,time)

      !Arret du programme si tolérence non atteinte
      If (k == kmax) Then
        Print*,"Tolérence non atteinte :", sqrt(Dot_Product(rk,rk))
        Print*,"Arret du programme"
        Print*,"!--------------------------------------------------!"
        Stop
      End If

      Deallocate(xk); Deallocate(rk); Deallocate(list); Deallocate(time)

    End Function GMRes

    !----------------------------------------------------------------------!

    !Variation d'Arnoldi si A est SDP
    Subroutine Var_Arnoldi(A,x0,V,H)
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
        If (j==1) Then
          H(1,j) = Dot_Product(wj, V(:,1))
          wj = wj - H(1,j) * V(:,1)
          H(j+1,j) = sqrt(Dot_Product(wj,wj))
          If (H(j+1,j) == 0) Then
            Exit
          End If
          V(:,j+1) = (1 / H(j+1,j)) * wj
        Else
          Do i = j-1, j
            H(i,j) = Dot_Product(wj,V(:,i))
            wj = wj - H(i,j) * V(:,i)
          End Do
          H(j+1,j) = sqrt(Dot_Product(wj,wj))
          If (H(j+1,j) == 0) Then
            Exit
          End If
          V(:,j+1) = (1 / H(j+1,j)) * wj
        End If
      End Do

      Deallocate(wj)
    End Subroutine Var_Arnoldi

    !Algorithme de décomposition de Cholesky
    Subroutine Decomp_Chol(A,L)
      !Déclaration des variables locales
      Real(PR),Dimension(:,:),Intent(In)::A
      Real(PR),Dimension(:,:),Allocatable,Intent(Out)::L

      Integer::n, i, j, k
      Real(PR)::sum1, sum2

      !Initialisation
      n=Size(A,1)
      sum1 = 0._PR; sum2 = 0._PR
      Allocate(L(n,n))
      L = 0._PR

      !Instructions
      Do i = 1, n
        Do k = 1, i-1
          sum1 = sum1 + (L(i,k)**2)
        End Do
        L(i,i) = sqrt(A(i,i) - sum1)
        sum1 = 0._PR
        Do j = i+1, n
          Do k = 1, i-1
            sum2 = sum2 + L(j,k) * L(i,k)
          End Do
          L(j,i) = (A(i,j) - sum2) / L(i,i)
          sum2 = 0._PR
        End Do
      End Do
    End Subroutine Decomp_Chol

    !Algorithme de résolution avec matrice de Cholesky
    Function ResolvChol(L,b) Result(x)
      !Déclaration des variables locales
      Real(PR),Dimension(:,:),Intent(In)::L
      Real(PR),Dimension(:),Intent(In)::b

      Integer::n, k, i
      Real(PR),Dimension(:,:),Allocatable::tL
      Real(PR),Dimension(:),Allocatable::x,y
      Real(PR)::sum1, sum2

      !Initialisation
      n = Size(L,1)
      Allocate(tL(n,n)); Allocate(x(n)); Allocate(y(n))
      tL = Transpose(L)

      !Instructions
      !Descente : Ly = b
      y=0._PR
      y(1) = b(1) / L(1,1)
      Do k = 2, n
        sum1 = 0._PR
        Do i = 1, (k-1)
          sum1 = sum1 + L(k,i) * y(i)
        End Do
        y(k) = (b(k) - sum1) / L(k,k)
        sum1 = 0._PR
      End Do

      !Remontée : tLx=y
      x = 0._PR
      x(n) = y(n) / tL(n,n)
      Do k = n-1, 1, -1
        sum2 = 0._PR
        Do i = k+1, n
          sum2 = sum2 + x(i) *tL(k,i)
        End Do
        x(k) = (y(k) - sum2) / tL(k,k)
        sum2 = 0._PR
      End Do
    End Function ResolvChol

    !Algorithme de résolution FOM avec Arnoldi amélioré et decomposition de Cholesky
    Function FOM_bis(A,b) Result(X)
      !Déclaration des variables locales
      Real(PR),Dimension(:,:),Allocatable,Intent(In)::A !Matrice A
      Real(PR),Dimension(:),Allocatable,Intent(In)::b !Vecteur b
      Real(PR),Dimension(:),Allocatable::X !Vecteur X

      Real(PR),Dimension(:,:),Allocatable::Vm, Hm, Lm
      Real(PR),Dimension(:),Allocatable::rk, yk, xk, list, time, e1
      Real(PR)::beta
      Integer::m, k
      Real(PR)::t_init, t_iteration

      !Initialisation
      Allocate(time(kmax))
      Allocate(xk(Size(A,2)))
      xk = 1._PR !Choix arbitraire de la condition initiale
      Allocate(rk(Size(b)))
      Call CPU_TIME(t_init)
      rk = b - Matmul(A,xk) !Erreur initiale
      beta = sqrt(Dot_Product(rk,rk))
      Call CPU_TIME(t_iteration)
      time(1) = t_iteration - t_init
      m = taille_K
      Print*,"La taille de l'espace de Krylov est : ",m
      k = 0
      If (Aff .eqv. .True.) Then
        Print*,"Le résidu à l'itération ",k," vaut : ",beta, " après : ", time(k+1), " s"
      End If
      Allocate(list(kmax))
      list(1) = beta
      Allocate(e1(Size(b)))
      e1 = 0._PR; e1(1) = 1._PR

      !Instructions
      Do While ((beta > eps) .and. (k < kmax))
        Call Var_Arnoldi(A,rk,Vm,Hm)
        Call Decomp_Chol(Hm(1:m,:),Lm)
        yk = ResolvChol(Lm, beta*e1)
        xk = xk + Matmul(Vm(:,1:m),yk)
        rk = -Hm(m+1,m) * yk(m) * Vm(:,m+1)
        beta = sqrt(Dot_Product(rk,rk))
        k = k+1
        Call CPU_time(t_iteration)
        time(k) = t_iteration - t_init
        If (Aff .eqv. .True.) Then
          Print*,"Le résidu à l'itération ",k," vaut : ",beta, " après : ", time(k), " s"
        End If
        list(k) = beta
        Deallocate(Hm); Deallocate(Vm); Deallocate(Lm); Deallocate(yk)
      End Do

      !Sortie
      Allocate(X(Size(xk)))
      X = xk
      list = Reshape(list,(/k/))
      time = Reshape(time,(/k/))
      Call SaveSol(list,time)

      !Arret du programme si tolérence non atteinte
      If (k == kmax) Then
        Print*,"Tolérence non atteinte :", sqrt(Dot_Product(rk,rk))
        Print*,"Arret du programme"
        Print*,"!--------------------------------------------------!"
        Stop
      End If

      Deallocate(xk); Deallocate(rk); Deallocate(list)

    End Function FOM_bis


End Module mod_algorithmes
