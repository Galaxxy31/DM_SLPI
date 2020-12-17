!--Programme principale
Program main
  Use mod_algorithmes
  Use mod_donnees
  Implicit None

  !Instructions
  !Demande du choix d'un problème
  Print*,"!--------------------------------------------------!"
  Print*,"Choisissez le type d'exo :"
  Print*,"1 - Test sur Matrice 3x3 SDP (Q1)"
  Print*,"2 - Matrice An = In + alpha * BnT Bn (Q2-Q3-Q4)"
  Print*,"3 - Matrice S3RMT3M3 (Q5)"
  Read*,userChoice
  Print*,"!--------------------------------------------------!"
  Select Case(userChoice)
  Case(1) !Test sur Matrice 3x3
      Print*,"Initialisation Matrice"
      sys = "Test"
      Call Init_Test(Mat_A,Vect_b)
      Print*,"Fin de l'initialisation"
    Case(2) !Matrice An = In + alpha BnT Bn
      Print*,"Initialisation Matrice"
      sys = "An"
      Call Init_An(Mat_A,Vect_b)
      Print*,"Fin de l'initialisation"
    Case(3) !Matrice S3RMT3M3
      Print*,"Initialisation Matrice"
      sys = "S3RMT3M3"
      Print*, "Stationnaire = 0, Instationnaire = 1 :"
      Read*,userChoice
      Call Init_S3RMT3M3(Mat_A,Vect_b,userChoice)
      Print*,"Fin de l'initialisation"
    Case Default
      Print*,"Ce choix n'est pas disponible, veuillez recommencer."
      Stop
  End Select

  !Demande de choix d'une méthode de résolution
  Print*,"!--------------------------------------------------!"
  Print*,"Choisissez une méthode de résolution :"
  Print*,"1 - Gradient à Pas Optimal (GPO)"
  Print*,"2 - Gradient Conjugué (GC)"
  Print*,"3 - Résidu Minimum"
  Print*,"4 - Full Orthogonalization Method (FOM)"
  Print*,"5 - Generalized Minimal Residual method (GMRes)"
  Read*,userChoice
  Print*,"!--------------------------------------------------!"

  Select Case(userChoice)
  Case(1) !Gradient à Pas Optimal
      Print*,"Début du calcul"
      meth = "GPO"
      Vect_X = GPO(Mat_A,Vect_b)
      Print*,"Fin du calcul"
    Case(2) !Gradient Conjugué
      Print*,"Début du calcul"
      meth = "GC"
      Vect_X = GC(Mat_A,Vect_b)
      Print*,"Fin du calcul"
    Case(3) !Résidu Minimum
      Print*,"Début du calcul"
      meth = "RM"
      Vect_X = ResMin(Mat_A,Vect_b)
      Print*,"Fin du calcul"
    Case(4) !FOM
      Print*,"Début du calcul"
      meth = "FOM"
      Vect_X = FOM(Mat_A,Vect_b)        !FOM initiale
      !Vect_X = FOM_bis(Mat_A,Vect_b)   !FOM si A est SDP
      Print*,"Fin du calcul"
    Case(5) !GMRes
      Print*,"Début du calcul"
      meth = "GMRes"
      Vect_X = GMRes(Mat_A,Vect_b)
      Print*,"Fin du calcul"
    Case Default
      Print*,"Ce choix n'est pas disponible, veuillez recommencer."
      Stop
  End Select
  Print*,"!--------------------------------------------------!"

  !Affichage de la solution
  Print*,"La solution du système linéaire est :"
  Print*,"Vect_X = ",Vect_X
  Print*,"!--------------------------------------------------!"

  !Vérification de la solution
  !Call Validation(Mat_A,Vect_X,Vect_b)

  Print*,"Fin du programme"
  Print*,"!--------------------------------------------------!"

End Program main
