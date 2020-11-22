!--Programme principale
Program main
  Use mod_algorithmes
  Use mod_donnees
  Implicit None

  !Initialisation
  Allocate(Mat_A(3,3)); Allocate(Vect_b(3))
  Vect_b(1) = 1._PR; Vect_b(2) = 2._PR; Vect_b(3) = 3._PR

  Mat_A(1,1) = 0._PR; Mat_A(1,2) = 2._PR; Mat_A(1,3) = 2._PR
  Mat_A(2,1) = 2._PR; Mat_A(2,2) = 0._PR; Mat_A(2,3) = -2._PR
  Mat_A(3,1) = 2._PR; Mat_A(3,2) = -2._PR; Mat_A(3,3) = 0._PR

  !Instructions

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
      Vect_X = GPO(Mat_A,Vect_b)
      Print*,"Fin du calcul"
    Case(2) !Gradient Conjugué
      Print*,"Début du calcul"
      Vect_X = GC(Mat_A,Vect_b)
      Print*,"Fin du calcul"
    Case(3) !Résidu Minimum
      Print*,"Début du calcul"
      Vect_X = ResMin(Mat_A,Vect_b)
      Print*,"Fin du calcul"
    Case(4) !FOM
      Print*,"Début du calcul"
      Vect_X = FOM(Mat_A,Vect_b)
      Print*,"Fin du calcul"
    Case(5) !GMRes
      Print*,"Début du calcul"
      Vect_X = GMRes(Mat_A,Vect_b)
      Print*,"Fin du calcul"
    Case Default
      Print*,"Ce choix n'est pas disponible, veuillez recommencer."
      Stop
  End Select

  !Affichage de la solution
  Print*,"La solution du système linéaire est :"
  Print*,"Vect_X = ",Vect_X
  Print*,"!--------------------------------------------------!"

  !Vérification de la solution
  Print*,"Vérification de la solution :"
  Call Validation(Mat_A,Vect_X,Vect_b)
  Print*,"!--------------------------------------------------!"

  Print*,"Fin du programme"
  Print*,"Bonne journée !"
  Print*,"!--------------------------------------------------!"

End Program main
