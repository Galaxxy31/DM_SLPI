!--Module contenant les données de notre problème
Module mod_donnees
  Implicit None

  !Déclaration des variables globales
  Integer, Parameter::PR = 8
  Real(PR), Parameter::pi = 4._PR*Atan(1._PR)
  Real(PR), Parameter::eps = 10._PR**(-15)
  Integer, Parameter::kmax = 1000000

  !Définition du système
  Real(PR),Dimension(:,:),Allocatable::Mat_A
  Real(PR),Dimension(:),Allocatable::Vect_b, Vect_X

  !Choix de l'utilisateur
  Integer::userChoice

  !Taille de la matrice An pour la Q2-Q3-Q4
  Integer,Parameter::taille_n = 200

  !Taille de l'espace de Krylov
  Integer,Parameter::taille_K = 200

  !Variables de nom de fichier
  Character(len=30)::sys, meth

  !Variable d'affichage pour afficher le résidu sur la console
  Logical,Parameter::Aff = .True.


  !Variables définies si besoin pour test

  !Décomposition QR
  Real(PR),Dimension(:,:),Allocatable::Mat_Q, Mat_R

  !Arnoldi
  Real(PR),Dimension(:,:),Allocatable::Mat_V, Mat_H

End Module mod_donnees
