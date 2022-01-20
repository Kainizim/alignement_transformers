"""
Programme d'alignement de sequence par paires a partir d'embedding obtenus par des methodes basees sur des transformers
"""

#import-------------------------------------------------------------------------------------------------------------------------------------
import numpy as np
import matplotlib.pyplot as plt
import argparse

#Arguments ligne de commande---------------------------------------------------------------------------------------------------------------
parser = argparse.ArgumentParser(description = "Alignement par paires de 2 sequences a l'aide de 2 embedding obtenu a l'aide de Transformers.")
parser.add_argument('-s1', '--sequence1', type=str, metavar='', help='Fichier fasta de la sequence 1')
parser.add_argument('-s2', '--sequence2', type=str, metavar='', help='Fichier fasta de la sequence 2')
parser.add_argument('-e1', '--embedding1', type=str, metavar='', help='Embedding de la sequence 1')
parser.add_argument('-e2', '--embedding2', type=str, metavar='', help='Embedding de la sequence 2')
parser.add_argument('-g', '--gap_penality', type=int, metavar='', help='penalite de gap')
args = parser.parse_args()

#Parametres---------------------------------------------------------------------------------------------------------------------------------
gap = args.gap_penality
#emb_1 = np.loadtxt("emb/1bina.t5emb")
#emb_2 = np.loadtxt("emb/1lh1.t5emb")
emb_1 = np.loadtxt("emb/" + args.embedding1)
emb_2 = np.loadtxt("emb/" + args.embedding2)


#Fonctions----------------------------------------------------------------------------------------------------------------------------------

#Chargement des fasta----------------------------------
def fasta(fasta_file):
    """
    Recuperation des fichiers fasta
    parametre : fichier fasta
    """
    seq = []
    with open(fasta_file, "r") as filin:
        for line in filin:
            seq.append(line)
        seq = seq[1]
    return seq

#matrice de score (DP)---------------------------------

def dp_matrix(emb1, emb2): #calcul du dot product
    """
    Matrice de dot product entre la sequence 1 et la sequence 2 ==> SCORE

    parametres : 
    emb 1 : embedding de la sequence 1
    emb 2 : embedding de la sequence 2 
    """

    matrix_dp = np.zeros((len(emb1), len(emb2)))
    #print(matrix_dp.shape)
    for i in range(len(emb1)):
        for j in range(len(emb2)):
            matrix_dp[i][j] = np.dot(emb1[i], emb2[j])
    return matrix_dp

#matrice de programmation dynamique--------------------

def matrix_dyna(seq_1, seq_2, matrix_dp): #matrice de programmation dynamique
    """
    Matrice de programmation dynamique

    parametres : 
    seq_1 : fichier fasta de la sequence 1
    seq_2 : fichier fasta de la sequence 2
    matrix_dp : matrice de score (dot product)
    """
    matrix_prog_1 = [0]
    matrix_prog_2 = [0]
    i = 1
    while i != len(seq_1):
        #print(matrix_score[-i])
        matrix_prog_1.append(matrix_prog_1[i-1] + gap)
        i += 1
    i = 1
    while i != len(seq_2):
        matrix_prog_2.append(matrix_prog_2[i-1] + gap)
        i += 1
    #matrix_score = []
    matrix_dyna = np.zeros((len(seq_1),len(seq_2)))
    for i in range(len(seq_1)):
        matrix_dyna[i][0] = matrix_prog_1[i]
    for i in range(len(seq_2)):
        matrix_dyna[0][i] = matrix_prog_2[i]     
    #print(matrix_score.shape)  ----> Matrice de prog dynamique : Ajout de + gap ou + dp_matrix en fonction de maximum
    for i in range(len(seq_1)):
        for j in range(len(seq_2)):
            if i > 0 and j > 0:
                matrix_dyna[i][j] = max(matrix_dyna[i - 1][j] + gap, matrix_dyna[i][j - 1] + gap, matrix_dyna[i - 1][j - 1] + matrix_dp[i-1][j-1])

    
    return matrix_prog_1, matrix_prog_2, matrix_dyna

#Matrice de chemin-------------------------------------------------------------

def chemin_matrix(matrix_dynamique, matrix_dp):
    """
    Recuperation du meilleur chemin a l'aide de flag : (0 = match/ 1 et 2 = gap)
    parametres:
    matrix_dynamique = matre de programation dynamique entre les 2 sequences
    matrix_dp : matrice de dot product ==> SCORE
    """
    chemin = []
    chemin_flag = []

    i = 0
    j = 0

    print("- Taille Sequence 1 (ligne) : " + str(len(matrix_dynamique)-1) + "\n- Taille Sequence 2 (colonne): " + str(len(matrix_dynamique[0])-1))
    #print("\n")
    #print("- Derniere case dela matrice de DP: ", matrix_dp[i-1][j-1])
    #print("\n")

    #max_value = max(matrix[i-1][j] + gap, matrix[i][j-1] + gap, matrix[i-1][j-1] + dp_matrix)
    #print(max_value)

    while i < len(matrix_dynamique)-1 and j < len(matrix_dynamique[0])-1:
        max_value = max(matrix_dynamique[i+1][j], matrix_dynamique[i][j+1], matrix_dynamique[i+1][j+1])
        
        if max_value == matrix_dynamique[i+1][j+1]:
            #print("Correspondance" , i)
            chemin_flag.append(0)
            chemin.append(matrix_dp[i][j])
            i = i + 1
            j = j + 1

        elif max_value == matrix_dynamique[i+1][j]:
            chemin_flag.append(1)
            chemin.append(matrix_dp[i][j-1])
            #print("GAP_BAS")
            i = i + 1
        
        elif max_value == matrix_dynamique[i][j+1]:
            chemin_flag.append(2)
            chemin.append(matrix_dp[i-1][j])
            #print("GAP_DROIT")
            j = j + 1
        else:
            print("ERREUR")
            break

    #print(chemin_flag)
    #print(chemin)

    #Score de l'alignement------------------------------------------------------------
    somme_score = sum(chemin)
    #print("- Score de l'alignement : " + str(somme_score))
    print("\n")
    return chemin, somme_score, chemin_flag

#stockage et affichage de l'alignement-------------------------------------------------------------------------------------------------------

def alignement(seq_1, seq_2, chemin_flag):
    """
    Alignement des sequence en fonction des evenement(match = lettre, gap = '-')
    parametres :
    seq_1 : fichier fasta de la sequence 1
    seq_2 : fichier fasta de la sequence 2
    chemin_flag = liste du chemin a l'aide de flag (0 = match, 1 et 2 = gap)
    """
    print("Sequence 1 = " + seq_1)
    print("Sequence 2 = " + seq_2)
    #print(chemin)
    #print(len(chemin))
    max_seq = max(len(seq_1), len(seq_2))
    #print(max_seq)
    sequence_1 = []
    sequence_2 = []

    for i in range(len(chemin_flag)):
        if i < len(seq_1):
            if chemin_flag[i] == 0: #match
                sequence_1.append(seq_1[i])
            elif chemin_flag[i] == 1: #gap vers le bas donc i = i+1
                sequence_1.append("-")
            elif chemin_flag[i] == 2: #gap vers le droite donc j = j+1
                sequence_1.append(seq_1[i])
        else:
            sequence_1.append("-")

    for i in range(len(chemin_flag)):
        if i < len(seq_2):
            if chemin_flag[i] == 0:
                sequence_2.append(seq_2[i])
            elif chemin_flag[i] == 1: #gap vers le bas donc i = i+1
                sequence_2.append(seq_2[i]) 
            elif chemin_flag[i] == 2:
                sequence_2.append("-")
        else:
            sequence_2.append("-")

    if seq_1 > seq_2:
        sequence_1.remove('\n')
    else:
        sequence_2.remove('\n')

    sequence_1 = "".join(sequence_1)
    sequence_2 = "".join(sequence_2)

    print("ALIGNEMENT PAR PAIRE DES EMBEDDING T5 OBTENUS PAR TRANSFORMERS : \n")
    print("Sequence 1 = " + sequence_1)
    print("Sequence 2 = " + sequence_2)

    return sequence_1, sequence_2


"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
main-------------------------------------------------------------------------------------------------------------------------------------
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

if __name__ == '__main__':

    # 1) recuperation des Fastas---------------------------------------------------------------------------------------------------------
    seq_1 = fasta("fasta/" + args.sequence1)
    seq_2 = fasta("fasta/" + args.sequence2)

    # 2) Recuperation de la matrice de score (DP)----------------------------------------------------------------------------------------
    matrix_dp = dp_matrix(emb_1, emb_2)
    #np.savetxt("dossier/dp_matrix.txt", matrix_dp)


    # 3) Recuperation de la matrice de programmation dynamique---------------------------------------------------------------------------
    matrix_dynamique_gap_1, matrix_dynamique_gap_2, matrix_dynamique = matrix_dyna(seq_1, seq_2, matrix_dp)
    #np.savetxt("dossier/matrice_dynamique.txt", matrix_dynamique)

    # 4) Recuperation du chemin optimal--------------------------------------------------------------------------------------------------
    chemin, score, chemin_flag = chemin_matrix(matrix_dynamique, matrix_dp)
    #print(chemin_flag)
    #print(chemin)

    # 5) Stockage et affichage de l'alignement-------------------------------------------------------------------------------------------
    seq1_align, seq2_align = alignement(seq_1, seq_2, chemin_flag)
    print("\n- Score de l'alignement : " + str(score))
    #np.savetxt("chemin_1bina_1lh1.txt", chemin_score, fmt="%1.17f")

    # 6) histogramme------------------------------------------------------------------------------------------------------------------------

    """

    #affiche l'histogramme permettant de determiner la penalite de gap------------------------------------------------------------------

    print(np.mean(matrix_dp))
    plt.hist(matrix_dp)
    plt.xlabel('dot product')
    plt.ylabel('y')
    plt.title('Histogramme de la repartition du DP : 1bina avec 1a7s')
    plt.xlim(-13, 13)
    plt.ylim(0, 0.07)
    plt.grid(True)
    plt.show()

    #affiche un representation de la matrice de DP afin d'avoir une vision d'ensmble des matchs et des gap

    plt.matshow(matrix_dp)
    plt.show()
    """