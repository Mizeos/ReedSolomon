/*
*
*   Thomas MICHELET
*   Gabriel PERON
*
*   Ceci est notre code de Reed Solomon. C'est à la fois un encodeur et un décodeur de Reed Solomon, 
*   L'encodage est systématique et le décodage se fait avec l'algorithme de Berlekamp.
*
*   Nous avons basé notre code sur un code de Simon Rockliff, University of Adelaide publié le 21/09/89
*
*/

//import java.util.*;
//import java.math.*;

public class ReedSolomon {

    /* déclaration des variables gloables utilisées dans notre algorithme */

    public static int mm = 8; // Nombre de bits par symboles GF(256) = GF(2^8)
    public static int nn = 255; // Longeur du code nn=2**mm-1
    public static int tt = 3; // Capacité de correction (ici 3)
    public static int kk = 249; // Dimension du code kk = nn-2*tt

    public static int [] pp= {1,1,0,0,1,0,0,0,0}; // coefficients du polynome primitif  0 à la fin en mode yolo
    public static int [] alpha_to = new int[nn+1];
    public static int [] index_of = new int[nn+1];
    public static int [] gg= new int [nn-kk+1]; // Definition des tableaux de coefficients de polynome
    public static int [] recd = new int[nn];
    public static int [] data = new int[kk];
    public static int [] bb = new int[nn-kk]; // Definition des tableaux de coefficients de polynome

    public static void main(String[] args) {
        int i;
        /* génération du champ de Gallois GF(2^8) */
        generate_gf();
        System.out.println("Look-up tables for GF(2** "+mm +")");
        System.out.println("  i   alpha_to[i]  index_of[i]");
        for (i = 0; i <= nn; i++) {
            System.out.println(""+i +"      "+alpha_to[i] +"          "+index_of[i]);
        }
        System.out.print("\n\n");
        /* génération du polynome générateur */
        gen_poly();

        /*
         * Création d'un tableau permettant permettant de stocker le message à encoder
         * sous forme polynomiale. Ses coefficients sont tout d'abord mis à zéro.
         */

        for (i = 0; i < kk; i++) {
            data[i] = 0;
        }

        /* On code ici le mot MATHS dans un tableau grace au code ASCII */
        data[0] = 77;
        data[1] = 41;
        data[2] = 84;
        data[3] = 72;
        data[4] = 83;

        /* encodage du mot */
        encode_rs();

        /* stockage du mot transmis dans recd[] */
        for (i = 0; i < nn - kk; i++)
            recd[i] = bb[i];
        for (i = 0; i < kk; i++)
            recd[i + nn - kk] = data[i];

        /*
         * Modification des données transmises pour forcer une erreur et tester la
         * robustesse de l'algorithme.
         */
        data[nn - nn / 2] = 3;

        /* Indexation de recd[] */
        for (i = 0; i < nn; i++)
            recd[i] = index_of[recd[i]];

        /* decodage du mot */
        decode_rs(); // Le programme renvoie un tableau avec les coefficients du polynome
                     // représentant le code décodé

        /* Affichage de tout ce qui a été obtenu précédemment */

        System.out.print("Results for Reed-Solomon code (n="+nn +", k="+kk +", t="+tt +")\n\n");
        System.out.println("  i  data[i]   recd[i](decoded)   (data, recd in polynomial form)");

        for (i = 0; i < nn - kk; i++) {
            System.out.println(""+i +"    "+bb[i] +"      "+recd[i]);
        }
        for (i = nn - kk; i < nn; i++) {
            System.out.println(""+i +"    "+data[i - nn + kk] +"      "+recd[i]);
        }
    }

    /*
     * fonction génératrice du champ de Gallois avec le polynome primitif afin de
     * pouvoir travailler en binaire
     */
    public static void generate_gf() {
        int i, mask;

        mask = 1;
        alpha_to[mm] = 0;
        for (i = 0; i < mm; i++) {
            alpha_to[i] = mask;
            index_of[alpha_to[i]] = i;
            if (pp[i] != 0)
                alpha_to[mm] ^= mask;
            mask <<= 1;
        }
        index_of[alpha_to[mm]] = mm;
        mask >>= 1;
        for (i = mm + 1; i < nn; i++) {
            if (alpha_to[i - 1] >= mask)
                alpha_to[i] = alpha_to[mm] ^ ((alpha_to[i - 1] ^ mask) << 1);
            else
                alpha_to[i] = alpha_to[i - 1] << 1;
            index_of[alpha_to[i]] = i;
        }
        index_of[0] = -1;
    }

    /* fonction pour obtenir le polynome générateur */
    public static void gen_poly() {
        int i, j;

        gg[0] = 2; // Le premier élément vaut 2 dans GF(2^mm)
        gg[1] = 1; // g(x) = (X+alpha) au début donc gg[1] doit être égal à 1
        for (i = 2; i <= nn - kk; i++) // initialisation du polynome générateur de degrè nn-kk = 6
        {
            gg[i] = 1;
            for (j = i - 1; j > 0; j--)
                if (gg[j] != 0)
                    gg[j] = gg[j - 1] ^ alpha_to[(index_of[gg[j]] + i) % nn];
                else
                    gg[j] = gg[j - 1];
            gg[0] = alpha_to[(index_of[gg[0]] + i) % nn]; // On sait que gg[0] ne peut pas valoir 0
        }
        // Pour des raisons de vitesse on convertit le gg[] en index[]
        for (i = 0; i <= nn - kk; i++)
            gg[i] = index_of[gg[i]];
    }

    /* fonction permettant d'encoder le messsage */
    public static void encode_rs() {

        /*
         * Prend l'ensemble des éléments de data[] et les encode un à un pour sortir
         * bb[] qui est la forme polynomiale du code encodé Ici le mot de code est c(X)
         * = data(X)*X**(nn-kk)+ b(X)
         */

        int i, j;
        int feedback;// Encodage systématique utilisant la division polynomiale

        for (i = 0; i < nn - kk; i++)
            bb[i] = 0; // Initialisation du tableau à 0
        for (i = kk - 1; i >= 0; i--) // Parcours le tableau à l'envers de 248 à 0 soit 249 éléments
        {
            feedback = index_of[data[i] ^ bb[nn - kk - 1]];
            if (feedback != -1) {
                for (j = nn - kk - 1; j > 0; j--)
                    if (gg[j] != -1)
                        bb[j] = bb[j - 1] ^ alpha_to[(gg[j] + feedback) % nn];
                    else
                        bb[j] = bb[j - 1];
                bb[0] = alpha_to[(gg[0] + feedback) % nn];
            } else {
                for (j = nn - kk - 1; j > 0; j--)
                    bb[j] = bb[j - 1];
                bb[0] = 0;
            }
            ;
        }
        ;
    };

    /* fonction permettant de décoder le message */
    public static void decode_rs(){

    /* 

    Nous recevons des bits regroupés en symboles de mm-bits dans recd[i],
    i=0...(nn-1), et recd[i] est sous forme d'index (c'est-à-dire sous forme de puissances de alpha).
    Nous calculons d'abord les syndromes 2*tt en substituant alpha**i dans rec(X) et en évaluant, 
    en stockant les syndromes dans s[i], i=1..2tt (laisser s[0] zéro). Ensuite, nous utilisons 
    l'itération de Berlekamp pour trouver le polynôme de localisation d'erreur elp[i]. 
    Si le degré de l'elp est >tt, nous ne pouvons pas corriger toutes les erreurs et nous 
    nous contentons donc de mettre les symboles d'information non corrigés.  Si le degré de
    l'elp est <=tt, nous substituons alpha**i, i=1..n dans l'elp pour obtenir les racines, 
    donc les racines inverses, les numéros d'emplacement d'erreur. Si le nombre d'erreurs 
    localisées n'est pas égal au degré de l'elp, nous avons plus d'erreurs que tt et ne pouvons 
    pas les corriger. Sinon, nous résolvons alors la valeur de l'erreur à l'emplacement de l'erreur
    et nous corrigeons l'erreur. 
    Dans les cas où l'on sait que le nombre d'erreurs est trop important pour être corrigé, 
    les symboles d'information tels qu'ils sont reçus sont émis (l'avantage du codage systématique
    est que l'on espère que certains des symboles d'information seront corrects et que si l'on a 
    de la chance, les erreurs se trouvent dans la partie de parité du mot de code transmis).  
    Bien entendu, ces cas insolubles peuvent être renvoyés comme indicateurs d'erreur à la 
    routine d'appel si on le souhaite. 

     */

    
        int i,j,u,q;
        int [][] elp = new int[nn-kk+2][nn-kk];
        int [] d = new int[nn-kk+2];
        int [] l = new int[nn-kk+2];
        int [] u_lu = new int[nn-kk+2];
        int [] s = new int[nn-kk+1] ;
        int count=0;
        int syn_error=0;
        int [] root = new int[tt];
        int [] loc = new int[tt];
        int [] z = new int[tt+1];
        int [] err = new int[nn];
        int [] reg = new int[tt+1] ;

        /* création des syndromes */ 
        for (i=1; i<=nn-kk; i++)
        { s[i] = 0 ;
            for (j=0; j<nn; j++)
                if (recd[j]!=-1)
                    s[i] ^= alpha_to[(recd[j]+i*j)%nn] ;   

            /* conversion des syndromes de la forme polynomiale à une forme d'indexs */     

            if (s[i]!=0)  syn_error=1 ;        
            s[i] = index_of[s[i]] ;
        } ;

        if (syn_error==1)   // si une erreur est détectée, le processus de correction est lancé 

        /* Lancement de l'algorithme de Berlekamp de manière à détecter les erreurs de transmission */

        {
            /* Initialisation des variables */
            d[0] = 0 ;    //forme index        
            d[1] = s[1] ;  //forme index       
            elp[0][0] = 0 ; //forme index  
            elp[1][0] = 1 ;  //forme polynomiale     
            for (i=1; i<nn-kk; i++)
            { elp[0][i] = -1 ;  //forme index  
                elp[1][i] = 0 ;   //forme polynomiale 
            }
            l[0] = 0 ;
            l[1] = 0 ;
            u_lu[0] = -1 ;
            u_lu[1] = 0 ;
            u = 0 ;

            do
            {
                u++ ;
                if (d[u]==-1)
                { l[u+1] = l[u] ;
                    for (i=0; i<=l[u]; i++)
                    {  elp[u+1][i] = elp[u][i] ;
                        elp[u][i] = index_of[elp[u][i]] ;
                    }
                }
                else

                { q = u-1 ;
                    while ((d[q]==-1) && (q>0)) q-- ;

                    if (q>0) // trouver le premier d[q] différent de 0
                    { j=q ;
                        do
                        { j-- ;
                            if ((d[j]!=-1) && (u_lu[q]<u_lu[j]))
                                q = j ;
                        }while (j>0) ;
                    } ;

                    //On trouve un q tel que d[u]différent de 0 et u_lu[q] est maximum afin de conserver les nouveaux degrès
                    if (l[u]>l[q]+u-q)  l[u+1] = l[u] ;
                    else  l[u+1] = l[q]+u-q ;

                    // crée un nouveau elp(x)
                    for (i=0; i<nn-kk; i++)    elp[u+1][i] = 0 ;
                    for (i=0; i<=l[q]; i++)
                        if (elp[q][i]!=-1)
                            elp[u+1][i+u-q] = alpha_to[(d[u]+nn-d[q]+elp[q][i])%nn] ;
                    for (i=0; i<=l[u]; i++)
                    { elp[u+1][i] ^= elp[u][i] ;
                        elp[u][i] = index_of[elp[u][i]] ;   // Pour convertir l'ancienne valeur de l'elp en index 
                    }
                }
                u_lu[u+1] = u-l[u+1] ;

                //Distance (u+1)
                if (u<nn-kk)    
                {
                    if (s[u+1]!=-1)
                        d[u+1] = alpha_to[s[u+1]] ;
                    else
                        d[u+1] = 0 ;
                    for (i=1; i<=l[u+1]; i++)
                        if ((s[u+1-i]!=-1) && (elp[u+1][i]!=0))
                            d[u+1] ^= alpha_to[(s[u+1-i]+index_of[elp[u+1][i]])%nn] ;
                    d[u+1] = index_of[d[u+1]] ; // On met d[u+1] en index     
                }
            } while ((u<nn-kk) && (l[u+1]<=tt)) ;

            u++ ;
            if (l[u]<=tt)     // L'erreur peut être corrigée    
            {
                // On met elp sous la forme index
                for (i=0; i<=l[u]; i++)   elp[u][i] = index_of[elp[u][i]] ;

                // On cherche les racines polynomiales pour la localisation de l'erreur 
                for (i=1; i<=l[u]; i++)
                    reg[i] = elp[u][i] ;
                count = 0 ;
                for (i=1; i<=nn; i++)
                {  q = 1 ;
                    for (j=1; j<=l[u]; j++)
                        if (reg[j]!=-1)
                        { reg[j] = (reg[j]+j)%nn ;
                            q ^= alpha_to[reg[j]] ;
                        } ;
                    if (q==0)        //Stock l'indice de la racine et et la location de l'erreur 
                    { root[count] = i;
                        loc[count] = nn-i ;
                        count++ ;
                    };
                } ;
                if (count==l[u])     // Pas de racines = le degré de l'elp donc <= tt erreurs 
                {
                    /* forme polynomiale de z(x) */
                    for (i=1; i<=l[u]; i++)        /* Z[0] = 1 toujours */
                    { if ((s[i]!=-1) && (elp[u][i]!=-1))
                            z[i] = alpha_to[s[i]] ^ alpha_to[elp[u][i]] ;
                        else if ((s[i]!=-1) && (elp[u][i]==-1))
                            z[i] = alpha_to[s[i]] ;
                        else if ((s[i]==-1) && (elp[u][i]!=-1))
                            z[i] = alpha_to[elp[u][i]] ;
                        else
                            z[i] = 0 ;
                        for (j=1; j<i; j++)
                            if ((s[j]!=-1) && (elp[u][i-j]!=-1))
                                z[i] ^= alpha_to[(elp[u][i-j] + s[j])%nn] ;
                        z[i] = index_of[z[i]] ;  // Met sous forme index  
                    } ;

                    /* Evalue les erreurs des locations données par loc[i] */
                    for (i=0; i<nn; i++)
                    { err[i] = 0 ;
                        if (recd[i]!=-1)        /* conversion de recd[] à la forme polynomiale */
                            recd[i] = alpha_to[recd[i]] ;
                        else  recd[i] = 0 ;
                    }
                    for (i=0; i<l[u]; i++)    /* execute le numérateur d'erreur en premier */
                    { err[loc[i]] = 1;       /* accounts for z[0] */
                        for (j=1; j<=l[u]; j++)
                            if (z[j]!=-1)
                                err[loc[i]] ^= alpha_to[(z[j]+j*root[i])%nn] ;
                        if (err[loc[i]]!=0)
                        { err[loc[i]] = index_of[err[loc[i]]] ;
                            q = 0 ;          /* le dominateur du terme d'erreur */
                            for (j=0; j<l[u]; j++)
                                if (j!=i)
                                    q += index_of[1^alpha_to[(loc[j]+root[i])%nn]] ;
                            q = q % nn ;
                            err[loc[i]] = alpha_to[(err[loc[i]]-q+nn)%nn] ;
                            recd[loc[i]] ^= err[loc[i]] ;  //recd[i] doit être sous forme polynominale
                        }
                    }
                }
                else    /* Pas de racine != degrès d'elp => >tt les erreurs sont trop nombreuses pour être résolues*/
                    for (i=0; i<nn; i++)  // retourne une notification d'erreur si désiré      
                        if (recd[i]!=-1)   //convertit recd[] sous forme poynomiale     
                            recd[i] = alpha_to[recd[i]] ;
                        else  recd[i] = 0 ;     // retourne le code reçu sans le modifier 
            }
            else         // elp a un degré >tt donc ne peut pas être résolu  
                for (i=0; i<nn; i++)      // retourne une notification d'erreur si désiré 
                    if (recd[i]!=-1)        //convertit recd[] sous forme poynomiale  
                        recd[i] = alpha_to[recd[i]] ;
                    else  recd[i] = 0 ;    // retourne le code reçu sans le modifier 
        }
        else      // pas de syndromes différents de zéro, il y a donc absence d'erreurs 
            for (i=0; i<nn; i++)
                if (recd[i]!=-1)        
                    recd[i] = alpha_to[recd[i]] ; // conversion de recd[] sous forme polynomiale
                else  recd[i] = 0 ;
    }
}