
//Llibreries necessàries 
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define N 512 //Definició de la mida dels vectors i matrius
#define MAX_ITER 1000 //Nombre màxim d'iteracions per al mètode Jacobi

//Definició dels vectors i matrius globals
float V1[N], V2[N], V3[N];
float Mat [N][N], MatDD[N][N];

//Macros per accedir a les matrius més fàcils
#define Mat(i, j) Mat[i][j]
#define MatDD(i, j) MatDD[i][j]

//Funció per imprimir un vector en un rang d'índex
void PrintVect(float vect[N],int from, int numel){
    int i; //Declarem  i per iterar 
    for (i = from; i<from + numel; i++){ //Bucle que recórrer i per iterar sobre els elements inicial from i final from+numel
        printf("%.6f  ", vect[i]); //Imprimeix els elements amb 6 decimals
    }
    printf("\n"); //Imprimeix un salt de línia per organitzar millor la sortida
}
//Funció per imprimir una fila de la matriu en un rang de columnes
void PrintRow( float Mat[N][N], int row, int from, int numel){
    int i; //Declarem variable i per iterar
    for(i = from; i < from + numel;i++){ //Bucle for per recórrer els elements de la fila seleccionada
        printf("%.6f  ", Mat[row][i]); // Imprimeix els elements amb 6 decimals
    }
    printf("\n");
}
//La funció multiplica cada element del vector vect per el valor alfa i emmagatzema el resultat en el vector vectres.
void MultEscalar(float vect[N], float vectres[N], float alfa){
    int i; //Declarem variable i per iterar sobre els elements del vector
    for (i=0; i < N; i++){ //Bucle per recórrer i per iterar sobre els elements 
        vectres[i] = vect[i] * alfa; //Multiplicació de cada element del vector per l'escala alfa i emmagatzemem el resultat en el vector vectres a la mateixa posició
    }
}
//Funció per calcular el producte escalar entre dos vectors
float Scalar(float vect1[N], float vect2[N]){
    int i; //Declarem la variable i per iterar sobre els elements dels vectors
    float resultat = 0; //Inicialitzem la variable resultat on s'emmatzemarà la suma dels productes 
	for(i = 0; i < N; i++ ){ //Bucle for per recórrer tots els elements dels vectors
        resultat += vect1[i] * vect2[i]; //Per cada element, multipliquem els elements corresponents a vect1 i vect2 i afegim a resultat
	}
    return resultat; //Retornem el valor final del producte escalar
}
//Funció per calcular la magnitud d'un vector
float Magnitude( float vect[N] ){
    int i; //Declarem la variable i per iterar sobre els elements de vector 
    float sum = 0; //Inicialitzem la variable sum per acumular la suma dels quadrats dels elements
    for(i = 0; i < N; i++){ //Bucle for per recórrer tots els elements del vector
        sum += vect[i] * vect[i]; //Suma el quadrat de cada element del vector vect a la variable sum
    }
    return sqrt(sum); //Retornem l'arrel de la suma dels quadrats dels elements del vector
}
//Funció per provocar si dos vectors són ortogonals
int Ortogonal( float vect1[N], float vect2[N] ){
    //Es calcula el producte escalar entre els dos vectors
    if(Scalar(vect1, vect2)== 0){
        //Si el producte escalar és 0, els vectors són ortogonals
        return 1;
    } else {
        //Si el producte escalar no és 0, els vectors no són ortogonals
        return 0;
    }
}
//Funció per calcular la projecció d'un vector sobre un altre vector
void Projection( float vect1[N], float vect2[N], float vectres[N] ){
    int i; //Declarem la variable i per iterar sobre els elements dels vectors
    float scalar_producte = Scalar(vect1, vect2); //Calculem el producte escalar entre els vectors 
    float magnitud_elevat = Magnitude(vect2); //Calculem la magnitud del vector vect2.
    for(i = 0; i < N; i++){ //Bucle for per recórrer tots els elements 
        //Per a cada element i, calculem la projecció de vect1 sobre vect2
        //Es calcula el producte escalar dividit per la magnitud de vect2 i multiplicat per cada component de vect2
        vectres[i]= (scalar_producte / magnitud_elevat) * vect2[i]; 
    } 
}
//Funció per calcular la norma infinita d'una matriu
float Infininorm(float M[N][N]){
    int i; //Declarem variables i per iterar sobre les files de la matriu
    int j; //Declarem variable j per iterar sobre les columnes de la matriu
    float max_sum = 0; //Inicialitzem max_sum per emmagatzema la màxima suma de les files
    
    for (i = 0; i < N; i++){//Bucle per iterar sobre cada fila de la matriu
        float fila_sum = 0; //Iniciliatzem fila_sum per acumular la suma absoluta de la fila actual 
        for(j = 0; j<N;j++){ //Bucle per iterar sobre cada element de la fila i
            fila_sum += fabs(M[i][j]); //Afegim el valor absolut de cada element de la fila a fila_sum
        }
        if (fila_sum > max_sum){ //Compovem si la suma absoluta de la fila actual és major que la màxima trobada fins ara
            max_sum = fila_sum;   //Actualitzem max_sum amb el valor de la suma més gran trobada
        }  
    }
    return max_sum; //Retornem la màxima suma absoluta de les files
}
//Funció per calcular la norma 1 d'una matriu
float Onenorm( float M[N][N] ){
    int i; //Declarem i per iterar sobre les files
    int j; //Declarem j per iterar sobre les columnes de la matriu
    float max_sum = 0; //Inicialitzem max_sum per emmagatzemar la màxima suma de les columnes
    
    for (j = 0; j < N; j++){ //Bucle per iterar cada columna de la matriu
        float colum_sum = 0; //Inicialitzem colum_sum per acumular la suma absoluta de la columna actual
        for(i = 0; i<N;i++){ //Bucle per iterar sobre cada element de la columna j
            colum_sum += fabs(M[i][j]); //Afegim el valor absolut de cada element de la columna a colum_sum
        }
        //Comprovem si la suma absoluta de la columna actual és major que la màxima trobada fins ara
        if (colum_sum > max_sum){
            max_sum = colum_sum;   //Actualitzem max_sum amb el valor de la suma més gran trobada
        }
    }
    return max_sum; //Retornem la màxima suma absoluta de les columnes
}
//Funció per calcular la norma de Frobenius d'una matriu
float NormFrobenius( float M[N][N] ){
    int i; //Declarem variable i per iterar sobre les files
    int j; //Declarem variable j per iterar sobre les columnes
    float sum_cuadrados = 0; //Inicialitzem sum_cuadrados per acumular la suma dels quadrats dels elements de la matriu
    for (i = 0; i < N; i++){ //Bucle per iterar sobre cada fila de la matriu
        for(j = 0; j<N;j++){ //Bucle per iterar sobre cada element de la fila i
            sum_cuadrados += (M[i][j]) * (M[i][j]); //Afegim el quadrat de l'element M[i][j] a la suma total
        }
    }
    return sqrt(sum_cuadrados); //Retornem l'arrel quadrada de la suma dels quadrats de tots els elements de la matriu
}
//Funció per verificar si una matriu és diagonal dominant
int DiagonalDom( float M[N][N] ){
    int i; //Declarem i per iterar sobre les files
    int j; //Declarem j per iterar sobre les columnes
    for(i=0;i<N;i++){ //Iterem sobre cada fila de la matriu 
        float valor_diagonal = fabs(M[i][i]); //Obtenim el valor absolut de l'element de la diagonal M[i][i]
        float sum = 0; //Inicialitzem la variable sum per acumular la suma dels elements fora de la diagonal 
        for(j=0;j<N;j++){ //Iterem sobre cada columna de la fila actual
            if(j != i){ //Si la columna no és la diagonal 
                sum += fabs(M[i][j]); //Afegim el valor absolut de l'element M[i][j] a la suma
            }
        }
        //Si el valor de la diagonal no és més gran que la suma dels altres elements de la fila 
        if (valor_diagonal < sum){ 
            return 0; //Retornem 0, la matriu no és diagonal dominant
        }
    }
    return 1; //Si totes les files compleixen la condició, retornem 1, la matriu és diagonal dominant 
}
//Funció per multiplicar una matriu per un vector
void Matriu_x_Vector(float Mat[N][N], float vect[N], float vectres[N] ){
    int i, j; //Declarem les variable i, j per iterar sobre les files i columnes de la matriu
    for(i=0; i < N; i++){ //Iterem sobre cada fila de la matriu
        vectres[i] = 0; //Inicialitzem l'element i del vector resultat a 0 per acumular el resultat
        for(j=0;j<N;j++){ //Iterem sobre cada columna de la matriu per realitzar la multiplicació de la matriu pel vector 
            vectres[i] += Mat[i][j] * vect[j]; //Sumar el producte de l'element de la matriu Mat[i][j] pel vector [j]
        } 
    }
}
//Funció de Jacobi per a la resolució d'un sistema d'equacions lineals de la forma Ax = b
int Jacobi( float Mat[N][N] , float vect[N], float vectres[N], unsigned iter ){
    int i,j; //Variables per a les iteracions sobres les files i, columnes j
    float vectres_nou[N];  //Vector auxiliar per guardar els nous valors de vectres en cada iteració
    for (i=0;i<N;i++){ //Comprovem que la matriu és diagonal dominant
        float sum = 0;
        //Calculem la suma dels valors absoluts de tots els elements de la fila, excepte la diagonal
        for (j=0;j<N;j++){ 
            if(i != j )
                sum += fabs(Mat[i][j]); //Afegim els valors absoluts de la matriu a la variable sum
            }
        //Comprovem si la matriu és diaogonal dominant    
        if (fabs(Mat[i][i]) <= sum) {
            printf("La matriu M no és diagonal dominant, no es pot aplicar Jacobi\n");
            return 0;  //Si no és diagonal dominant, sortim de la funció amb un valor 0
        }
    }
    //Realitzem el mètode de Jacobi durant un nombre d'iteracions indicat per iter
    for (unsigned k = 0; k < iter; k++){
        //Per cada iteració, calculem els nous valors de vectres per cada element 
        for(i=0; i<N; i++){
            vectres_nou[i] = vect[i]; //Copiem el valor actual de vect a vectres_nou per evitar modificar vectres durant la iteració
            //Calculem el valor de vectres_nou[i] sumant els valor corresponents 
            for (j = 0; j<N; j++){
                if(i != j ){
                    vectres_nou[i] -= Mat[i][j] * vectres[j]; //Restem la multiplicació dels elements de la matriu per vectres
                }
            }
            //Dividim pel valor diagonal per obtenir el nou valor de la solució
            vectres_nou[i] /= Mat[i][i];
        }
        //Copiem els resultats calculats a vectres per a la següent iteració
        for (i=0; i<N; i++){
            vectres[i] = vectres_nou[i];
        }
    }
    return 1; //Retornem 1 per identificar la funció a acabat amb èxit
}

//Funció per inicialitzar les dades de les matrius i vectors amb valors aleatoris
void InitData(){
    int i,j; //Variables d'iteració per recórrer matrius i vectors
    srand(334411); //S'inicialitza el generador de números aleatoris amb fixa 334411 
    //Això garanteix que els resultats siguin sempre els mateixos per a cada execució del programa
    //Inicialització de la matriu Mat amb valors aleatoris
    for( i = 0; i < N; i++ )
        for( j = 0; j < N; j++ ){
        //Els valors de la matriu Mat són aleatoris entre -100 i 100
        //El signe es determina per i*j %3, on el valor serà 1 si és 0, -1 si no ho és
            Mat[i][j]=(((i*j)%3)?-1:1)*(100.0*(rand()/(1.0*RAND_MAX)));
            //Si la distància entre i,j és menor o igual a 3(excepte la diagonal),es posa un valor aleatori
            if ( (abs(i - j) <= 3) && (i != j))
            //La diagonal té valors molt més grans, per això s'assigna un valor molt més alt 
                MatDD[i][j] = (((i*j)%3) ? -1 : 1)*(rand()/(1.0*RAND_MAX));
            else if ( i == j )
                MatDD[i][j]=(((i*j)%3)?-1:1)*(10000.0*(rand()/(1.0*RAND_MAX)));
            //Els valors que no estan a la diagonal ni propers són iguals a 0
            else MatDD[i][j] = 0.0;
    }
    //Inicialització dels vectors, V1,V2,V3 amb valors aleatoris
    for( i = 0; i < N; i++ ){
        //V1 conté valors aleatoris per la primera meitat dels seus elements (i < N/2 )
        V1[i]=(i<N/2)?(((i*j)%3)?-1:1)*(100.0*(rand()/(1.0*RAND_MAX))):0.0;
        //V2 conté valors aleatoris per la segona meitat dels seus elements (i >= N/2)
        V2[i]=(i>=N/2)?(((i*j)%3)?-1:1)*(100.0*(rand()/(1.0*RAND_MAX))):0.0;
        //V3 conté valors aleatoris per tots els elements del vector
        V3[i]=(((i*j)%5)?-1:1)*(100.0*(rand()/(1.0*RAND_MAX)));
    }
}

//Funció principal del programa
int main(){
    //Inicialitza les matrius i vectors amb valors aleatoris
    InitData();
    printf("Apartat A:\n");
    printf("V1 del 0 al 9 i del 256 al 265:\n");
    PrintVect(V1,0,10); //Imprimeix els elements de V1 de l'index 0 al 9
    PrintVect(V1,256,10); //Imprimeix els elements de V1 de líndex 256 al 265
    printf("V2 del 0 al 9 i del 256 al 265:\n");
    PrintVect(V2,0,10); //Imprimeix els elements de V2 de l'índex 0 al 9
    PrintVect(V2,256,10);//Imprimeix els elements de V2 de l'índex 256 al 265
    printf("V3 del 0 al 9 i del 256 al 265:\n");
    PrintVect(V3,0,10); //Imprimeix els elements de V3 de l'índex 0 al 9
    PrintVect(V3,256,10);  //Imprimeix els elements de V3 de l'índex 256 al 265
    printf("\n");

    printf("Apartat B:\n");
    printf("Mat fila 0 i fila 100 del 0 al 9:\n");
    PrintRow(Mat, 0, 0, 10); //Imprimeix la fila 0 de la matriu Mat des de l'índex 0 fins al 9
    PrintRow(Mat, 100, 0, 10); //Imprimeix la fila 100 de la matriu Mat des de l'índex 0 fins al 9
    printf("\n");

    printf("Apartat C:\n");
    printf("MatDD fila 0 del 0 al 9 i fila 100 del 95 al 104:");
    PrintRow(MatDD, 0, 0, 10); //Imprimeix la fila 0 de la matriu MatDD des de l'índex 0 fins al 9
    PrintRow(MatDD, 100, 95, 10); //Imprimeix la fila 100 de la matriu MatDD des de l'índex 0 fins al 9
    printf("\n");

    printf("Apartat D:\n");
    float uni_f =Infininorm(Mat); //Calcula la norma infinita de la matriu Mat
    printf("Infininorma de Mat %.3f \n", uni_f);
    float uni_col = Onenorm(Mat); //Calcula la norma 1 de la matriu Mat
    printf("Norma ú de Mat %.3f  \n", uni_col);
    float uni_fro= NormFrobenius(Mat); //Calcula la norma de Frobenius de la matriu Mat
    printf("Norma de Frobenius de Mat %.3f \n", uni_fro);
    int res_Diagonal = DiagonalDom(Mat); //Comprova si la matriu Mat és diagonal dominant
    if (res_Diagonal){
        printf("La matriu és diagonal dominant\n");
    }else{
        printf("La matriu no és diagonal dominant\n");
    }
    printf("\n");

    float dd_uni = Infininorm(MatDD);  //Calcula la norma infinita de la matriu MatDD
    printf("Infininorma de MatDD %.3f \n", dd_uni);
    float uni_col_dd = Onenorm(MatDD); //Calcula la norma 1 de la matriu MatDD
    printf("Norma ú de MatDD %.3f  \n",uni_col_dd);
    float uni_fro_dd= NormFrobenius(MatDD); //Calcula la norma de Frobenius de la matriu MatDD
    printf("Norma de Frobenius de MatDD %.3f \n", uni_fro_dd);
    int resDD_Diagonal = DiagonalDom(MatDD); //Comprova si la matriu MatDD és diagonal dominant
    if (resDD_Diagonal){
        printf("La matriu és diagonal dominant\n");
    }else{
        printf("La matriu no és diagonal dominant\n");
    }
    printf("\n");    

    printf("Apartat E:\n");
    float resultat_1 = Scalar(V1, V2); //Calcula el producte escalar de V1 i V2
    printf("El producte escalar de los vectores V1, V2: %f\n", resultat_1);
    float resultat_2 = Scalar(V1, V3); //Calcula el producte escalar de V1 i V3
    printf("El producte escalar de los vectores V1, V3: %f\n", resultat_2);
    float resultat_3 = Scalar(V2, V3);  //Calcula el producte escalar de V2 i V3
    printf("El producte escalar de los vectores V2, V3: %f\n", resultat_3);
    printf("\n");

    printf("Apartat F:\n");
    float magnitud_v1 = Magnitude(V1);//Calcula la magnitud de V1
    printf("La magnitud del vector és: %f\n", magnitud_v1);
    float magnitud_v2 = Magnitude(V2); //Calcula la magnitud de V2
    printf("La magnitud del vector és: %f\n", magnitud_v2);
    float magnitud_v3 = Magnitude(V3); //Calcula la magnitud de V3
    printf("La magnitud del vector és: %f\n", magnitud_v3);
    printf("\n");

    printf("Apartat G:\n");
    //Comprovació si els vectors V1 i V2 són ortogonals
    if (Ortogonal(V1,V2)){; 
        printf("Els vectors V1, V2 són ortogonals.\n");
    } else {
        printf("Els vectors V1, V2 no són ortogonals.\n");
    }
    //Comprovació si els vectors V1 i V3 són ortogonals
      if (Ortogonal(V1,V3)){;
        printf("Els vectors V1, V3 són ortogonals.\n");
    } else {
        printf("Els vectors V1, V3 no són ortogonals.\n");
    }
    // Comprovació si els vectors V2 i V3 són ortogonals
      if (Ortogonal(V2,V3)){;
        printf("Els vectors V2, V3 són ortogonals.\n");
    } else {
        printf("Els vectors V2, V3 no són ortogonals.\n");
    }    
    printf("\n");

    printf("Apartat H:\n");
    float V3_resultat[N];
    MultEscalar(V3, V3_resultat, 2.0); //Multiplica el vector V3 per 2.0
    printf("Elements 0 al 9 del resultat de multiplicar V3 x 2.0 són:\n");
    PrintVect(V3_resultat, 0, 10); //Imprimeix els elements de V3_resultat del 0 al 9
    printf("Elements 256 al 265 del resultat de multiplicar V3 x 2.0 són:\n");
    PrintVect(V3_resultat, 256, 10);//Imprimeix els elements de V3_resultat de 256 a 265
    printf("\n");

    printf("Apartat I:\n");
    float V_Project[N]; //Definim un vector per emmagatzemar el resultat de la projecció de V2 sobre V3
    Projection(V2,V3,V_Project); //Cridem la funció Projection per calcular la projecció del vector V2 sobre el vector V3
    printf("Els elements 0 a 9 del resultat de la projecció de V2 sobre V3 són:\n");
    PrintVect(V_Project, 0, 10);  //Imprimeix els elements 0 a 9 del resultat de la projecció
    //// Definim un altre vector per emmagatzemar el resultat de la projecció de V1 sobre V2
    float Vn_Project[N];
    Projection(V1,V2,Vn_Project); //Projecció de V1 sobre V2
    printf("Els elements 0 a 9 del resultat de la projecció de V1 sobre V2 són:\n");
    PrintVect(Vn_Project, 0, 10); //Imprimeix els elements 0 a 9 del resultat de la projecció
    printf("\n");

    printf("Apartat J:\n");
    //// Definim un vector 'vectres' de longitud N per emmagatzemar el resultat de la multiplicació Mat * V2
    float vectres[N];
    //Cridem la funció Matriu_x_Vector per multiplicar la matriu Mat pel vector V2 i emmagatzemar el resultat en vectres
    Matriu_x_Vector(Mat,V2,vectres);
    //Imprimim els primers 10 elements del vector vectres, que conté el resultat de la multiplicació Mat * V2
    printf("Els elements 0 a 9 del resultat de la multiplicació de Mat per v2 són:\n");
    PrintVect(vectres, 0, 10);
    printf("\n");
    
    printf("Apartat K:\n");
    //Definim un vector X de longitud N inicialitzat a zero. Aquest vector emmagatzemarà la solució del sistema d'equacions
    float X[N] = {0};

    //Cridem la funció Jacobi per resoldre el sistema d'equacions MatDD * X = V3, amb 1 iteració.
    //El vector X s'actualitza amb la solució en cada iteració.
    if(Jacobi(MatDD, V3, X, 1)){
        //Si la funció Jacobi retorna 1, vol dir que la solució ha estat trobada després de 1 iteració
        //Imprimim els primers 10 elements de la solució X
        printf("Els elements 0 a 9 de la solució (1 iter) del sistema d'equacions són:\n");
        PrintVect(X, 0, 10);
    } 
    //Cridem la funció Jacobi per resoldre el sistema d'equacions MatDD * X = V3, però amb 1000 iteracions.
    //El vector X s'actualitza després de cada iteració.
    if(Jacobi(MatDD, V3, X, 1000)){
        //Si la funció Jacobi retorna 1, vol dir que la solució ha estat trobada després de 1000 iteracions
        //Imprimim els primers 10 elements de la solució X
        printf("Els elements 0 a 9 de la solució (1000 iters) del sistema d'equacions són:\n");
        PrintVect(X, 0, 10);
    } 
    //Imprimim un missatge per mostrar el resultat de multiplicar Mat per la solució X, que ha de donar el vector V3
    printf("Mat * X = V3:\n");
    //Tornem a cridar la funció Jacobi per resoldre el sistema d'equacions Mat * X = V3, amb 1000 iteracions.
    if (Jacobi(Mat, V3, X, 1000)) {
        //Si la funció Jacobi retorna 1, vol dir que la solució ha estat trobada després de 1000 iteracions
        //Imprimim els primers 10 elements de la solució X
        printf("Els elements 0 a 9 de la solució (1000 iters) del sistema d'equacions són:\n");
        PrintVect(X, 0, 10);
    } 
}
