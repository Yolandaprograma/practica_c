#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define N 512
#define MAX_ITER 1000
float V1[N], V2[N], V3[N];
float Mat [N][N], MatDD[N][N];

#define Mat(i, j) Mat[i][j]
#define MatDD(i, j) MatDD[i][j]

void PrintVect(float vect[N],int from, int numel){
    int i;
    for (i = from; i<from + numel; i++){
        printf("%.6f  ", vect[i]);
    }
    printf("\n");
}

void PrintRow( float Mat[N][N], int row, int from, int numel){
    int i;
    for(i = from; i < from + numel;i++){
        printf("%.6f  ", Mat[row][i]);
    }
    printf("\n");
}

void MultEscalar(float vect[N], float vectres[N], float alfa){
    int i;
    for (i=0; i < N; i++){
        vectres[i] = vect[i] * alfa;
        //printf("%.6f ", vectres[i]);
    }
}

float Scalar(float vect1[N], float vect2[N]){
    int i;
    float resultat = 0;
	for(i = 0; i < N; i++ ){
        resultat += vect1[i] * vect2[i];
	}
    return resultat;
}

float Magnitude( float vect[N] ){
    int i;
    float sum = 0;
    for(i = 0; i < N; i++){
        sum += vect[i] * vect[i];
    }
    return sqrt(sum);
}

int Ortogonal( float vect1[N], float vect2[N] ){
    if(Scalar(vect1, vect2)== 0){
        return 1;
    } else {
        return 0;
    }
}

void Projection( float vect1[N], float vect2[N], float vectres[N] ){
    int i;
    float scalar_producte = Scalar(vect1, vect2);
    float magnitud_elevat = Magnitude(vect2);
    

    for(i = 0; i < N; i++){
        vectres[i]= (scalar_producte / magnitud_elevat) * vect2[i]; 
    } 
}

//i filas j columnas
float Infininorm(float M[N][N]){
    int i;
    int j;
    float max_sum = 0; 
    
    for (i = 0; i < N; i++){
        float fila_sum = 0;
        for(j = 0; j<N;j++){
            fila_sum += fabs(M[i][j]);
        }
        if (fila_sum > max_sum){
            max_sum = fila_sum;   
        }  
    }
    return max_sum;
}

float Onenorm( float M[N][N] ){
    int i;
    int j;
    float max_sum = 0; 
    
    for (j = 0; j < N; j++){
        float colum_sum = 0;
        for(i = 0; i<N;i++){
            colum_sum += fabs(M[i][j]);
        }
        if (colum_sum > max_sum){
            max_sum = colum_sum;   
        }
    }
    return max_sum;
}

float NormFrobenius( float M[N][N] ){
    int i;
    int j;
    float sum_cuadrados = 0; 
    
    for (i = 0; i < N; i++){
        for(j = 0; j<N;j++){
            sum_cuadrados += (M[i][j]) * (M[i][j]);
        }
    }
    return sqrt(sum_cuadrados);
}

int DiagonalDom( float M[N][N] ){
    int i;
    int j;
    for(i=0;i<N;i++){
        float valor_diagonal = fabs(M[i][i]);
        float sum = 0;

        for(j=0;j<N;j++){
            if(j != i){
                sum += fabs(M[i][j]);
            }
        }
        if (valor_diagonal < sum){
            return 0;
        }
    }
    return 1;
}
void Matriu_x_Vector(float Mat[N][N], float vect[N], float vectres[N] ){
    int i, j;
    for(i=0; i < N; i++){
        vectres[i] = 0;
        for(j=0;j<N;j++){
            vectres[i] += Mat[i][j] * vect[j];
        } 
    }
}

int Jacobi( float Mat[N][N] , float vect[N], float vectres[N], unsigned iter ){
    int i,j;
    float vectres_nou[N];
    for (i=0;i<N;i++){
        float sum = 0;
        for (j=0;j<N;j++){
            if(i != j )
                sum += fabs(Mat[i][j]);
            }
        if (fabs(Mat[i][i]) <= sum) {
            printf("La matriu M no és diagonal dominant, no es pot aplicar Jacobi\n");
            return 0; 
                
        }
    }
    for (unsigned k = 0; k < iter; k++){
        for(i=0; i<N; i++){
            vectres_nou[i] = vect[i];
            for (j = 0; j<N; j++){
                if(i != j ){
                    vectres_nou[i] -= Mat[i][j] * vectres[j];
                }
            }
            vectres_nou[i] /= Mat[i][i];
        }
        for (i=0; i<N; i++){
            vectres[i] = vectres_nou[i];
        }
    }
    return 1;
}

void InitData(){
    int i,j;
    srand(334411);
    for( i = 0; i < N; i++ )
        for( j = 0; j < N; j++ ){
            Mat[i][j]=(((i*j)%3)?-1:1)*(100.0*(rand()/(1.0*RAND_MAX)));
            if ( (abs(i - j) <= 3) && (i != j))
                MatDD[i][j] = (((i*j)%3) ? -1 : 1)*(rand()/(1.0*RAND_MAX));
            else if ( i == j )
                MatDD[i][j]=(((i*j)%3)?-1:1)*(10000.0*(rand()/(1.0*RAND_MAX)));
            else MatDD[i][j] = 0.0;
    }
    for( i = 0; i < N; i++ ){
        V1[i]=(i<N/2)?(((i*j)%3)?-1:1)*(100.0*(rand()/(1.0*RAND_MAX))):0.0;
        V2[i]=(i>=N/2)?(((i*j)%3)?-1:1)*(100.0*(rand()/(1.0*RAND_MAX))):0.0;
        V3[i]=(((i*j)%5)?-1:1)*(100.0*(rand()/(1.0*RAND_MAX)));
    }
}

int main(){
    InitData();
    printf("Apartat A:\n");
    printf("V1 del 0 al 9 i del 256 al 265:\n");
    PrintVect(V1,0,10);
    PrintVect(V1,256,10);
    printf("V2 del 0 al 9 i del 256 al 265:\n");
    PrintVect(V2,0,10);
    PrintVect(V2,256,10);
    printf("V3 del 0 al 9 i del 256 al 265:\n");
    PrintVect(V3,0,10);
    PrintVect(V3,256,10);
    printf("\n");
    printf("Apartat B:\n");
    printf("Mat fila 0 i fila 100 del 0 al 9:\n");
    PrintRow(Mat, 0, 0, 10);
    PrintRow(Mat, 100, 0, 10);
    printf("\n");
    printf("Apartat C:\n");
    printf("MatDD fila 0 del 0 al 9 i fila 100 del 95 al 104:");
    PrintRow(MatDD, 0, 0, 10);
    PrintRow(MatDD, 100, 95, 10);
    //MultEscalar(V1, V2, 2.5);
    printf("\n");    
    printf("Apartat D:\n");
    float uni_f =Infininorm(Mat);
    printf("Infininorma de Mat %.3f \n", uni_f);

    float uni_col = Onenorm(Mat);
    printf("Norma ú de Mat %.3f  \n", uni_col);

    float uni_fro= NormFrobenius(Mat);
    printf("Norma de Frobenius de Mat %.3f \n", uni_fro);

    int res_Diagonal = DiagonalDom(Mat);
    if (res_Diagonal){
        printf("La matriu és diagonal dominant\n");
    }else{
        printf("La matriu no és diagonal dominant\n");
    }
    float dd_uni = Infininorm(MatDD);
    printf("Infininorma de MatDD %.3f \n", dd_uni);

    float uni_col_dd = Onenorm(MatDD);
    printf("Norma ú de MatDD %.3f  \n",uni_col_dd);
    
    float uni_fro_dd= NormFrobenius(MatDD);
    printf("Norma de Frobenius de MatDD %.3f \n", uni_fro_dd);

    int resDD_Diagonal = DiagonalDom(MatDD);
    if (resDD_Diagonal){
        printf("La matriu és diagonal dominant\n");
    }else{
        printf("La matriu no és diagonal dominant\n");
    }
    printf("\n");    
    printf("Apartat E:\n");
    float resultat_1 = Scalar(V1, V2);
    printf("El producte escalar de los vectores V1, V2: %f\n", resultat_1);
    float resultat_2 = Scalar(V1, V3);
    printf("El producte escalar de los vectores V1, V3: %f\n", resultat_2);
    float resultat_3 = Scalar(V2, V3);
    printf("El producte escalar de los vectores V2, V3: %f\n", resultat_3);
    printf("\n");
    printf("Apartat F:\n");
    float magnitud_v1 = Magnitude(V1);
    printf("La magnitud del vector és: %f\n", magnitud_v1);
    float magnitud_v2 = Magnitude(V2);
    printf("La magnitud del vector és: %f\n", magnitud_v2);
    float magnitud_v3 = Magnitude(V3);
    printf("La magnitud del vector és: %f\n", magnitud_v3);
    printf("\n");
    printf("Apartat G:\n");
    if (Ortogonal(V1,V2)){;
        printf("Els vectors V1, V2 són ortogonals.\n");
    } else {
        printf("Els vectors V1, V2 no són ortogonals.\n");
    }
      if (Ortogonal(V1,V3)){;
        printf("Els vectors V1, V3 són ortogonals.\n");
    } else {
        printf("Els vectors V1, V3 no són ortogonals.\n");
    }
      if (Ortogonal(V2,V3)){;
        printf("Els vectors V2, V3 són ortogonals.\n");
    } else {
        printf("Els vectors V2, V3 no són ortogonals.\n");
    }    
    printf("\n");
    printf("Apartat H:\n");
    float V3_resultat[N];
    MultEscalar(V3, V3_resultat, 2.0);
    printf("Elements 0 al 9 del resultat de multiplicar V3 x 2.0 són:\n");
    PrintVect(V3_resultat, 0, 10);
    printf("Elements 256 al 265 del resultat de multiplicar V3 x 2.0 són:\n");
    PrintVect(V3_resultat, 256, 10);
    printf("\n");
    printf("Apartat I:\n");
    float V_Project[N];
    Projection(V2,V3,V_Project);
    printf("Els elements 0 a 9 del resultat de la projecció de V2 sobre V3 són:\n");
    PrintVect(V_Project, 0, 10);
    float Vn_Project[N];
    Projection(V1,V2,Vn_Project);
    printf("Els elements 0 a 9 del resultat de la projecció de V1 sobre V2 són:\n");
    PrintVect(Vn_Project, 0, 10);
    printf("\n");
    printf("Apartat J:\n");
    float vectres[N];
    Matriu_x_Vector(Mat,V2,vectres);
    printf("Els elements 0 a 9 del resultat de la multiplicació de Mat per v2 són:\n");
    PrintVect(vectres, 0, 10);
    printf("\n");
    printf("Apartat K:\n");
    float X[N] = {0};
    if(Jacobi(MatDD, V3, X, 1)){
        printf("Els elements 0 a 9 de la solució (1 iter) del sistema d'equacions són:\n");
        PrintVect(X, 0, 10);
    } 
    if(Jacobi(MatDD, V3, X, 1000)){
        printf("Els elements 0 a 9 de la solució (1000 iters) del sistema d'equacions són:\n");
        PrintVect(X, 0, 10);
    } 

      printf("Mat * X = V3:\n");
    if (Jacobi(Mat, V3, X, 1000)) {
        printf("Els elements 0 a 9 de la solució (1000 iters) del sistema d'equacions són:\n");
        PrintVect(X, 0, 10);
    } 
    
}
