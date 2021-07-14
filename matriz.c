/*
 * matriz.c
 * 
 * 
 * Copyright 2021 gcp19 <gcp19@inf.ufpr.br> jcld14 <jcld14@inf.ufpr.br>
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 * MA 02110-1301, USA.
 * 
 * 
 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "utils.h"
#include "matriz.h"

/*!
  \brief Calcula a norma L2 do resíduo de uma multiplicação entre a matriz e sua inversa
  \details A norma L2 do resíduo de retornada é sqrt(|I-M*M'|), em que I é a identidade, M e M' são matriz e sua inversa. O retorno tambem pode ser descrito como raiz da soma da diferença dos quadrados   
  \param Matriz_t *matriz, Ponteiro para a matriz original
  \param Matriz_t *inversa, Ponteiro para a matriz inversa
  \return Norma L2 do resíduo.
*/
real_t normaL2Residuo(Matriz_t *matriz, Matriz_t *inversa){
    Matriz_t *C;
    real_t norma;

    C = alocaMatriz(matriz->n);
    for (int i=0; i < matriz->n; ++i)
        for (int j=0; j < matriz->n; ++j)
            for (int k=0; k < matriz->n; ++k){
                C->A[i][j] += matriz->A[i][k] * inversa->A[k][j];
                //printf("C[%i][%i](%f) =+ matriz[%i][%i](%f) * inversa[%i][%i](%f) \n",i,j,C->A[i][j],i,k,matriz->A[i][k],k,j,inversa->A[k][j]);
            }

    //prnMatriz(C,stdout);
    norma = 0;
    for (int i=0; i < matriz->n; ++i)
        for (int j=0; j < matriz->n; ++j){
            if(i == j){
                C->A[i][j] -= 1;
                norma += fabs(C->A[i][j])*fabs(C->A[i][j]);
            }
            else{
                norma += fabs(C->A[i][j])*fabs(C->A[i][j]);
            }
    }
    norma = sqrt(norma);
    liberaMatriz(C);
    return norma;
}

/*!
  \brief Decompoe uma matriz em um sistema LU e depois encontra a matriz inversa.
  \param matriz Ponteiro para a matriz original
  \param inversa Ponteiro para a matriz inversa
  \param Ttriangulo tempo para ser feito a triangularizacão da matriz (e talvez o pivoteamento)
  \param Tcaly Tempo para calcular Ly = b 
  \param Tcalx tempo para calcular Ux = y
  \return 0 caso sucesso e não zero em caso de erro
*/
int decomposicaoLU(Matriz_t *matriz,Matriz_t *inversa,real_t *Ttriangulo,real_t *Tcaly,real_t *Tcalx,int p){
    
    Matriz_t *AUX;
    AUX = copiaMatriz(matriz); 
    if (AUX == NULL)
        return 2; //falha, não foi possivel alocar matriz auxiliar


    int P[matriz->n]; //vetor de permutacão
    int i, j, k, maior; 
    double maxA, absA;
    real_t *ptr;

    for (i = 0; i <= matriz->n; i++)
        P[i] = i; //preenche o vetor de permutacão

    //---------------------------triangularizacão de matriz ----------------//
    *Ttriangulo = timestamp();
    for (i = 0; i < matriz->n; i++) {
        
        maxA = 0.0;
        maior = i;

        for (k = i; k < matriz->n; k++)                   
            if ((absA = fabs(AUX->A[k][i])) > maxA) {  //verifica se o elementoda diagonal e nulo
                maxA = absA;
                maior = k;
            }

        if (maxA == 0) 
            return 1; //falha, determinante da matriz nulo

        //--------------pivoteamento parcial opcional-------------------//
        if (p){                                                         //
            if (maior != i) {                                           //
                //pivoteamento em P                                     //
                j = P[i];                                               //
                P[i] = P[maior];                                        //
                P[maior] = j;                                           //
                                                                        //
                //pivoteamento nas linhas da matriz                     //
                ptr = AUX->A[i];                                        //
                AUX->A[i] = AUX->A[maior];                              //
                AUX->A[maior] = ptr;                                    //
                                                                        //
            }                                                           //
        }                                                               //
        //--------------------------------------------------------------//

        //elimianacao de gauss
        for (j = i+1; j < matriz->n; j++) {
            AUX->A[j][i] /= AUX->A[i][i];
            for (k = i+1; k < matriz->n; k++)
                AUX->A[j][k] -= AUX->A[j][i] * AUX->A[i][k];
        }
    }
    *Ttriangulo = timestamp() - *Ttriangulo;
    //------------------fim da traingularizacao -----------------------------//

    //------------------Calculo de Ly = b------------------------------------//
    *Tcaly = timestamp();
    for (j = 0; j < matriz->n; j++) {
        for (i = 0; i < matriz->n; i++) {
            inversa->A[i][j] = P[i] == j ? 1.0 : 0.0;

            for (k = 0; k < i; k++)
               inversa->A[i][j] -= AUX->A[i][k] *inversa->A[k][j];
        }
    }    
    *Tcaly = timestamp() - *Tcaly;
    //-----------------------------------------------------------------------//

    //------------------Calculo de Ux = y------------------------------------//
    *Tcalx = timestamp();

    for (j = 0; j < matriz->n; j++) {
        for (i = matriz->n-1; i >= 0; i--) {
            for (k = i+1; k < matriz->n; k++)
               inversa->A[i][j] -= AUX->A[i][k] *inversa->A[k][j];

           inversa->A[i][j] /= AUX->A[i][i];
        }
    }
    *Tcalx = timestamp() - *Tcalx;
    //-----------------------------------------------------------------------//

    liberaMatriz(AUX); //libera matriz auxiliar utilizada para preservar a matriz original
    return 0;
}



/*!
  \brief Alocaçao de memória para uma nova estrutura matriz
  \param n tamanho da matriz
  \return ponteiro para matriz. NULL se houve erro de alocação
  */
Matriz_t* alocaMatriz (int n){
    Matriz_t *novo;
    novo = (Matriz_t*) malloc (sizeof(Matriz_t));
    if (novo == NULL)
        return NULL;
    novo->n = n;
    novo->A = malloc (n * sizeof (real_t*)) ;
    if (novo->A == NULL)
        return NULL;
    for (int i = 0; i < n; i++){
        novo->A[i] = (real_t*) malloc (n * sizeof (real_t)) ;  
        if (novo->A[i] == NULL)
        return NULL;
    }
    for ( int i = 0; i < novo->n; i++)
        for (int j = 0; j < novo->n; j++)
            novo->A[i][j] = 0;
    return(novo);    
}

/*!
  \brief Liberação de memória
  \param matriz Ponteiro para matriz
  */
void liberaMatriz (Matriz_t *matriz){
    for (int i = 0 ; i < matriz->n ; i++)
        free(matriz->A[i]);
    free(matriz->A);
    free(matriz);
    return;
}

/*!
  \brief Leitura de matriz a partir de Entrada padrão (stdin).
  \return ponteiro para nova matriz. NULL se houve erro (leitura ou alocação)
  */
Matriz_t *lerMatriz (){
    Matriz_t *novo;
    int n;
    scanf("%i",&n);
    if(n <= 0)
        return NULL;
    novo = alocaMatriz(n);
    if(novo == NULL){
        return NULL;
    }
    for (int i = 0; i < novo->n; i++)
        for (int j = 0; j < novo->n; j++)
            scanf("%lf",&novo->A[i][j]);
    return(novo);
}


/*!
  \brief Imprime matriz "SL" no arquivo "saida"
  \param SL Ponteiro para matriz.
  \param saida Ponteiro para arquivo.
  \return void, não há valor de retorno.
*/
void prnMatriz (Matriz_t *SL,FILE *saida){
    for (int i = 0; i < SL->n; i++){
        for(int j = 0; j < SL->n; j++)
            fprintf(saida,"%0.2e ",SL->A[i][j]);
        fprintf(saida,"\n");
    }  
}

/*!
  \brief Copia uma matriz em uma nova matriz
  uso nova_matriz = copiaMatriz(matriz_velha);
  \return ponteiro para nova matriz. NULL se houve erro de alocação
*/
Matriz_t* copiaMatriz(Matriz_t *SL){
    Matriz_t *novo;
    novo = alocaMatriz(SL->n);
    if (novo == NULL){
        fprintf(stderr,"erro ao alocar a matriz\n");
        return NULL;
    }
    novo->n = SL->n;
    for (int  i = 0; i < novo->n; i++)
        for (int j = 0; j < novo->n; j++)
            novo->A[i][j]= SL->A[i][j]; 
    return(novo);
}
