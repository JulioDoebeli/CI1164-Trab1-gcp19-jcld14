# CI1164-Trab1-gcp19-jcld14

## Alunos
gcp19 (gcp19@inf.ufpr.br) e jcld14 (jcld14@inf.ufpr.br).


## Estruturas de Dados
### real_t
```c
typedef double real_t;
```
### Matriz_t
```c
typedef struct {
  int n; // tamanho da matriz quadrada
  real_t **A; // coeficientes
} Matriz_t;
```
## Funções

### Alocação e Desalocação de matrizes
#### alocaMatriz
```c
/*!
  DESCRIÇÃO Alocaçao de memória para uma nova estrutura matriz.
  PARAMETRO n tamanho da matriz.
  RETORNO   ponteiro para matriz. NULL se houve erro de alocação.
*/
struct Matriz_t* alocaMatriz (int n);
```
#### copiaMatriz
```c
/*!
  DESCRIÇÃO Copia uma matriz em uma nova matriz
  USO       nova_matriz = copiaMatriz(matriz_velha);
  RETORNO   ponteiro para nova matriz. NULL se houve erro de alocação
*/
Matriz_t* copiaMatriz(Matriz_t *SL);
```


#### liberaMatriz

```c
/*!
  DESCRIÇÃO Liberação de memória
  PARAMETRO matriz Ponteiro para matriz
  */
void liberaMatriz (Matriz_t *SL);
```


### Leitura, Impressão de matrizes
#### lerMatriz
```c
/*!
  DESCRIÇÃO Leitura de matriz a partir de Entrada padrão (stdin).
  RETORNO   ponteiro para nova matriz. NULL se houve erro (leitura ou alocação)
*/
Matriz_t *lerMatriz ();
```

### Aritmética de Matrizes
#### normaL2Residuo
```c
/*!
  DESCRIÇÃO Calcula a norma L2 do resíduo de uma multiplicação entre a matriz e sua inversa
  DETALHES  A norma L2 do resíduo de retornada é sqrt(|I-M*M'|), em que I é a identidade, M e M' são matriz e sua inversa.
  DETALHES  O retorno tambem pode ser descrito como raiz da soma da diferença dos quadrados.
  PARAMETRO Matriz_t *matriz, Ponteiro para a matriz original
  PARAMETRO Matriz_t *inversa, Ponteiro para a matriz inversa
  RETORNO   Norma L2 do resíduo.
*/
real_t normaL2Residuo(Matriz_t *matriz, Matriz_t *inversa);
```

#### decomposicaoLU
```c
/*!
  DESCRIÇÃO Decompoe uma matriz em um sistema LU e depois encontra a matriz inversa 
  PARAMETRO matriz Ponteiro para a matriz original
  PARAMETRO inversa Ponteiro para a matriz inversa
  PARAMETRO Ttriangulo tempo para ser feito a triangularizacão da matriz (e talvez o pivoteamento)
  PARAMETRO Tcaly Tempo para calcular Ly = b 
  PARAMETRO Tcalx tempo para calcular Ux = y
  RETORNO   0 caso sucesso e não zero em caso de erro
*/
int decomposicaoLU(Matriz_t *matriz,Matriz_t *inversa,real_t *Ttriangulo,real_t *Tcaly,real_t *Tcalx,int p);
```

Segue abaixo destaque para os nucleos do algoritmo da função de decompusição e fatoração LU (triangularizacão de matriz, pivoteamento, elimianacao de gauss e Calculo de Ly-b). 
```c
int decomposicaoLU(Matriz_t *matriz,Matriz_t *inversa,real_t *Ttriangulo,real_t *Tcaly,real_t *Tcalx,int p){
    // AUX é uma copia da Matriz 
    for (i = 0; i <= matriz->n; i++)
        P[i] = i; //preenche o vetor de permutacão
    //---------------------------triangularizacão de matriz ----------------//
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
    //------------------fim da traingularizacao -----------------------------//

    //------------------Calculo de Ly = b------------------------------------//
    for (j = 0; j < matriz->n; j++) {
        for (i = 0; i < matriz->n; i++) {
            inversa->A[i][j] = P[i] == j ? 1.0 : 0.0;

            for (k = 0; k < i; k++)
               inversa->A[i][j] -= AUX->A[i][k] *inversa->A[k][j];
        }
    }    
    //-----------------------------------------------------------------------//

    /*...*/
}

```
