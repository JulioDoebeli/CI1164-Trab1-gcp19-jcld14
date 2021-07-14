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
