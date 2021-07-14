/*
 * matriz.h
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

#ifndef __MATRIZ_H__
#define __MATRIZ_H__


// Estruturas de Dados
typedef double real_t;

typedef struct {
  int n; // tamanho da matriz quadrada
  real_t **A; // coeficientes
} Matriz_t;

// Alocaçao e desalocação de memória
/*!
  \brief Alocaçao de memória para uma nova estrutura matriz
  \param n tamanho da matriz
  \return ponteiro para matriz. NULL se houve erro de alocação
  */
Matriz_t* alocaMatriz (int n);

/*!
  \brief Liberação de memória
  \param matriz Ponteiro para matriz
  */
void liberaMatriz (Matriz_t *SL);

//parte aritmética 
/*!
  \brief Calcula a norma L2 do resíduo de uma multiplicação entre a matriz e sua inversa
  \details A norma L2 do resíduo de retornada é sqrt(|I-M*M'|), em que I é a identidade, M e M' são matriz e sua inversa. O retorno tambem pode ser descrito como raiz da soma da diferença dos quadrados   
  \param Matriz_t *matriz, Ponteiro para a matriz original
  \param Matriz_t *inversa, Ponteiro para a matriz inversa
  \return Norma L2 do resíduo.
*/
real_t normaL2Residuo(Matriz_t *matriz, Matriz_t *inversa);

/*!
  \brief Decompoe uma matriz em um sistema LU e depois encontra a matriz inversa 
  \param matriz Ponteiro para a matriz original
  \param inversa Ponteiro para a matriz inversa
  \param Ttriangulo tempo para ser feito a triangularizacão da matriz (e talvez o pivoteamento)
  \param Tcaly Tempo para calcular Ly = b 
  \param Tcalx tempo para calcular Ux = y
  \return 0 caso sucesso e não zero em caso de erro
*/
int decomposicaoLU(Matriz_t *matriz,Matriz_t *inversa,real_t *Ttriangulo,real_t *Tcaly,real_t *Tcalx,int p);

// Leitura e impressão de matrizes
/*!
  \brief Leitura de matriz a partir de Entrada padrão (stdin).
  \return ponteiro para nova matriz. NULL se houve erro (leitura ou alocação)
  */
Matriz_t *lerMatriz ();

/*!
  \brief Imprime matriz "SL" no arquivo "saida"
  \param SL Ponteiro para matriz.
  \param saida Ponteiro para arquivo.
  \return void, não há valor de retorno.
*/
void prnMatriz (Matriz_t *SL,FILE *saida);

/*!
  \brief Copia uma matriz em uma nova matriz
  uso nova_matriz = copiaMatriz(matriz_velha);
  \return ponteiro para nova matriz. NULL se houve erro de alocação
*/
Matriz_t* copiaMatriz(Matriz_t *SL);

#endif // __MATRIZ_H__

