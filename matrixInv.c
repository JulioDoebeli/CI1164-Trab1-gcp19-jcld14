/*
 * matrixInv.c
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

#include <stdio.h>	//leitura e escrita 
#include <string.h>	//manipulacao de string
#include <stdlib.h>	//alocacao dinamica (malloc /free)           
#include <math.h>	// calculos com ponto flutuante 
#include "matriz.h" //funcoes com matrizes
#include "utils.h"	//medidores de tempo

/*	
	args:
	-p faz o pivoteamento parcial antes de traingularizar a matriz
	-o <arquivo> modifica o arquivo de saida
	
	codigos de saida:
	0 saida com sucesso
	1 erro: arquivo de saida nao pode ser encontrado ou nao pode ser criado
	2 erro: falha ao locar memoria ou no preenchimento dos dados 
	3 erro: falha no calculo da inversa (matriz nao admite inversa)
	4 erro: falha no calculo da norma L2 do residuo 
*/
int main (int argc , char **argv){

	int p = 0; 

    FILE *saida = stdout;     //coloca o valor de SAIDA em stdout 

	for(int countparam = 0; countparam < argc ; countparam++){ //analiza os parametros do programa ignorando os que nao se encaixam
        if(!strcmp(argv[countparam],"-o")){// caso a entrada seja "-o" a proxima informacao sera o diretorio de saida 
            saida = fopen(argv[countparam + 1],"w");
            countparam++;
            if (!saida){ //teste para ver se o arquivode saida existe
                perror("ERRO: arquivo de saida nao encontrado ou nao pode ser criado\n");  
                return 1; 
            }
        }    
        if(!strcmp(argv[countparam],"-p")){// caso a entrada seja "-p" o tag do pivoteamento se torna 1 
            p = 1; 
        }    
    }

    //declaracao de variaveis
    Matriz_t *matriz, *inversa;
    real_t Ttriangulo = 0,Tcalx = 0,Tcaly = 0,res = 0;


    do{
    	matriz = lerMatriz(); //aloca e le matriz inicial
    	if(matriz == NULL){
    		perror("erro ao alocar memoria ou inserir dados de 'matriz'\n");
    		return 2;
    	}
    	//aloca memoria para a matriz inversa
    	inversa = alocaMatriz(matriz->n);
    	if(inversa == NULL){
    		perror("erro ao alocar memoria ou inserir dados de 'inversa'\n");
    		return 2;
    	}

    	//encontra a matriz inversao pelo metodo de fatoracao LU 
    	if(decomposicaoLU(matriz,inversa,&Ttriangulo,&Tcalx,&Tcaly,p)){
    		perror("erro ao calcular a inversa\n");
    		return 3;
    	}

    	//descobre a norma L2 do residuo da operacao
    	res = normaL2Residuo(matriz,inversa);
    	if (isnan(res)){
    		perror("erro ao calcular a normaL2\n");
    		return 4;
    	}

    	//imprime os resultados 
    	fprintf(saida,"%i \n",matriz->n);
    	prnMatriz(matriz,saida);
    	fprintf(saida,"# \n");
    	prnMatriz(inversa,saida);
    	fprintf(saida,"###############\n");
    	fprintf(saida,"# Tempo Triangularização: %0.6e ms\n",Ttriangulo);
    	fprintf(saida,"# Tempo cálculo de Y: %0.6e ms\n",Tcaly);
    	fprintf(saida,"# Tempo cálculo de X: %0.6e ms\n",Tcalx);
    	fprintf(saida,"# Norma L2 do Residuo: %0.6e\n",res);
    	fprintf(saida,"################\n");
    	fprintf(saida,"\n");
    	
    	//libera memoria 
    	liberaMatriz(matriz);
    	liberaMatriz(inversa);
    
    }while(!feof(stdin)); //para o laco caso encontre o fim do arquivo (crtl+d para simular o fim de um arquivo no terminal)

    return 0; //return com sucesso
}
