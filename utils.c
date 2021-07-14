/*
 * utils.c
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

#include "utils.h"

/*!  
	\brief Retorna UnixTimeStamp em milisegundos,     
    \details Retorna número de milisegundos que se passaram des do dia 1 de janeiro de 1970 às 00:00:00 hora do meridiano de greenwich.
    \param void, nenhum parametro deve ser passado.
    \return Retorna tempo em milisegundos;
    tempo = timestamp();
    <trecho de programa do qual se deseja medir tempo>
    tempo = timestamp() - tempo;
*/
double timestamp(void){
  struct timeval tp;
  gettimeofday(&tp, NULL);
  return((double)(tp.tv_sec*1000.0 + tp.tv_usec/1000.0));
}

