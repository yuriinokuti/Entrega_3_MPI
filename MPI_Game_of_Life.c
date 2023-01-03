#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

#define PRIN_PROC 0

int getNeighbors(int **grid, int i, int j); 
int **grid;
int **newGrid;
int **aux;
#define SRAND_VALUE 1985 // Seed para srand()
#define TAMANHO_TOTAL 2048
#define GERACAO_TOTAL 2000
int bufSnd[TAMANHO_TOTAL];
int bufRcv[TAMANHO_TOTAL]; 

typedef struct {
  int secs;
  int usecs;
} TIME_DIFF;
TIME_DIFF *my_difftime(struct timeval *start, struct timeval *end) {
  TIME_DIFF *diff = (TIME_DIFF *)malloc(sizeof(TIME_DIFF));

  if (start->tv_sec == end->tv_sec) {
    diff->secs = 0;
    diff->usecs = end->tv_usec - start->tv_usec;
  } else {
    diff->usecs = 1000000 - start->tv_usec;
    diff->secs = end->tv_sec - (start->tv_sec + 1);
    diff->usecs += end->tv_usec;
    if (diff->usecs >= 1000000) {
      diff->usecs -= 1000000;
      diff->secs += 1;
    }
  }
  return diff;
}

int getNeighbors(int **grid, int i, int j) {
  int totalNeighbors = 0;
  int linhaMin, colunaMin, linhaMax, colunaMax;

  if (i - 1 < 0) {
    linhaMin = TAMANHO_TOTAL - 1;
  } else {
    linhaMin = i - 1;
  }
  if (i + 1 == TAMANHO_TOTAL) {
    linhaMax = 0;
  } else {
    linhaMax = i + 1;
  }
  if (j - 1 < 0) {
    colunaMin = TAMANHO_TOTAL - 1;
  } else {
    colunaMin = j - 1;
  }
  if (j + 1 == TAMANHO_TOTAL) {
    colunaMax = 0;
  } else {
    colunaMax = j + 1;
  }

  totalNeighbors += grid[linhaMin][j];
  totalNeighbors += grid[linhaMin][colunaMin];
  totalNeighbors += grid[linhaMin][colunaMax];
  totalNeighbors += grid[i][colunaMin];
  totalNeighbors += grid[i][colunaMax];
  totalNeighbors += grid[linhaMax][j];
  totalNeighbors += grid[linhaMax][colunaMin];
  totalNeighbors += grid[linhaMax][colunaMax];
  return totalNeighbors;
}

int inicializarCelula(int **grid) {
  int lin, col;
  // GLIDER
  lin = 1;
  col = 1;
  grid[lin][col + 1] = 1;
  grid[lin + 1][col + 2] = 1;
  grid[lin + 2][col] = 1;
  grid[lin + 2][col + 1] = 1;
  grid[lin + 2][col + 2] = 1;

  // R-pentomino
  lin = 10;
  col = 30;
  grid[lin][col + 1] = 1;
  grid[lin][col + 2] = 1;
  grid[lin + 1][col] = 1;
  grid[lin + 1][col + 1] = 1;
  grid[lin + 2][col + 1] = 1;
  return 10;
}

void atualizaGrid(int numProc, int part, int rank) {
  int companheiro;
  int ini = part * rank;

  for (int i = 0; i < (ini + part); i++) {
    for (int j = 0; j < TAMANHO_TOTAL; j++) {
      companheiro = getNeighbors(grid, i, j);
      if (companheiro == 2 && grid[i][j] == 1) {
        newGrid[i][j] = 1;
      } else {
        if (companheiro == 3) {
          newGrid[i][j] = 1;
        } else {
          newGrid[i][j] = 0;
        }
      }
    }
  }
  aux = grid;
  grid = newGrid;
  newGrid = aux;
}

int contarTotal() {
  int totalCelulas = 0;
  for (int i = 0; i < TAMANHO_TOTAL; i++) {
    for (int j = 0; j < TAMANHO_TOTAL; j++) {
      if (grid[i][j] == 1) {
        totalCelulas++;
      }
    }
  }
  return totalCelulas;
}

void prinProc(int numProc) {
  int i, j, origem, tag = 0, div = 0, rank, gen;
  int part = TAMANHO_TOTAL / numProc;
  TIME_DIFF *time;
  struct timeval start, end;

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  printf("Sistema principal iniciado! [%d]\n", rank);
  gettimeofday(&start, NULL);

  grid = malloc(sizeof(int *) * TAMANHO_TOTAL);
  newGrid = malloc(sizeof(int *) * TAMANHO_TOTAL);
  for (int cont = 0; cont < TAMANHO_TOTAL; cont++) {
    grid[cont] = (int *)calloc(TAMANHO_TOTAL, sizeof(int));
    newGrid[cont] = (int *)calloc(TAMANHO_TOTAL, sizeof(int));
  }


  int lin = 0, col = 0;
  for (i = 0; i < (part); i++) {
    for (j = 0; j < TAMANHO_TOTAL; j++) {
      if (i == 1 && j == 1) {
        lin = i;
        col = j;
        grid[lin][col + 1] = 1;
        grid[lin + 1][col + 2] = 1;
        grid[lin + 2][col] = 1;
        grid[lin + 2][col + 1] = 1;
        grid[lin + 2][col + 2] = 1;
      }
      if (i == 10 && j == 30) {
        lin = i;
        col = j;
        grid[lin][col + 1] = 1;
        grid[lin][col + 2] = 1;
        grid[lin + 1][col] = 1;
        grid[lin + 1][col + 1] = 1;
        grid[lin + 2][col + 1] = 1;
      }
      
    }
  }
  div = i;
  for (origem = 1; origem < numProc; origem++) {
    for (i = div; i < (origem + 1) * part; i++) {
      MPI_Recv(bufRcv, TAMANHO_TOTAL, MPI_INT, origem, tag, MPI_COMM_WORLD,
               MPI_STATUS_IGNORE);
      for (j = 0; j < TAMANHO_TOTAL; j++) {
        grid[i][j] = bufRcv[j];
      }
    }
    div = i;
  }

  printf("Condicao Inicial: %d Celulas Vivas\n", contarTotal());

  for (gen = 0; gen < GERACAO_TOTAL; gen++) {
    for (i = 0; i < TAMANHO_TOTAL; i++) { 
      for (j = 0; j < TAMANHO_TOTAL; j++)
        bufSnd[j] = grid[i][j];
      MPI_Barrier(MPI_COMM_WORLD);
      MPI_Bcast(bufSnd, TAMANHO_TOTAL, MPI_INT, PRIN_PROC, MPI_COMM_WORLD);
      MPI_Barrier(MPI_COMM_WORLD);
    }

    atualizaGrid(numProc, part, rank);

    for (origem = 1; origem < numProc; origem++) {
      for (i = (part * origem); i < (part * (origem + 1)); i++) {
        MPI_Recv(bufRcv, TAMANHO_TOTAL, MPI_INT, origem, tag, MPI_COMM_WORLD,
                 MPI_STATUS_IGNORE);
        for (j = 0; j < TAMANHO_TOTAL; j++) {
          grid[i][j] = bufRcv[j];
        }
      }
    }
  }
  printf("Ultima Geracao: %d Celulas Vivas\n", contarTotal());

  gettimeofday(&end, NULL);
  time = my_difftime(&start, &end);
  printf("Tempo: %d,%d s\n", time->secs, time->usecs);
}

void secProc(int numProc) {
  int tag = 0, dest = 0, i, j, gen, rank, part = TAMANHO_TOTAL / numProc;

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  printf("Sistema secundario iniciado! [%d]\n", rank);

  grid = malloc(sizeof(int *) * TAMANHO_TOTAL);
  newGrid = malloc(sizeof(int *) * TAMANHO_TOTAL);
  for (i = 0; i < TAMANHO_TOTAL; i++) {
    grid[i] = malloc(sizeof(int) * TAMANHO_TOTAL);
    newGrid[i] = malloc(sizeof(int) * TAMANHO_TOTAL);
  }


  for(i=0;i<TAMANHO_TOTAL; i++){
      for(j = 0; j<TAMANHO_TOTAL; j++){
          if(i >= rank*part && i < (rank+1)*part){
              grid[i][j] = 0;
          }else{
              grid[i][j] = 0;
          }

      }
  }

  for (i = rank * part; i < (rank + 1) * part; i++) {
    for (j = 0; j < TAMANHO_TOTAL; j++)
      bufSnd[j] = grid[i][j];
    MPI_Send(bufSnd, TAMANHO_TOTAL, MPI_INT, dest, tag, MPI_COMM_WORLD);
  }

  for (gen = 0; gen < GERACAO_TOTAL; gen++) {
    for (i = 0; i < TAMANHO_TOTAL; i++) { 
      MPI_Barrier(MPI_COMM_WORLD);
      MPI_Bcast(bufRcv, TAMANHO_TOTAL, MPI_INT, PRIN_PROC, MPI_COMM_WORLD);
      MPI_Barrier(MPI_COMM_WORLD);
      for (j = 0; j < TAMANHO_TOTAL; j++)
        grid[i][j] = bufRcv[j];
    }

    atualizaGrid(numProc, part, rank);

    for (i = (part * rank); i < part * (rank + 1); i++) { 
      for (j = 0; j < TAMANHO_TOTAL; j++)
        bufSnd[j] = grid[i][j];
      MPI_Send(bufSnd, TAMANHO_TOTAL, MPI_INT, dest, tag, MPI_COMM_WORLD);
    }
  }
}

int main(int argc, char *argv[]) {
  int rank;   
  int numProc; 
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &numProc);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (rank == 0) { 
    prinProc(numProc);
  } else {
    secProc(numProc);
  }
  MPI_Finalize();
  return 0;
}