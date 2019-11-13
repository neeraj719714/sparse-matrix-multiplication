

#include <stdio.h>
#include<stdlib.h>
#include "mpi.h"
#include <omp.h>


struct sparse_row {
	int* val, * colNum, * numNonZerosIn;
    int n_rows,cnt;
	
};
void sparse_row_initialize(struct sparse_row* sm,int** A, int n, int m)
	{
		int ptr = 0;
        sm->cnt=0;
        sm->n_rows=n;
		sm->numNonZerosIn = (int*)malloc((n+1) * sizeof(int));
		for (int i = 0; i < n; ++i)
			for (int j = 0; j < m; ++j)
				if (A[i][j] != 0)
					sm->cnt++;
		sm->val = (int*)malloc(sm->cnt * sizeof(int));
		sm->colNum = (int*)malloc(sm->cnt * sizeof(int));
		sm->numNonZerosIn[0] = 0;
		for (int i = 0; i < n; ++i)
		{
			int totNonZeros = 0;
			for (int j = 0; j < m; j++)
			{
				if (A[i][j] != 0)
				{
					totNonZeros++;
					sm->val[ptr] = A[i][j];
					sm->colNum[ptr++] = j;
				}
			}
			sm->numNonZerosIn[i + 1] = sm->numNonZerosIn[i] + totNonZeros;
		}
	}


struct sparse_col {
	int* val, * rowNum, * numNonZerosIn;
	int n_cols,cnt;
	
};

void sparse_col_initialize(struct sparse_col* sm,int ** A, int n, int m)
	{
        int ptr=0;
        sm->cnt=0;
        sm->n_cols=m;
		sm->numNonZerosIn = (int*)malloc((m + 1) * sizeof(int));
		for (int i = 0; i < n; ++i)
			for (int j = 0; j < m; ++j)
				if (A[i][j] != 0)
					sm->cnt++;
		sm->val = (int*)malloc(sm->cnt * sizeof(int));
		sm->rowNum = (int*)malloc(sm->cnt * sizeof(int));
		sm->numNonZerosIn[0] = 0;
		for (int i = 0; i < m; ++i)
		{
			int totNonZeros = 0;
			for (int j = 0; j < n; j++)
			{
				if (A[j][i] != 0)
				{
					totNonZeros++;
					sm->val[ptr] = A[i][j];
					sm->rowNum[ptr++] = j;
				}
			}
			sm->numNonZerosIn[i + 1] = sm->numNonZerosIn[i] + totNonZeros;
		}
	}

int **alloc_2d_int(int rows, int cols) {
    int *data = (int *)malloc(rows*cols*sizeof(int));
    int **array= (int **)malloc(rows*sizeof(int*));
    for (int i=0; i<rows; i++)
        array[i] = &(data[cols*i]);

    return array;
}	

int** sparse_mult(struct sparse_row a, struct sparse_col b){
    int *arsum=a.numNonZerosIn;
	int *acid=a.colNum;
	int *aval=a.val;
	int *bcsum=b.numNonZerosIn;
	int *brid=b.rowNum;
	int *bval=b.val; 
	int RA=a.n_rows;
	int CB=b.n_cols;
    
   // for(int i=0;i<=a.n_rows;i++) printf("%d ",arsum[i]); 
    //printf("\n");
   // for(int i=0;i<a.cnt;i++) printf("%d ",acid[i]); 
   // printf("\n");
   // for(int i=0;i<b.cnt;i++) printf("%d ",brid[i]);
    //printf("%d ",a.n_rows);
	int** c= alloc_2d_int(RA,CB);
    
    # pragma parallel for collapse(2) 
	for(int i=1;i<=RA;i++){
	    for(int j=1;j<=CB;j++){
	        int k=arsum[i-1],l=bcsum[j-1],sum=0;
            
	        while(k<arsum[i] && l<bcsum[j]){
        
                if(acid[k]==brid[l]){
	                sum+=aval[k]*bval[l];
                    l++;
                    k++;
	            }
	            else if(acid[k]>brid[l]){
	                l++;
	            }
	            else k++;
            }
	        c[i-1][j-1]=sum;
	    }
	   
	}

	return c;
}



int main(int argc, char *argv[]) {
  int numprocs, rank, namelen;
  char processor_name[MPI_MAX_PROCESSOR_NAME];
  //int iam = 0, np = 1;

  MPI_Init(&argc, &argv);

  MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Get_processor_name(processor_name, &namelen);
/*
  #pragma omp parallel default(shared) private(iam, np)
  {
    np = omp_get_num_threads();
    iam = omp_get_thread_num();
    printf("Hello from thread %d out of %d from process %d out of %d on %s\n",
           iam, np, rank, numprocs, processor_name);
  }
*/

double t1 = MPI_Wtime();
    
int count;
int *val_array,*col_array,*sum_array;
int numrows;
MPI_Status status;

int b_rows=1000,b_cols=1000;
int **b=alloc_2d_int(b_rows,b_cols);

    //time_t t;
    //srand((unsigned) time(&t));

if(rank==0)
{

for(int i=0;i<b_rows;i++)
{
    for(int j=0;j<b_cols;j++)
    {
        if(rand()%100==0)
            b[i][j]=5;
        else b[i][j]=0;
    }
}

}

MPI_Bcast(&b[0][0], b_rows*b_cols , MPI_INT, 0, MPI_COMM_WORLD);

if(rank==0)
{
    int a_rows=1000,a_cols=1000; 

    int **a=alloc_2d_int(a_rows,a_cols);
    int **c=alloc_2d_int(a_rows,b_cols);

    
    
    for(int i=0;i<a_rows;i++)
    {
        for(int j=0;j<a_cols;j++)
        {
            
         if(rand()%100==0)
            a[i][j]=5;
         else a[i][j]=0;
        }
    }
    struct sparse_row sa;
    sparse_row_initialize(&sa,a,a_rows,a_cols);

   

    int send_counts[numprocs];
    int displs[numprocs];
    int n_rows[numprocs];
    int start_row[numprocs];
    int row_each=sa.cnt/numprocs;
    int cur=0,rn=0;
    for(int i=0;i<numprocs;i++)
    {
        if(i==numprocs-1) {
            
            send_counts[i]=sa.cnt-sa.numNonZerosIn[cur];
            displs[i]=sa.numNonZerosIn[cur];
            n_rows[i]=a_rows-cur;
            start_row[i]=cur;
            continue;
            
        }
        
        while(rn<a_rows && (sa.numNonZerosIn[rn]-sa.numNonZerosIn[cur])<row_each)
        {
            rn++;
        }
        send_counts[i]=sa.numNonZerosIn[rn]-sa.numNonZerosIn[cur];
        displs[i]=sa.numNonZerosIn[cur];
        n_rows[i]=rn-cur;
        start_row[i]=cur;
        cur=rn;
    }


    for(int i=1;i<numprocs;i++)
    {
        
        
        MPI_Send( &send_counts[i], 1, MPI_INT,i,0, MPI_COMM_WORLD);

        MPI_Send(sa.val+displs[i] , send_counts[i], MPI_INT,i,0, MPI_COMM_WORLD);

        MPI_Send(sa.colNum+displs[i] , send_counts[i], MPI_INT,i,0, MPI_COMM_WORLD);
        
        MPI_Send( &n_rows[i], 1, MPI_INT,i,0, MPI_COMM_WORLD);
        
        MPI_Send(sa.numNonZerosIn+start_row[i] , n_rows[i]+1, MPI_INT,i,0, MPI_COMM_WORLD);
        
    }
    struct sparse_row sparse_a;
    sparse_a.val=sa.val;
    sparse_a.colNum=sa.colNum;
    sparse_a.numNonZerosIn=sa.numNonZerosIn;
    sparse_a.n_rows=n_rows[0];
    sparse_a.cnt=send_counts[0];
    
    struct sparse_col sparse_b;
    sparse_col_initialize(&sparse_b,b,b_rows,b_cols);
    //printf("%d ",sparse_a.n_rows);
    
    int** result = sparse_mult(sparse_a,sparse_b);
    
    for(int i=0;i<sparse_a.n_rows;i++)
    {
        for(int j=0;j<sparse_b.n_cols;j++)
        {
           // printf("%d ",result[i][j]);
            
            c[i][j]=result[i][j];
            
        }
        //printf("\n");
    }
   for(int i=1;i<numprocs;i++)
        MPI_Recv(&c[start_row[i]][0], n_rows[i]*sparse_b.n_cols , MPI_INT,i,0, MPI_COMM_WORLD,&status);
   
   
   for(int i=0;i<10;i++)
    {
        for(int j=0;j<10;j++)
        {
        printf("%d ",c[i][j]);  
        }
        printf("\n");
    }
}
else
{
    MPI_Recv( &count, 1, MPI_INT,0,0, MPI_COMM_WORLD,&status);
    
    val_array=(int *)malloc(count*sizeof(int));
    col_array=(int *)malloc(count*sizeof(int));

    MPI_Recv(val_array, count, MPI_INT,0,0, MPI_COMM_WORLD,&status);

    MPI_Recv(col_array , count, MPI_INT,0,0, MPI_COMM_WORLD,&status);
    
    MPI_Recv( &numrows, 1, MPI_INT,0,0, MPI_COMM_WORLD,&status);
    
    sum_array=(int *)malloc((numrows+1)*sizeof(int));
    
    MPI_Recv(sum_array, numrows+1, MPI_INT,0,0, MPI_COMM_WORLD,&status);
    
    struct sparse_col sparse_b;
    sparse_col_initialize(&sparse_b,b,b_rows,b_cols);
    
    //printf("process : %d \n",rank);
    
    
    for(int i=1;i<=numrows; i++)
    {
        sum_array[i]-=sum_array[0];
        //printf("%d ",sum_array[i]);
    }
    sum_array[0]=0;
    
    struct sparse_row sparse_a;
    sparse_a.val=val_array;
    sparse_a.colNum=col_array;
    sparse_a.numNonZerosIn=sum_array;
    sparse_a.n_rows=numrows;
    sparse_a.cnt=count;
    //printf("%d ",sparse_a.n_rows);
    
    int** result = sparse_mult(sparse_a,sparse_b);
    
    MPI_Send(&result[0][0] , sparse_a.n_rows*sparse_b.n_cols, MPI_INT,0,0, MPI_COMM_WORLD);
} 
double t2 = MPI_Wtime();
    
    printf( "Elapsed time is %f\n", t2 - t1 );
    
    
MPI_Finalize();
}

  
  
  


