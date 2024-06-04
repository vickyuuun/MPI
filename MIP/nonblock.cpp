//该文件中，包含了普通的MPI优化，优化后的MPI优化，循环划分和块划分的MPI优化
//非阻塞通信的MPI优化
#include<iostream>
#include <stdio.h>
#include<typeinfo>
#include<arm_neon.h>
#include <stdlib.h>
#include<cmath>
#include<mpi.h>
#include <cstring>
using namespace std;
#define N 3000
#define NUM_THREADS 3
float** A = NULL;

struct timespec sts, ets;
time_t dsec;
long dnsec;


struct threadParam_t {    //参数数据结构
    int k;
    int t_id;
};

void A_init() {     //未对齐的数组的初始化
    A = new float* [N];
    for (int i = 0; i < N; i++) {
        A[i] = new float[N];
    }
    for (int i = 0; i < N; i++) {
        A[i][i] = 1.0;
        for (int j = i + 1; j < N; j++) {
            A[i][j] = rand() % 5000;
        }

    }
    for (int k = 0; k < N; k++) {
        for (int i = k + 1; i < N; i++) {
            for (int j = 0; j < N; j++) {
                A[i][j] += A[k][j];
                A[i][j] = (int)A[i][j] % 5000;
            }
        }
    }
}
void A_initAsEmpty() {
    A = new float* [N];
    for (int i = 0; i < N; i++) {
        A[i] = new float[N];
        memset(A[i], 0, N*sizeof(float));
    }

}

void deleteA() {
    for (int i = 0; i < N; i++) {
        delete[] A[i];
    }
    delete A;
}

void print(float** a) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            cout << a[i][j] << " ";
        }
        cout << endl;
    }
}

void LU() {    //普通消元算法
    for (int k = 0; k < N; k++) {
        for (int j = k + 1; j < N; j++) {
            A[k][j] = A[k][j] / A[k][k];
        }
        A[k][k] = 1.0;

        for (int i = k + 1; i < N; i++) {
            for (int j = k + 1; j < N; j++) {
                A[i][j] = A[i][j] - A[i][k] * A[k][j];
            }
            A[i][k] = 0;
        }
    }
}


double LU_mpi_async(int argc, char* argv[]) {  //非阻塞通信
    double start_time = 0;
    double end_time = 0;
    MPI_Init(&argc, &argv);
    int total = 0;
    int rank = 0;
    int i = 0;
    int j = 0;
    int k = 0;
   // MPI_Status status;
    MPI_Comm_size(MPI_COMM_WORLD, &total);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int begin = N / total * rank;
    int end = (rank == total - 1) ? N : N / total * (rank + 1);

    if (rank == 0) {  //0号进程初始化矩阵
        A_init();
        MPI_Request* request = new MPI_Request[N-end];
        for (j = 1; j < total; j++) {
            int b = j * (N / total), e = (j == total - 1) ? N : (j + 1) * (N / total);

            for (i = b; i < e; i++) {
                MPI_Isend(&A[i][0], N, MPI_FLOAT, j, 1, MPI_COMM_WORLD,&request[i-end]);//非阻塞传递矩阵数据
            }

        }
        MPI_Waitall(N-end, request, MPI_STATUS_IGNORE); //等待传递

    }
    else {
        A_initAsEmpty();
        MPI_Request* request = new MPI_Request[end - begin];
        for (i = begin; i < end; i++) {
            MPI_Irecv(&A[i][0], N, MPI_FLOAT, 0, 1, MPI_COMM_WORLD, &request[i-begin]);  //非阻塞接收
        }
        MPI_Waitall(end-begin,request, MPI_STATUS_IGNORE);

    }

    MPI_Barrier(MPI_COMM_WORLD);
    start_time = MPI_Wtime();
    for (k = 0; k < N; k++) {
        if ((begin <= k && k < end)) {
            for (j = k + 1; j < N; j++) {
                A[k][j] = A[k][j] / A[k][k];
            }
            A[k][k] = 1.0;
            MPI_Request* request = new MPI_Request[total - 1-rank];  //非阻塞传递
            for (j = rank + 1; j < total; j++) { //块划分中，已经消元好且进行了除法置1的行向量仅

                MPI_Isend(&A[k][0], N, MPI_FLOAT, j, 0, MPI_COMM_WORLD,&request[j-rank-1]);//0号消息表示除法完毕
            }
            MPI_Waitall(total-1-rank, request, MPI_STATUS_IGNORE);
            if (k == end - 1)
                break; //若执行完自身的任务，可直接跳出
        }
        else {
            int src = k / (N / total);
            MPI_Request request;
            MPI_Irecv(&A[k][0], N, MPI_FLOAT, src, 0, MPI_COMM_WORLD, &request);
            MPI_Wait(&request, MPI_STATUS_IGNORE);         //实际上仍然是阻塞接收，因为接下来的操作需要这些数据
        }
        for (i = max(begin, k + 1); i < end; i++) {
            for (j = k + 1; j < N; j++) {
                A[i][j] = A[i][j] - A[i][k] * A[k][j];
            }
            A[i][k] = 0;
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);	//各进程同步
    if (rank == total - 1) {
        end_time = MPI_Wtime();
        printf("平凡MPI，块划分+非阻塞耗时：%.4lf ms\n", 1000 * (end_time - start_time));
    //    print(A);
    }
    MPI_Finalize();
    return end_time - start_time;
}

void cal(void(*func)()) {
    A_init();
    timespec_get(&sts, TIME_UTC);
    func();
    timespec_get(&ets, TIME_UTC);
    dsec = ets.tv_sec - sts.tv_sec;
    dnsec = ets.tv_nsec - sts.tv_nsec;
    if (dnsec < 0) {
        dsec--;
        dnsec += 1000000000ll;
    }
}

int main(int argc, char* argv[]) {

   LU_mpi_async(argc, argv);

  

}
