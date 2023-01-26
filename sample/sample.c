#include <stdio.h>
#include <lbfgs.h>

static lbfgsfloatval_t evaluate(
    void *instance,
    const lbfgsfloatval_t *x,
    lbfgsfloatval_t *g,
    const int n,
    const lbfgsfloatval_t step
    )
{
    int i;
    lbfgsfloatval_t fx = 0.0;
    lbfgsfloatval_t d = 0.0;
    lbfgsfloatval_t t1 = 0.0;
    lbfgsfloatval_t t2 = 0.0;
    lbfgsfloatval_t t3 = 0.0;
        for (i = 0; i < (n / 2) - 2; i += 2) {
            g[i] = (12.0 / 5.0) * x[i] + (1.0 / 5.0) * x[i + 1] - (12.0 / 5.0) * x[i + 2] + (1.0 / 5.0) * x[i + 3];
            g[i + 1] = (1.0 / 5.0) * x[i] + (4.0 / 15.0) * x[i + 1] - (1.0 / 5.0) * x[i + 2] - (1.0 / 15.0) * x[i + 3];
            g[i + 2] = -(12.0 / 5.0) * x[i] - (1.0 / 5.0) * x[i + 1] + (12.0 / 5.0) * x[i + 2] - (1.0 / 5.0) * x[i + 3];
            g[i + 3] = (1.0 / 5.0) * x[i] - (1.0 / 15.0) * x[i + 1] - (1.0 / 5.0) * x[i + 2] + (4.0 / 15.0) * x[i + 3];



            t1 += (6.0 / 5.0) * x[i] * x[i] + (2.0 / 15.0) * x[i + 1] * x[i + 1] + (6.0 / 5.0) * x[i + 2] * x[i + 2] + (2.0 / 15.0) * x[i + 3] * x[i + 3]
                + (1.0 / 5.0) * x[i] * x[i + 1] - (12.0 / 5.0) * x[i] * x[i + 2] + (1.0 / 5.0) * x[i] * x[i + 3]
                - (1.0 / 5.0) * x[i + 2] * x[i + 1] - (1.0 / 15.0) * x[i + 1] * x[i + 3] - (1.0 / 5.0) * x[i + 2] * x[i + 3];


        }
        printf("aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa\n");
        for (int jjj = 0; jjj < n / 2; jjj = jjj + 2) {
            printf("x[%d]  %f x[%d]  %f g[%d]  %f g[%d]  %f\n",jjj, x[jjj],jjj + 1,x[jjj + 1],jjj,g[jjj],jjj + 1,g[jjj + 1]);
        }
        printf("                       %f\n", t1);
        for (i = n / 2; i < n - 2; i += 2) {
            g[i] = (12.0 / 5.0) * x[i] + (1.0 / 5.0) * x[i + 1] - (12.0 / 5.0) * x[i + 2] + (1.0 / 5.0) * x[i + 3];
            g[i + 1] = (1.0 / 5.0) * x[i] + (4.0 / 15.0) * x[i + 1] - (1.0 / 5.0) * x[i + 2] - (1.0 / 15.0) * x[i + 3];
            g[i + 2] = -(12.0 / 5.0) * x[i] - (1.0 / 5.0) * x[i + 1] + (12.0 / 5.0) * x[i + 2] - (1.0 / 5.0) * x[i + 3];
            g[i + 3] = (1.0 / 5.0)* x[i] - (1.0 / 15.0) * x[i + 1] - (1.0 / 5.0) * x[i + 2] + (4.0 / 15.0) * x[i + 3];



            t2 += (6.0 / 5.0) * x[i] * x[i] + (2.0 / 15.0) * x[i + 1] * x[i + 1] + (6.0 / 5.0) * x[i + 2] * x[i + 2] + (2.0 / 15.0) * x[i + 3] * x[i + 3]
                + (1.0 / 5.0) * x[i] * x[i + 1] - (12.0 / 5.0) * x[i] * x[i + 2] + (1.0 / 5.0) * x[i] * x[i + 3]
                - (1.0 / 5.0) * x[i + 2] * x[i + 1] - (1.0 / 15.0) * x[i + 1] * x[i + 3] - (1.0 / 5.0) * x[i + 2] * x[i + 3];


        }
        printf("bbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb\n");
        for (int jjj = n / 2; jjj < n; jjj = jjj + 2) {
            printf("x[%d]  %f x[%d]  %f g[%d]  %f g[%d]  %f\n", jjj, x[jjj], jjj + 1, x[jjj + 1], jjj, g[jjj], jjj + 1, g[jjj + 1]);
        }
        printf("                       %f\n", t2);
        d = ((x[8] - x[18]) / 2.2 )* ((x[8] - x[18]) / 2.2);
        if (d < 1) {
            t3 = 2.0 * (1 / (d * d) + (d * d) - 2.0);
        }
        else {
            t3 = 0.0;
        }
        printf("ccccccccccccccccccccccccccccccccccc\n");
        printf("%f   %f   %f\n", x[18], x[8], t3);

        fx += t1 + t2 + t3;
        printf("ddddddddddddddddddddddddddddddddddd\n");
        printf("                       %f\n", fx);
    
    return fx;
}

static int progress(
    void *instance,
    const lbfgsfloatval_t *x,
    const lbfgsfloatval_t *g,
    const lbfgsfloatval_t fx,
    const lbfgsfloatval_t xnorm,
    const lbfgsfloatval_t gnorm,
    const lbfgsfloatval_t step,
    int n,
    int k,
    int ls
    )
{
    printf("Iteration %d:\n", k);
    printf("  fx = %f, x[0] = %f, x[1] = %f, x[2] = %f, x[3] = %f\n", fx, x[0], x[1],x[2],x[3]);
    printf("  fx = %f, g[0] = %f, g[1] = %f, g[2] = %f, g[3] = %f\n", fx, g[0], g[1],g[2],g[3]);
    printf("  xnorm = %f, gnorm = %f, step = %f\n", xnorm, gnorm, step);
    printf("\n");
    return 0;
   

}

#define N   28

int main(int argc, char *argv[])
{
    int i, ret = 0;
    lbfgsfloatval_t fx;
    lbfgsfloatval_t *x = lbfgs_malloc(N);
    lbfgs_parameter_t param;
    lbfgsfloatval_t *kk = lbfgs_malloc(N);;
    int kj = 0;

    if (x == NULL) {
        printf("ERROR: Failed to allocate a memory block for variables.\n");
        return 1;
    }
    //printf("\n");

    /* Initialize the variables. 
       変数の初期化
    */

    lbfgsfloatval_t takasaI1[7] = { 0.0,0.0,1.0,0.0,-1.0,0.0,0.0 };
   
    lbfgsfloatval_t takasaI2[7] = { 0.0,0.0,-1.0,0.0,1.0,0.0,0.0 };

    lbfgsfloatval_t takasaJ1[7] = { 0.0,0.0,-1.0,0.0,1.0,0.0,0.0 };

    lbfgsfloatval_t takasaJ2[7] = { 0.0,0.0,1.0,0.0,-1.0,0.0,0.0 };

    lbfgsfloatval_t katamukiI1[7] = { 0.0,1.0,-1.0,-1.0,1.0,0.0,1.0 };

    lbfgsfloatval_t katamukiI2[7] = {0.0,-1.0,1.0,1.0,-1.0,0.0,-1.0};

    lbfgsfloatval_t katamukiJ1[7] = { 0.0,-1.0,1.0,1.0,-1.0,0.0,-1.0 };

    lbfgsfloatval_t katamukiJ2[7] = { 0.0,1.0,-1.0,-1.0,1.0,0.0,1.0 };

    lbfgsfloatval_t takasa1I2[7] = { 0.0,0.0,-1.0,0.0,0.0,0.0,0.0 };

    lbfgsfloatval_t takasa1J1[7] = { 0.0,0.0,0.0,0.0,1.0,0.0,0.0 };

    lbfgsfloatval_t katamuki1I2[7] = { 0.0,-1.0,1.0,0.0,0.0,0.0,0.0 };

    lbfgsfloatval_t katamuki1J1[7] = { 0.0,0.0,0.0,1.0,-1.0,0.0,0.0 };
     
    int j = 0;
    for (i = 0;i < N / 2;i += 2) {
        x[i] = takasaJ1[j];
        x[i + 1] = katamukiJ1[j];

        j = j + 1;
    }
    int k = 0;
    for (i = N / 2; i < N; i += 2) {
        x[i] = takasaI2[k];
        x[i + 1] = katamukiI2[k];
 
        k = k + 1;
    }
   
    /* Initialize the parameters for the L-BFGS optimization.
       L-BFGS最適化のためのパラメータの初期化 
    */
    lbfgs_parameter_init(&param);
    /*param.linesearch = LBFGS_LINESEARCH_BACKTRACKING;*/

    /*
        Start the L-BFGS optimization; this will invoke the callback functions
        evaluate() and progress() when necessary.
        L-BFGS最適化を開始します．必要に応じてコールバック関数
        evaluate()とprogress()を呼び出す
     */
    
    ret = lbfgs(N, x, &fx, evaluate, progress, NULL, &param);

    
    /* Report the result.         &paramにはデフォルトのパラメータがはいっている
       結果を報告する
    */
    printf("L-BFGS optimization terminated with status code = %d\n", ret);
    printf("  fx = %f, x[0] = %f, x[1] = %f, x[2] = %f, x[3] = %f\n", fx, x[0], x[1],x[2],x[3]);

    lbfgs_free(x);
    return 0;
}
