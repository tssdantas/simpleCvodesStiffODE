#include <iostream>
#include <cvodes/cvodes.h>               
#include <nvector/nvector_serial.h>      
#include <sunmatrix/sunmatrix_dense.h>   // Biblioteca correta para matriz densa
#include <sunlinsol/sunlinsol_dense.h>   // Biblioteca correta para solver linear
#include <sundials/sundials_types.h>     
#include <sundials/sundials_math.h>      

// Definição da função de EDOs
int f(realtype t, N_Vector y, N_Vector ydot, void *user_data) {
    realtype y1 = NV_Ith_S(y, 0);
    realtype y2 = NV_Ith_S(y, 1);
    realtype y3 = NV_Ith_S(y, 2);

    NV_Ith_S(ydot, 0) = -0.04 * y1 + 1.0e4 * y2 * y3;
    NV_Ith_S(ydot, 1) = 0.04 * y1 - 1.0e4 * y2 * y3 - 3.0e7 * y2 * y2;
    NV_Ith_S(ydot, 2) = 3.0e7 * y2 * y2;
    
    return 0;
}

int main() {
    realtype t0 = 0.0;  
    realtype t1 = 0.4;
    realtype reltol = 1e-4;
    realtype abstol = 1e-8;

    N_Vector y = N_VNew_Serial(3);
    NV_Ith_S(y, 0) = 1.0;
    NV_Ith_S(y, 1) = 0.0;
    NV_Ith_S(y, 2) = 0.0;

    void *cvode_mem = CVodeCreate(CV_BDF);
    CVodeInit(cvode_mem, f, t0, y);
    CVodeSStolerances(cvode_mem, reltol, abstol);
    CVodeSetUserData(cvode_mem, NULL);

    // Correção: uso das funções corretas para matriz densa e solver linear
    SUNMatrix A = SUNDenseMatrix(3, 3);
    SUNLinearSolver LS = SUNLinSol_Dense(y, A);

    CVodeSetLinearSolver(cvode_mem, LS, A);

    realtype t = t0;
    for (int i = 0; i < 10; i++) {
        CVode(cvode_mem, t1, y, &t, CV_NORMAL);
        std::cout << "t = " << t
                  << ", y1 = " << NV_Ith_S(y, 0)
                  << ", y2 = " << NV_Ith_S(y, 1)
                  << ", y3 = " << NV_Ith_S(y, 2) << std::endl;
    }

    N_VDestroy(y);
    CVodeFree(&cvode_mem);
    SUNLinSolFree(LS);
    SUNMatDestroy(A);

    return 0;
}
