// @useDynLib RAINBOWR
// @importFrom Rcpp evalCpp
//
// [[Rcpp::depends(Rcpp)]]
// [[Rcpp::depends(RcppEigen)]]
#include <Rcpp.h>
#include <RcppEigen.h>

using namespace Eigen;

typedef Map<MatrixXd> MapMat;
typedef Map<VectorXd> MapVec;

MatrixXd solve(MatrixXd A, MatrixXd B){

  Rcpp::Function f("solve");

  Rcpp::NumericMatrix sol = f(Rcpp::wrap(A), Rcpp::wrap(B));
  MapMat Sol = Rcpp::as<MapMat>(sol);

  return Sol;
}

MatrixXd  inv(MatrixXd A){

  Rcpp::Environment pkg = Rcpp::Environment::namespace_env("MASS");
  Rcpp::Function f = pkg["ginv"];

  Rcpp::NumericMatrix sol = f(Rcpp::wrap(A));
  MapMat Sol = Rcpp::as<MapMat>(sol);

  return Sol;
}


Rcpp::List  svd(MatrixXd A, int nu, int nv){

  Rcpp::Function f("svd");

  Rcpp::List sol = f(Rcpp::wrap(A), nu, nv);

  return sol;
}

Rcpp::NumericVector matching(VectorXd A, VectorXd B){
  const int n1(A.rows()), n2(B.rows());

  VectorXd xs = -VectorXd::Ones(n1, 1);

  for(int i = 0; i < n1; i++){
    double trial = A[i];

    for(int j = 0; j < n2; j++){
      if(B[j] == trial){
        xs[i] = j;
        break;
      }
    }
  }

  Rcpp::NumericVector xs2 = Rcpp::wrap(xs);

  for(int i = 0; i < n1; i ++){
    if(xs2[i] == -1){
      xs2[i] = NA_REAL;
    }
  }
  return xs2;
}

MatrixXd extract(MatrixXd A, VectorXd B, int method = 0){
  const int n1(A.rows()), m1(B.rows()), n2(B.rows());
  MatrixXd A2 = MatrixXd::Zero(n2, m1);


  if(method == 0){
    for(int i = 0; i < n2; i++){
      A2.row(i) = A.row(B[i]);
    }
  }


  if(method == 1){
    MatrixXd A2 = MatrixXd::Zero(n1, n2);
    for(int i = 0; i < n2; i++){
      A2.col(i) = A.col(B[i]);
    }
  }
  return A2;
}

MatrixXd power(MatrixXd A, int num = 2){
  const MatrixXd Aseg = A;
  MatrixXd Ares = A;
  for(int i = 0; i < (num - 1); i++){
    Ares = (Ares.array() * Aseg.array()).matrix();
  }

  return Ares;
}

MatrixXd elepro(MatrixXd A, MatrixXd B, int method = 0){
  const int n1(A.rows()), n2(B.rows()), m1(A.cols()), m2(B.cols());

  MatrixXd Bs = MatrixXd::Zero(n1, m1);
  if(method == 0){
    MatrixXd Bs = B;
  }

  if(method == 1){
    if((m2 != 1) || (n1 != n2)){
      Rcpp::stop("The dimension of input matrix or method type is uncorrect!");
    }

    for(int i = 0; i < m1; i++){
      Bs.col(i) = B;
    }
  }

  if(method == 2){
    if((n2 != 1) || (m1 != m2)){
      Rcpp::stop("The dimension of input matrix or method type is uncorrect!");
    }

    for(int i = 0; i < n1; i++){
      Bs.row(i) = B;
    }
  }

  MatrixXd AB = (A.array() * Bs.array()).matrix();

  return AB;
}

MatrixXd elediv(MatrixXd A, MatrixXd B, int method = 0){
  const int n1(A.rows()), n2(B.rows()), m1(A.cols()), m2(B.cols());

  MatrixXd Bs = MatrixXd::Zero(n1, m1);
  if(method == 0){
    Bs = B;
  }

  if(method == 1){
    if((m2 != 1) || (n1 != n2)){
      Rcpp::stop("The dimension of input matrix or method type is uncorrect!");
    }

    for(int i = 0; i < m1; i++){
      Bs.col(i) = B;
    }
  }

  if(method == 2){
    if((n2 != 1) || (m1 != m2)){
      Rcpp::stop("The dimension of input matrix or method type is uncorrect!");
    }

    for(int i = 0; i < n1; i++){
      Bs.row(i) = B;
    }
  }

  MatrixXd AdivB = (A.array() / Bs.array()).matrix();

  return AdivB;
}

MatrixXd diagelements(MatrixXd A){
  const int n(A.rows()), m(A.cols());
  int nmmin = 0;

  if(n >= m){
    nmmin = m;
  }else{
    nmmin = n;
  }
  MatrixXd Adiag = MatrixXd::Zero(nmmin, 1);
  for(int i = 0; i < nmmin; i++){
    Adiag(i, 0) = A(i, i);
  }
  return Adiag;
}

MatrixXd cbind(MatrixXd A, MatrixXd B){
  const int n1(A.rows()), n2(B.rows()), m1(A.cols()), m2(B.cols());

  if(n1 != n2){
    Rcpp::stop("The numbers of rows don't match between two matrices!");
  }

  MatrixXd ABcbind = MatrixXd::Zero(n1, m1 + m2);
  ABcbind.block(0, 0, n1, m1) = A;
  ABcbind.block(0, m1, n1, m2) = B;

  return(ABcbind);
}

MatrixXd rbind(MatrixXd A, MatrixXd B){
  const int n1(A.rows()), n2(B.rows()), m1(A.cols()), m2(B.cols());

  if(m1 != m2){
    Rcpp::stop("The numbers of rows don't match between two matrices!");
  }

  MatrixXd ABrbind = MatrixXd::Zero(n1 + n2, m1);
  ABrbind.block(0, 0, n1, m1) = A;
  ABrbind.block(n1, 0, n2, m1) = B;

  return(ABrbind);
}


MatrixXd reshape(MatrixXd A, int nrow, int ncol, int method = 0){
  const int n(A.rows()), m(A.cols());
  int n0 = 0, m0 = 0;

  MatrixXd A2 = MatrixXd::Zero(nrow, ncol);
  int nm = n * m;
  int nm2 = nrow * ncol;

  if(nm != nm2){
    Rcpp::stop("Size of the reshaped matrix is not same as that of the original Matrix!");
  }

  if(method == 0){
    n0 = 1;
    m0 = n * m;
  }else{
    n0 = n * m;
    m0 = 1;
  }

  MatrixXd A_vec = MatrixXd::Zero(n0, m0);
  if(method == 0){
    for(int i = 0; i < n; i++){
      A_vec.block(0, m * i, 1, m) = A.row(i);
    }
  }else{
    for(int i = 0; i < m; i++){
      A_vec.block(n * i, 0, n, 1) = A.col(i);
    }
  }

  // reshape vector -> Matrix

  if(method == 0){
    for(int i = 0; i < nrow; i++){
      A2.row(i) = A_vec.block(0, ncol * i, 1, ncol);
    }
  }else{
    for(int i = 0; i < ncol; i++){
      A2.col(i) = A_vec.block(nrow * i, 0, nrow, 1);
    }
  }

  return(A2);
}

MatrixXd crossprod(MatrixXd A, MatrixXd B){
  const int n1(A.rows()), n2(B.rows());

  if(n1 != n2){
    Rcpp::stop("Crossproduct cannnot be calculated! Check the dimension of two matrices!");
  }

  MatrixXd At = A.transpose();
  MatrixXd AtB = At * B;

  return AtB;
}



MatrixXd tcrossprod(MatrixXd A, MatrixXd B){
  const int m1(A.cols()), m2(B.cols());

  if(m1 != m2){
    Rcpp::stop("Tcrossproduct cannnot be calculated! Check the dimension of two matrices!");
  }

  MatrixXd Bt = B.transpose();
  MatrixXd ABt = A * Bt;

  return ABt;
}




Rcpp::List eigen(MatrixXd A, bool symmetric = true, bool ainv = false){
  Rcpp::Function f("eigen");

  Rcpp::List sol = f(Rcpp::wrap(A), symmetric);
  MapMat D = Rcpp::as<MapMat>(sol[0]);
  MapMat U = Rcpp::as<MapMat>(sol[1]);

  const int m(U.cols());
  if(!ainv){
    return Rcpp::List::create(Rcpp::Named("values") = Rcpp::wrap(D),
                              Rcpp::Named("vectors") = Rcpp::wrap(U));
  }else{
    MatrixXd D2 = MatrixXd::Zero(1, m);
    D2.row(0) = D.col(0);
    MatrixXd UdivD2 = elediv(U, D2, 2);

    MatrixXd Ainv = tcrossprod(UdivD2, U);

    return Rcpp::List::create(Rcpp::Named("values") = Rcpp::wrap(D),
                              Rcpp::Named("vectors") = Rcpp::wrap(U),
                              Rcpp::Named("inverse") = Rcpp::wrap(Ainv));
  }
}


Rcpp::List eigen2(MatrixXd A, bool ainv = false){
  SelfAdjointEigenSolver<MatrixXd> es(A);


  MatrixXd U = es.eigenvectors();
  MatrixXd D = es.eigenvalues();
  const int m(U.cols());
  MatrixXd D2 = MatrixXd::Zero(m, 1);

  for(int i = 0; i < m; i++){
    if(i < m / 2){
      U.col(i).swap(U.col(m - i - 1));
    }

    D2.row(i) = D.row(m - i - 1);
  }

  if(!ainv){
    return Rcpp::List::create(Rcpp::Named("values") = Rcpp::wrap(D2),
                              Rcpp::Named("vectors") = Rcpp::wrap(U));
  }else{
    MatrixXd D4 = MatrixXd::Zero(1, m);
    D4.row(0) = D2.col(0);
    MatrixXd UdivD3 = elediv(U, D4, 2);

    MatrixXd Ainv = tcrossprod(UdivD3, U);

    return Rcpp::List::create(Rcpp::Named("values") = Rcpp::wrap(D2),
                              Rcpp::Named("vectors") = Rcpp::wrap(U),
                              Rcpp::Named("inverse") = Rcpp::wrap(Ainv));
  }
}


// aHinv : the function to calculate aHinvb
MatrixXd aHinvb(MatrixXd A, MatrixXd B,
                MatrixXd U, MatrixXd EV)
{
  MatrixXd va = U.adjoint() * A;
  MatrixXd vb = U.adjoint() * B;

  MatrixXd EV2 = power(EV, 2);
  MatrixXd EV3 = power(EV, 3);


  MatrixXd aH1b = va.adjoint() * elediv(vb, EV, 1);
  MatrixXd aH2b = va.adjoint() * elediv(vb, EV2, 1);
  MatrixXd aH3b = va.adjoint() * elediv(vb, EV3, 1);

  MatrixXd aHbs = rbind(rbind(aH1b, aH2b), aH3b);

  return aHbs;
}

// aPb_series : the function to calculate aPb, aPPb, aPPPPb, tr(P), tr(PP) (and P)
Rcpp::List aPb_series(MatrixXd A, MatrixXd B,
                      MatrixXd U, MatrixXd EV,
                      MatrixXd X)
{
  const int p(X.cols()), m(U.cols());
  const MatrixXd ones = MatrixXd::Ones(m, 1);

  const MatrixXd aPb_0 = aHinvb(A, B, U, EV);
  const MatrixXd aPX_0 = aHinvb(A, X, U, EV);
  const MatrixXd bPX_0 = aHinvb(B, X, U, EV);
  const MatrixXd XPX_0 = aHinvb(X, X, U, EV);

  const MatrixXd XPXs_0 = reshape(XPX_0.block(0, 0, p, p), 1, p * p, 1);
  const MatrixXd XPPXs_0 = reshape(XPX_0.block(p, 0, p, p), 1, p * p, 1);
  const MatrixXd XPPPXs_0 = reshape(XPX_0.block(2 * p, 0, p, p), 1, p * p, 1);


  VectorXd aPbs = VectorXd::Zero(p + 1, 1);
  VectorXd aPPbs = VectorXd::Zero(p + 1, 1);
  VectorXd aPPPbs = VectorXd::Zero(p + 1, 1);

  VectorXd tr_Ps = VectorXd::Zero(p + 1, 1);
  VectorXd tr_PPs = VectorXd::Zero(p + 1, 1);

  MatrixXd aPXs = MatrixXd::Zero(p + 1, p);
  MatrixXd aPPXs = MatrixXd::Zero(p + 1, p);
  MatrixXd aPPPXs = MatrixXd::Zero(p + 1, p);
  MatrixXd bPXs = MatrixXd::Zero(p + 1, p);
  MatrixXd bPPXs = MatrixXd::Zero(p + 1, p);
  MatrixXd bPPPXs = MatrixXd::Zero(p + 1, p);

  MatrixXd XPXs = MatrixXd::Zero(p + 1, p * p);
  MatrixXd XPPXs = MatrixXd::Zero(p + 1, p * p);
  MatrixXd XPPPXs = MatrixXd::Zero(p + 1, p * p);

  aPbs[0] = aPb_0(0, 0);
  aPXs.row(0) = aPX_0.row(0);
  bPXs.row(0) = bPX_0.row(0);
  XPXs.row(0) = XPXs_0.row(0);

  aPPbs[0] = aPb_0(1, 0);
  aPPXs.row(0) = aPX_0.row(1);
  bPPXs.row(0) = bPX_0.row(1);
  XPPXs.row(0) = XPPXs_0.row(0);

  aPPPbs[0] = aPb_0(2, 0);
  aPPPXs.row(0) = aPX_0.row(2);
  bPPPXs.row(0) = bPX_0.row(2);
  XPPPXs.row(0) = XPPPXs_0.row(0);

  tr_Ps[0] = elediv(ones, EV).sum();
  tr_PPs[0] = elediv(ones, power(EV, 2)).sum();

  for(int i = 0; i < p; i++){
    // Calculate aPb
    double aPb_next = aPbs[i] - (aPXs(i, i) * bPXs(i, i)) / XPXs(i, i * (p + 1));
    aPbs[i + 1] = aPb_next;

    VectorXd aPX_next = aPXs.row(i) - (aPXs(i, i) * XPXs.block(i, i * p, 1, p)) / XPXs(i, i * (p + 1));
    aPXs.row(i + 1) = aPX_next;

    VectorXd bPX_next = bPXs.row(i) - (bPXs(i, i) * XPXs.block(i, i * p, 1, p)) / XPXs(i, i * (p + 1));
    bPXs.row(i + 1) = bPX_next;

    MatrixXd XPX_next = XPXs.row(i) - reshape(crossprod(XPXs.block(i, i * p, 1, p),
                                              XPXs.block(i, i * p, 1, p)), 1, p * p, 1) / XPXs(i, i * (p + 1));
    XPXs.row(i + 1) = XPX_next;


    // Calculate aPPbs
    double aPPb_next = aPPbs[i] + (aPXs(i, i) * bPXs(i, i) * XPPXs(i, i * (p + 1))) / pow(XPXs(i, i * (p + 1)), 2) -
      (aPXs(i, i) * bPPXs(i, i)) / XPXs(i, i * (p + 1)) - (aPPXs(i, i) * bPXs(i, i)) / XPXs(i, i * (p + 1));
    aPPbs[i + 1] = aPPb_next;

    VectorXd aPPXs_next = aPPXs.row(i) + (aPXs(i, i) * XPXs.block(i, i * p, 1, p) * XPPXs(i, i * (p + 1))) / pow(XPXs(i, i * (p + 1)), 2) -
      (aPXs(i, i) * XPPXs.block(i, i * p, 1, p)) / XPXs(i, i * (p + 1)) - (aPPXs(i, i) * XPXs.block(i, i * p, 1, p)) / XPXs(i, i * (p + 1));
    aPPXs.row(i + 1) = aPPXs_next;

    VectorXd bPPXs_next = bPPXs.row(i) + (bPXs(i, i) * XPXs.block(i, i * p, 1, p) * XPPXs(i, i * (p + 1))) / pow(XPXs(i, i * (p + 1)), 2) -
      (bPXs(i, i) * XPPXs.block(i, i * p, 1, p)) / XPXs(i, i * (p + 1)) - (bPPXs(i, i) * XPXs.block(i, i * p, 1, p)) / XPXs(i, i * (p + 1));
    bPPXs.row(i + 1) = bPPXs_next;

    MatrixXd XPPXs_next = XPPXs.row(i) +
      reshape(crossprod(XPXs.block(i, i * p, 1, p), XPXs.block(i, i * p, 1, p)), 1, p * p, 1) * XPPXs(i, i * (p + 1)) / pow(XPXs(i, i * (p + 1)), 2) -
      reshape(crossprod(XPXs.block(i, i * p, 1, p), XPPXs.block(i, i * p, 1, p)), 1, p * p, 1) / XPXs(i, i * (p + 1)) -
      reshape(crossprod(XPPXs.block(i, i * p, 1, p), XPXs.block(i, i * p, 1, p)), 1, p * p, 1) / XPXs(i, i * (p + 1));
    XPPXs.row(i + 1) = XPPXs_next;



    // Calculate aPPPb
    double aPPPb_next = aPPPbs[i] - (aPXs(i, i) * bPXs(i, i) * pow(XPPXs(i, i * (p + 1)), 2)) / pow(XPXs(i, i * (p + 1)), 3) -
      (aPXs(i, i) * bPPPXs(i, i)) / XPXs(i, i * (p + 1)) - (aPPPXs(i, i) * bPXs(i, i)) / XPXs(i, i * (p + 1)) -
      (aPPXs(i, i) * bPPXs(i, i)) / XPXs(i, i * (p + 1)) + (aPXs(i, i) * bPPXs(i, i) * XPPXs(i, i * (p + 1))) / pow(XPXs(i, i * (p + 1)), 2) +
      (aPPXs(i, i) * bPXs(i, i) * XPPXs(i, i * (p + 1))) / pow(XPXs(i, i * (p + 1)), 2) +
      (aPXs(i, i) * bPXs(i, i) * XPPPXs(i, i * (p + 1))) / pow(XPXs(i, i * (p + 1)), 2);
    aPPPbs[i + 1] = aPPPb_next;

    VectorXd aPPPX_next = aPPPXs.row(i) - (aPXs(i, i) * XPXs.block(i, i * p, 1, p) * pow(XPPXs(i, i * (p + 1)), 2)) / pow(XPXs(i, i * (p + 1)), 3) -
      (aPXs(i, i) * XPPPXs.block(i, i * p, 1, p)) / XPXs(i, i * (p + 1)) - (aPPPXs(i, i) * XPXs.block(i, i * p, 1, p)) / XPXs(i, i * (p + 1)) -
      (aPPXs(i, i) * XPPXs.block(i, i * p, 1, p)) / XPXs(i, i * (p + 1)) + (aPXs(i, i) * XPPXs.block(i, i * p, 1, p) * XPPXs(i, i * (p + 1))) / pow(XPXs(i, i * (p + 1)), 2) +
      (aPPXs(i, i) * XPXs.block(i, i * p, 1, p) * XPPXs(i, i * (p + 1))) / pow(XPXs(i, i * (p + 1)), 2) +
      (aPXs(i, i) * XPXs.block(i, i * p, 1, p) * XPPPXs(i, i * (p + 1))) / pow(XPXs(i, i * (p + 1)), 2);
    aPPPXs.row(i + 1) = aPPPX_next;

    VectorXd bPPPX_next = bPPPXs.row(i) - (bPXs(i, i) * XPXs.block(i, i * p, 1, p) * pow(XPPXs(i, i * (p + 1)), 2)) / pow(XPXs(i, i * (p + 1)), 3) -
      (bPXs(i, i) * XPPPXs.block(i, i * p, 1, p)) / XPXs(i, i * (p + 1)) - (bPPPXs(i, i) * XPXs.block(i, i * p, 1, p)) / XPXs(i, i * (p + 1)) -
      (bPPXs(i, i) * XPPXs.block(i, i * p, 1, p)) / XPXs(i, i * (p + 1)) + (bPXs(i, i) * XPPXs.block(i, i * p, 1, p) * XPPXs(i, i * (p + 1))) / pow(XPXs(i, i * (p + 1)), 2) +
      (bPPXs(i, i) * XPXs.block(i, i * p, 1, p) * XPPXs(i, i * (p + 1))) / pow(XPXs(i, i * (p + 1)), 2) +
      (bPXs(i, i) * XPXs.block(i, i * p, 1, p) * XPPPXs(i, i * (p + 1))) / pow(XPXs(i, i * (p + 1)), 2);
    bPPPXs.row(i + 1) = bPPPX_next;

    MatrixXd XPPPX_next = XPPPXs.row(i) -
      (reshape(crossprod(XPXs.block(i, i * p, 1, p), XPXs.block(i, i * p, 1, p)), 1, p * p, 1) * pow(XPPXs(i, i * (p + 1)), 2)) / pow(XPXs(i, i * (p + 1)), 3) -
      reshape(crossprod(XPXs.block(i, i * p, 1, p), XPPPXs.block(i, i * p, 1, p)), 1, p * p, 1) / XPXs(i, i * (p + 1)) -
      reshape(crossprod(XPPPXs.block(i, i * p, 1, p), XPXs.block(i, i * p, 1, p)), 1, p * p, 1)  / XPXs(i, i * (p + 1)) -
      reshape(crossprod(XPPXs.block(i, i * p, 1, p), XPPXs.block(i, i * p, 1, p)), 1, p * p, 1)  / XPXs(i, i * (p + 1)) +
      (reshape(crossprod(XPXs.block(i, i * p, 1, p), XPPXs.block(i, i * p, 1, p)), 1, p * p, 1) * XPPXs(i, i * (p + 1))) / pow(XPXs(i, i * (p + 1)), 2) +
      (reshape(crossprod(XPPXs.block(i, i * p, 1, p), XPXs.block(i, i * p, 1, p)), 1, p * p, 1)  * XPPXs(i, i * (p + 1))) / pow(XPXs(i, i * (p + 1)), 2) +
      (reshape(crossprod(XPXs.block(i, i * p, 1, p), XPXs.block(i, i * p, 1, p)), 1, p * p, 1)  * XPPPXs(i, i * (p + 1))) / pow(XPXs(i, i * (p + 1)), 2);
    XPPPXs.row(i + 1) = XPPPX_next;



    // Calculate trace(P)
    double tr_P_next = tr_Ps[i] - XPPXs(i, i * (p + 1)) / XPXs(i, i * (p + 1));
    tr_Ps[i + 1] = tr_P_next;

    // Calculate trace(PP)
    double tr_PP_next = tr_PPs[i] + pow(XPPXs(i, i * (p + 1)), 2) / pow(XPXs(i, i * (p + 1)), 2) -
      2 * XPPPXs(i, i * (p + 1)) / XPXs(i, i * (p + 1));
    tr_PPs[i + 1] = tr_PP_next;
  }

  return Rcpp::List::create(Rcpp::Named("aPb") = Rcpp::wrap(aPbs),
                            Rcpp::Named("aPPb") = Rcpp::wrap(aPPbs),
                            Rcpp::Named("aPPPb") = Rcpp::wrap(aPPPbs),
                            Rcpp::Named("tr.P") = Rcpp::wrap(tr_Ps),
                            Rcpp::Named("tr.PP") = Rcpp::wrap(tr_PPs),
                            Rcpp::Named("XtHinvX") = Rcpp::wrap(XPX_0.block(0, 0, p, p)));
}

// llik : the function to calculate log-likelihood
double llik(int n, int p, double logH, double yPy, double logXtX, double logXtHinvX, bool REML = true){
  double pi = 3.14159;
  double l = 0;

  if(REML){
    l = ((n - p) * std::log((n - p) / (2 * pi)) - (n - p) + logXtX -
      logH - logXtHinvX -  (n - p) * std::log(yPy)) / 2;
    Rcpp::Named("l.REML") = l;
  }else{
    l = (n * std::log(n / (2 * pi)) - n - logH - n * std::log(yPy)) / 2;
    Rcpp::Named("l.ML") = l;
  }
  return(l);
}

// score_func : the function to calculate the score function
double score_func(int n, int p, double tr_HinvG, double yPy,
                  double yPGPy, double tr_PG, bool REML = true){
  double l1 = 0;
  if(REML){
    l1 = (- tr_PG + (n - p) * yPGPy / yPy) / 2;
    Rcpp::Named("l1.REML") = l1;
  }else{
    l1 = (- tr_HinvG + n * yPGPy / yPy) / 2;
    Rcpp::Named("l1.ML") = l1;
  }
  return(l1);
}

// hess_func : the function to calculate the score function
double hess_func(int n, int p, double tr_HinvG2, double yPy, double yPGPy,
                 double yPGPGPy, double tr_PG2, bool REML = true){
  double l2 = 0;
  if(REML){
    l2 = (tr_PG2 - (n - p) * (2 * yPGPGPy * yPy - pow(yPGPy, 2)) / pow(yPy, 2)) / 2;
    Rcpp::Named("l2.ML") = l2;
  }else{
    l2 = (tr_HinvG2 - n * (2 * yPGPGPy * yPy - pow(yPGPy, 2)) / pow(yPy, 2)) / 2;
    Rcpp::Named("l2.ML") = l2;
  }
  return(l2);
}


Rcpp::List ml_est(double lambda, MatrixXd Y, MatrixXd X, MatrixXd U, MatrixXd D, MatrixXd Z,
                  MatrixXd K, double logXtX, double conv_param = 1e-06, int count_max = 15,
                  double bounds1 = 1e-06, double bounds2 = 1e06, bool REML = true){

  const int n(Y.rows()), p(X.cols()), ucol(U.cols());

  // Set initial values and prepare empty boxes.
  int count = 0;
  VectorXd lambdas = VectorXd::Zero(count_max + 1, 1);
  lambdas[0] = lambda;
  VectorXd ls = VectorXd::Zero(count_max, 1);
  VectorXd l1s = VectorXd::Zero(count_max, 1);
  VectorXd l2s = VectorXd::Zero(count_max, 1);
  bool l_conv = false;

  while(!(l_conv || count >= count_max)){
    MatrixXd EV = lambdas[count] * D + MatrixXd::Ones(ucol, 1);

    // Calculate aHinvb.series (The components of log-likelifood, score and hessian matrix)
    Rcpp::List aPb_series_res = aPb_series(Y, Y, U, EV, X);

    MapVec aPbs = Rcpp::as<MapVec>(aPb_series_res[0]);
    MapVec aPPbs = Rcpp::as<MapVec>(aPb_series_res[1]);
    MapVec aPPPbs = Rcpp::as<MapVec>(aPb_series_res[2]);
    MapVec tr_Ps = Rcpp::as<MapVec>(aPb_series_res[3]);
    MapVec tr_PPs = Rcpp::as<MapVec>(aPb_series_res[4]);
    MapMat XtHinvX = Rcpp::as<MapMat>(aPb_series_res[5]);

    double logH = (EV.array().log()).sum();
    SelfAdjointEigenSolver<MatrixXd> es(XtHinvX);
    double logXtHinvX = (es.eigenvalues().array().log()).sum();

    double tr_HinvG = (n - tr_Ps[0]) / lambdas[count];
    double tr_HinvG2 = (n + tr_PPs[0] - 2 * tr_Ps[0]) / pow(lambdas[count], 2);

    double tr_PG = (n - p - tr_Ps[p]) / lambdas[count];
    double tr_PG2 = (n - p + tr_PPs[p] - 2 * tr_Ps[p]) / pow(lambdas[count], 2);

    double yPy = aPbs[p];
    double yPGPy = (aPbs[p] - aPPbs[p]) / lambdas[count];
    double yPGPGPy = (aPbs[p] + aPPPbs[p] - 2 * aPPbs[p]) / pow(lambdas[count], 2);


    // Calculate log-likelihood, score and hessian ###
    double l = llik(n, p, logH, yPy, logXtX, logXtHinvX, REML);
    double l1 = score_func(n, p, tr_HinvG, yPy, yPGPy, tr_PG, REML);
    double l2 = hess_func(n, p, tr_HinvG2, yPy, yPGPy, yPGPGPy, tr_PG2, REML);

    ls[count] = l;
    l1s[count] = l1;
    l2s[count] = l2;

    // Update the parameter lambda
    double lambda_next = lambdas[count] - l1 / l2;
    if((lambda_next < bounds1) || (lambda_next > bounds2)){
      Rcpp::NumericVector lambda_next2 = Rcpp::runif(1, bounds1, (lambdas.segment(0, count + 1)).mean());
      lambda_next = lambda_next2[0];
    }

    if(count >= 1){
      double l_correction = ls[count] - ls[count - 1];
      l_conv = std::abs(l_correction) <= conv_param;
      if(l_correction < 0){
        Rcpp::NumericVector lambda_next2 = Rcpp::runif(1, bounds1, (lambdas.segment(0, count + 1)).mean());
        lambda_next = lambda_next2[0];
      }
    }
    lambdas[count + 1] = lambda_next;
    count += 1;
  }


  VectorXd max_lsvec = VectorXd::Zero(1, 1);
  double max_ls = (ls.segment(0, count - 1)).maxCoeff();
  max_lsvec[0] = max_ls;
  Rcpp::NumericVector max_no = matching(max_lsvec, ls);
  double lambda_opt = lambdas[max_no[0] + 1];

  return Rcpp::List::create(Rcpp::Named("lambda") = Rcpp::wrap(lambda_opt),
                            Rcpp::Named("max_l") = Rcpp::wrap(max_ls),
                            Rcpp::Named("count") = Rcpp::wrap(count));
}



// The start of the output functions
// Calculate aHinvb
// [[Rcpp::export]]
Rcpp::List aHinvb_out(Rcpp::NumericMatrix a, Rcpp::NumericMatrix b,
                      Rcpp::NumericMatrix u, Rcpp::NumericMatrix ev)
{
  const MapMat A = Rcpp::as<MapMat>(a);
  const MapMat B = Rcpp::as<MapMat>(b);
  const MapMat U = Rcpp::as<MapMat>(u);
  const MapMat EV = Rcpp::as<MapMat>(ev);

  MatrixXd va = U.adjoint() * A;
  MatrixXd vb = U.adjoint() * B;

  MatrixXd EV2 = power(EV, 2);
  MatrixXd EV3 = power(EV, 3);


  MatrixXd aH1b = va.adjoint() * elediv(vb, EV, 1);
  MatrixXd aH2b = va.adjoint() * elediv(vb, EV2, 1);
  MatrixXd aH3b = va.adjoint() * elediv(vb, EV3, 1);

  return Rcpp::List::create(Rcpp::Named("aH1b") = aH1b,
                            Rcpp::Named("aH2b") = aH2b,
                            Rcpp::Named("aH3b") = aH3b);
}

// Calculate aPb, aPPb, aPPPPb, tr(P), tr(PP) (and P)
//
// [[Rcpp::export]]
Rcpp::List aPb_series_out(Rcpp::NumericMatrix a, Rcpp::NumericMatrix b,
                          Rcpp::NumericMatrix u, Rcpp::NumericMatrix ev,
                          Rcpp::NumericMatrix x)
{
  const MapMat A = Rcpp::as<MapMat>(a);
  const MapMat B = Rcpp::as<MapMat>(b);
  const MapMat U = Rcpp::as<MapMat>(u);
  const MapMat EV = Rcpp::as<MapMat>(ev);
  const MapMat X = Rcpp::as<MapMat>(x);
  // const MapMat H = Rcpp::as<MapMat>(h);

  const int p(X.cols()), m(U.cols());
  const MatrixXd ones = MatrixXd::Ones(m, 1);

  const MatrixXd aPb_0 = aHinvb(A, B, U, EV);
  const MatrixXd aPX_0 = aHinvb(A, X, U, EV);
  const MatrixXd bPX_0 = aHinvb(B, X, U, EV);
  const MatrixXd XPX_0 = aHinvb(X, X, U, EV);

  const MatrixXd XPXs_0 = reshape(XPX_0.block(0, 0, p, p), 1, p * p, 1);
  const MatrixXd XPPXs_0 = reshape(XPX_0.block(p, 0, p, p), 1, p * p, 1);
  const MatrixXd XPPPXs_0 = reshape(XPX_0.block(2 * p, 0, p, p), 1, p * p, 1);


  VectorXd aPbs = VectorXd::Zero(p + 1, 1);
  VectorXd aPPbs = VectorXd::Zero(p + 1, 1);
  VectorXd aPPPbs = VectorXd::Zero(p + 1, 1);

  VectorXd tr_Ps = VectorXd::Zero(p + 1, 1);
  VectorXd tr_PPs = VectorXd::Zero(p + 1, 1);

  MatrixXd aPXs = MatrixXd::Zero(p + 1, p);
  MatrixXd aPPXs = MatrixXd::Zero(p + 1, p);
  MatrixXd aPPPXs = MatrixXd::Zero(p + 1, p);
  MatrixXd bPXs = MatrixXd::Zero(p + 1, p);
  MatrixXd bPPXs = MatrixXd::Zero(p + 1, p);
  MatrixXd bPPPXs = MatrixXd::Zero(p + 1, p);

  MatrixXd XPXs = MatrixXd::Zero(p + 1, p * p);
  MatrixXd XPPXs = MatrixXd::Zero(p + 1, p * p);
  MatrixXd XPPPXs = MatrixXd::Zero(p + 1, p * p);

  aPbs[0] = aPb_0(0, 0);
  aPXs.row(0) = aPX_0.row(0);
  bPXs.row(0) = bPX_0.row(0);
  XPXs.row(0) = XPXs_0.row(0);

  aPPbs[0] = aPb_0(1, 0);
  aPPXs.row(0) = aPX_0.row(1);
  bPPXs.row(0) = bPX_0.row(1);
  XPPXs.row(0) = XPPXs_0.row(0);

  aPPPbs[0] = aPb_0(2, 0);
  aPPPXs.row(0) = aPX_0.row(2);
  bPPPXs.row(0) = bPX_0.row(2);
  XPPPXs.row(0) = XPPPXs_0.row(0);

  tr_Ps[0] = elediv(ones, EV).sum();
  tr_PPs[0] = elediv(ones, power(EV, 2)).sum();

  for(int i = 0; i < p; i++){
    // Calculate aPb
    double aPb_next = aPbs[i] - (aPXs(i, i) * bPXs(i, i)) / XPXs(i, i * (p + 1));
    aPbs[i + 1] = aPb_next;

    VectorXd aPX_next = aPXs.row(i) - (aPXs(i, i) * XPXs.block(i, i * p, 1, p)) / XPXs(i, i * (p + 1));
    aPXs.row(i + 1) = aPX_next;

    VectorXd bPX_next = bPXs.row(i) - (bPXs(i, i) * XPXs.block(i, i * p, 1, p)) / XPXs(i, i * (p + 1));
    bPXs.row(i + 1) = bPX_next;

    MatrixXd XPX_next = XPXs.row(i) - reshape(crossprod(XPXs.block(i, i * p, 1, p),
                                              XPXs.block(i, i * p, 1, p)), 1, p * p, 1) / XPXs(i, i * (p + 1));
    XPXs.row(i + 1) = XPX_next;


    // Calculate aPPbs
    double aPPb_next = aPPbs[i] + (aPXs(i, i) * bPXs(i, i) * XPPXs(i, i * (p + 1))) / pow(XPXs(i, i * (p + 1)), 2) -
      (aPXs(i, i) * bPPXs(i, i)) / XPXs(i, i * (p + 1)) - (aPPXs(i, i) * bPXs(i, i)) / XPXs(i, i * (p + 1));
    aPPbs[i + 1] = aPPb_next;

    VectorXd aPPXs_next = aPPXs.row(i) + (aPXs(i, i) * XPXs.block(i, i * p, 1, p) * XPPXs(i, i * (p + 1))) / pow(XPXs(i, i * (p + 1)), 2) -
      (aPXs(i, i) * XPPXs.block(i, i * p, 1, p)) / XPXs(i, i * (p + 1)) - (aPPXs(i, i) * XPXs.block(i, i * p, 1, p)) / XPXs(i, i * (p + 1));
    aPPXs.row(i + 1) = aPPXs_next;

    VectorXd bPPXs_next = bPPXs.row(i) + (bPXs(i, i) * XPXs.block(i, i * p, 1, p) * XPPXs(i, i * (p + 1))) / pow(XPXs(i, i * (p + 1)), 2) -
      (bPXs(i, i) * XPPXs.block(i, i * p, 1, p)) / XPXs(i, i * (p + 1)) - (bPPXs(i, i) * XPXs.block(i, i * p, 1, p)) / XPXs(i, i * (p + 1));
    bPPXs.row(i + 1) = bPPXs_next;

    MatrixXd XPPXs_next = XPPXs.row(i) +
      reshape(crossprod(XPXs.block(i, i * p, 1, p), XPXs.block(i, i * p, 1, p)), 1, p * p, 1) * XPPXs(i, i * (p + 1)) / pow(XPXs(i, i * (p + 1)), 2) -
      reshape(crossprod(XPXs.block(i, i * p, 1, p), XPPXs.block(i, i * p, 1, p)), 1, p * p, 1) / XPXs(i, i * (p + 1)) -
      reshape(crossprod(XPPXs.block(i, i * p, 1, p), XPXs.block(i, i * p, 1, p)), 1, p * p, 1) / XPXs(i, i * (p + 1));
    XPPXs.row(i + 1) = XPPXs_next;



    // Calculate aPPPb
    double aPPPb_next = aPPPbs[i] - (aPXs(i, i) * bPXs(i, i) * pow(XPPXs(i, i * (p + 1)), 2)) / pow(XPXs(i, i * (p + 1)), 3) -
      (aPXs(i, i) * bPPPXs(i, i)) / XPXs(i, i * (p + 1)) - (aPPPXs(i, i) * bPXs(i, i)) / XPXs(i, i * (p + 1)) -
      (aPPXs(i, i) * bPPXs(i, i)) / XPXs(i, i * (p + 1)) + (aPXs(i, i) * bPPXs(i, i) * XPPXs(i, i * (p + 1))) / pow(XPXs(i, i * (p + 1)), 2) +
      (aPPXs(i, i) * bPXs(i, i) * XPPXs(i, i * (p + 1))) / pow(XPXs(i, i * (p + 1)), 2) +
      (aPXs(i, i) * bPXs(i, i) * XPPPXs(i, i * (p + 1))) / pow(XPXs(i, i * (p + 1)), 2);
    aPPPbs[i + 1] = aPPPb_next;

    VectorXd aPPPX_next = aPPPXs.row(i) - (aPXs(i, i) * XPXs.block(i, i * p, 1, p) * pow(XPPXs(i, i * (p + 1)), 2)) / pow(XPXs(i, i * (p + 1)), 3) -
      (aPXs(i, i) * XPPPXs.block(i, i * p, 1, p)) / XPXs(i, i * (p + 1)) - (aPPPXs(i, i) * XPXs.block(i, i * p, 1, p)) / XPXs(i, i * (p + 1)) -
      (aPPXs(i, i) * XPPXs.block(i, i * p, 1, p)) / XPXs(i, i * (p + 1)) + (aPXs(i, i) * XPPXs.block(i, i * p, 1, p) * XPPXs(i, i * (p + 1))) / pow(XPXs(i, i * (p + 1)), 2) +
      (aPPXs(i, i) * XPXs.block(i, i * p, 1, p) * XPPXs(i, i * (p + 1))) / pow(XPXs(i, i * (p + 1)), 2) +
      (aPXs(i, i) * XPXs.block(i, i * p, 1, p) * XPPPXs(i, i * (p + 1))) / pow(XPXs(i, i * (p + 1)), 2);
    aPPPXs.row(i + 1) = aPPPX_next;

    VectorXd bPPPX_next = bPPPXs.row(i) - (bPXs(i, i) * XPXs.block(i, i * p, 1, p) * pow(XPPXs(i, i * (p + 1)), 2)) / pow(XPXs(i, i * (p + 1)), 3) -
      (bPXs(i, i) * XPPPXs.block(i, i * p, 1, p)) / XPXs(i, i * (p + 1)) - (bPPPXs(i, i) * XPXs.block(i, i * p, 1, p)) / XPXs(i, i * (p + 1)) -
      (bPPXs(i, i) * XPPXs.block(i, i * p, 1, p)) / XPXs(i, i * (p + 1)) + (bPXs(i, i) * XPPXs.block(i, i * p, 1, p) * XPPXs(i, i * (p + 1))) / pow(XPXs(i, i * (p + 1)), 2) +
      (bPPXs(i, i) * XPXs.block(i, i * p, 1, p) * XPPXs(i, i * (p + 1))) / pow(XPXs(i, i * (p + 1)), 2) +
      (bPXs(i, i) * XPXs.block(i, i * p, 1, p) * XPPPXs(i, i * (p + 1))) / pow(XPXs(i, i * (p + 1)), 2);
    bPPPXs.row(i + 1) = bPPPX_next;

    MatrixXd XPPPX_next = XPPPXs.row(i) -
      (reshape(crossprod(XPXs.block(i, i * p, 1, p), XPXs.block(i, i * p, 1, p)), 1, p * p, 1) * pow(XPPXs(i, i * (p + 1)), 2)) / pow(XPXs(i, i * (p + 1)), 3) -
      reshape(crossprod(XPXs.block(i, i * p, 1, p), XPPPXs.block(i, i * p, 1, p)), 1, p * p, 1) / XPXs(i, i * (p + 1)) -
      reshape(crossprod(XPPPXs.block(i, i * p, 1, p), XPXs.block(i, i * p, 1, p)), 1, p * p, 1)  / XPXs(i, i * (p + 1)) -
      reshape(crossprod(XPPXs.block(i, i * p, 1, p), XPPXs.block(i, i * p, 1, p)), 1, p * p, 1)  / XPXs(i, i * (p + 1)) +
      (reshape(crossprod(XPXs.block(i, i * p, 1, p), XPPXs.block(i, i * p, 1, p)), 1, p * p, 1) * XPPXs(i, i * (p + 1))) / pow(XPXs(i, i * (p + 1)), 2) +
      (reshape(crossprod(XPPXs.block(i, i * p, 1, p), XPXs.block(i, i * p, 1, p)), 1, p * p, 1)  * XPPXs(i, i * (p + 1))) / pow(XPXs(i, i * (p + 1)), 2) +
      (reshape(crossprod(XPXs.block(i, i * p, 1, p), XPXs.block(i, i * p, 1, p)), 1, p * p, 1)  * XPPPXs(i, i * (p + 1))) / pow(XPXs(i, i * (p + 1)), 2);
    XPPPXs.row(i + 1) = XPPPX_next;



    // Calculate trace(P)
    double tr_P_next = tr_Ps[i] - XPPXs(i, i * (p + 1)) / XPXs(i, i * (p + 1));
    tr_Ps[i + 1] = tr_P_next;

    // Calculate trace(PP)
    double tr_PP_next = tr_PPs[i] + pow(XPPXs(i, i * (p + 1)), 2) / pow(XPXs(i, i * (p + 1)), 2) -
      2 * XPPPXs(i, i * (p + 1)) / XPXs(i, i * (p + 1));
    tr_PPs[i + 1] = tr_PP_next;
  }

  return Rcpp::List::create(Rcpp::Named("aPb") = Rcpp::wrap(aPbs),
                            Rcpp::Named("aPPb") = Rcpp::wrap(aPPbs),
                            Rcpp::Named("aPPPb") = Rcpp::wrap(aPPPbs),
                            Rcpp::Named("tr.P") = Rcpp::wrap(tr_Ps),
                            Rcpp::Named("tr.PP") = Rcpp::wrap(tr_PPs),
                            Rcpp::Named("XtHinvX") = Rcpp::wrap(XPX_0.block(0, 0, p, p)));
}


// Calculate log-likelihood
//
// [[Rcpp::export]]
double llik_out(int n, int p, double logH, double yPy, double logXtX, double logXtHinvX, bool REML = true){
  double pi = 3.14159;
  double l = 0;

  if(REML){
    l = ((n - p) * std::log((n - p) / (2 * pi)) - (n - p) + logXtX -
      logH - logXtHinvX -  (n - p) * std::log(yPy)) / 2;
    Rcpp::Named("l.REML") = l;
  }else{
    l = (n * std::log(n / (2 * pi)) - n - logH - n * std::log(yPy)) / 2;
    Rcpp::Named("l.ML") = l;
  }
  return(l);
}

// Calculate the score function
// [[Rcpp::export]]
double score_func_out(int n, int p, double tr_HinvG, double yPy,
                      double yPGPy, double tr_PG, bool REML = true){
  double l1 = 0;
  if(REML){
    l1 = (- tr_PG + (n - p) * yPGPy / yPy) / 2;
    Rcpp::Named("l1.REML") = l1;
  }else{
    l1 = (- tr_HinvG + n * yPGPy / yPy) / 2;
    Rcpp::Named("l1.ML") = l1;
  }
  return(l1);
}

// Calculate the score function
//
//
// [[Rcpp::export]]
double hess_func_out(int n, int p, double tr_HinvG2, double yPy, double yPGPy,
                     double yPGPGPy, double tr_PG2, bool REML = true){
  double l2 = 0;
  if(REML){
    l2 = (tr_PG2 - (n - p) * (2 * yPGPGPy * yPy - pow(yPGPy, 2)) / pow(yPy, 2)) / 2;
    Rcpp::Named("l2.ML") = l2;
  }else{
    l2 = (tr_HinvG2 - n * (2 * yPGPGPy * yPy - pow(yPGPy, 2)) / pow(yPy, 2)) / 2;
    Rcpp::Named("l2.ML") = l2;
  }
  return(l2);
}



// Maximize likelihood when given lambda_0
//
// [[Rcpp::export]]
Rcpp::List ml_est_out(double lambda, Rcpp::NumericMatrix y, Rcpp::NumericMatrix x,
                      Rcpp::NumericMatrix u, Rcpp::NumericMatrix delta, Rcpp::NumericMatrix z,
                      Rcpp::NumericMatrix k, double logXtX, double conv_param = 1e-06,
                      int count_max = 15, double bounds1 = 1e-06, double bounds2 = 1e06, bool REML = true){

  const MapMat Y = Rcpp::as<MapMat>(y);
  const MapMat X = Rcpp::as<MapMat>(x);
  const MapMat U = Rcpp::as<MapMat>(u);
  const MapMat D = Rcpp::as<MapMat>(delta);

  const int n(Y.rows()), p(X.cols()), ucol(U.cols());

  // Set initial values and prepare empty boxes.
  int count = 0;
  VectorXd lambdas = VectorXd::Zero(count_max + 1, 1);
  lambdas[0] = lambda;
  VectorXd ls = VectorXd::Zero(count_max, 1);
  VectorXd l1s = VectorXd::Zero(count_max, 1);
  VectorXd l2s = VectorXd::Zero(count_max, 1);
  bool l_conv = false;

  while(!(l_conv || count >= count_max)){
    MatrixXd EV = lambdas[count] * D + MatrixXd::Ones(ucol, 1);

    // Calculate aHinvb.series (The components of log-likelifood, score and hessian matrix)
    Rcpp::List aPb_series_res = aPb_series(Y, Y, U, EV, X);

    MapVec aPbs = Rcpp::as<MapVec>(aPb_series_res[0]);
    MapVec aPPbs = Rcpp::as<MapVec>(aPb_series_res[1]);
    MapVec aPPPbs = Rcpp::as<MapVec>(aPb_series_res[2]);
    MapVec tr_Ps = Rcpp::as<MapVec>(aPb_series_res[3]);
    MapVec tr_PPs = Rcpp::as<MapVec>(aPb_series_res[4]);
    MapMat XtHinvX = Rcpp::as<MapMat>(aPb_series_res[5]);

    double logH = (EV.array().log()).sum();
    SelfAdjointEigenSolver<MatrixXd> es(XtHinvX);
    double logXtHinvX = (es.eigenvalues().array().log()).sum();

    double tr_HinvG = (n - tr_Ps[0]) / lambdas[count];
    double tr_HinvG2 = (n + tr_PPs[0] - 2 * tr_Ps[0]) / pow(lambdas[count], 2);

    double tr_PG = (n - p - tr_Ps[p]) / lambdas[count];
    double tr_PG2 = (n - p + tr_PPs[p] - 2 * tr_Ps[p]) / pow(lambdas[count], 2);

    double yPy = aPbs[p];
    double yPGPy = (aPbs[p] - aPPbs[p]) / lambdas[count];
    double yPGPGPy = (aPbs[p] + aPPPbs[p] - 2 * aPPbs[p]) / pow(lambdas[count], 2);


    // Calculate log-likelihood, score and hessian ###
    double l = llik(n, p, logH, yPy, logXtX, logXtHinvX, REML);
    double l1 = score_func(n, p, tr_HinvG, yPy, yPGPy, tr_PG, REML);
    double l2 = hess_func(n, p, tr_HinvG2, yPy, yPGPy, yPGPGPy, tr_PG2, REML);

    ls[count] = l;
    l1s[count] = l1;
    l2s[count] = l2;

    // Update the parameter lambda
    double lambda_next = lambdas[count] - l1 / l2;
    if((lambda_next < bounds1) || (lambda_next > bounds2)){
      Rcpp::NumericVector lambda_next2 = Rcpp::runif(1, bounds1, (lambdas.segment(0, count + 1)).mean());
      lambda_next = lambda_next2[0];
    }

    if(count >= 1){
      double l_correction = ls[count] - ls[count - 1];
      l_conv = std::abs(l_correction) <= conv_param;
      if(l_correction < 0){
        Rcpp::NumericVector lambda_next2 = Rcpp::runif(1, bounds1, (lambdas.segment(0, count + 1)).mean());
        lambda_next = lambda_next2[0];
      }
    }
    lambdas[count + 1] = lambda_next;
    count += 1;
  }


  VectorXd max_lsvec = VectorXd::Zero(1, 1);
  double max_ls = (ls.segment(0, count - 1)).maxCoeff();
  max_lsvec[0] = max_ls;
  Rcpp::NumericVector max_no = matching(max_lsvec, ls);
  double lambda_opt = lambdas[max_no[0] + 1];

  return Rcpp::List::create(Rcpp::Named("lambda") = Rcpp::wrap(lambda_opt),
                            Rcpp::Named("max_l") = Rcpp::wrap(max_ls),
                            Rcpp::Named("count") = Rcpp::wrap(count));
}

// Eigen decomposition via Rcpp
//
// [[Rcpp::export]]
Rcpp::List eigen_out(Rcpp::NumericMatrix a, bool symmetric = true, bool ainv = false){
  const MapMat A = Rcpp::as<MapMat>(a);
  Rcpp::Function f("eigen");

  Rcpp::List sol = f(Rcpp::wrap(A), symmetric);
  MapMat D = Rcpp::as<MapMat>(sol[0]);
  MapMat U = Rcpp::as<MapMat>(sol[1]);

  const int m(U.cols());
  if(!ainv){
    return Rcpp::List::create(Rcpp::Named("values") = Rcpp::wrap(D),
                              Rcpp::Named("vectors") = Rcpp::wrap(U));
  }else{
    MatrixXd D2 = MatrixXd::Zero(1, m);
    D2.row(0) = D.col(0);
    MatrixXd UdivD2 = elediv(U, D2, 2);

    MatrixXd Ainv = tcrossprod(UdivD2, U);

    return Rcpp::List::create(Rcpp::Named("values") = Rcpp::wrap(D),
                              Rcpp::Named("vectors") = Rcpp::wrap(U),
                              Rcpp::Named("inverse") = Rcpp::wrap(Ainv));
  }
}



// Calculate some values when given lambda
//
// [[Rcpp::export]]
Rcpp::List ml_one_step(double lambda, Rcpp::NumericMatrix y, Rcpp::NumericMatrix x,
                       Rcpp::NumericMatrix u, Rcpp::NumericMatrix delta, Rcpp::NumericMatrix z,
                       Rcpp::NumericMatrix k, double logXtX, double conv_param = 1e-06,
                       int count_max = 15, double bounds1 = 1e-06, double bounds2 = 1e06,
                       bool REML = true, bool SE = false, bool return_Hinv = true){


  const MapMat Y = Rcpp::as<MapMat>(y);
  const MapMat X = Rcpp::as<MapMat>(x);
  const MapMat U = Rcpp::as<MapMat>(u);
  const MapMat D = Rcpp::as<MapMat>(delta);
  const MapMat Z = Rcpp::as<MapMat>(z);
  const MapMat K = Rcpp::as<MapMat>(k);

  const int n(Y.rows()), p(X.cols()), ucol(U.cols());


  MatrixXd EV = lambda * D + MatrixXd::Ones(ucol, 1);

  // Calculate aHinvb.series (The components of log-likelifood, score and hessian matrix)
  Rcpp::List aPb_series_res = aPb_series(Y, Y, U, EV, X);

  MapVec aPbs = Rcpp::as<MapVec>(aPb_series_res[0]);
  MapMat XtHinvX = Rcpp::as<MapMat>(aPb_series_res[5]);

  double logH = (EV.array().log()).sum();
  SelfAdjointEigenSolver<MatrixXd> es(XtHinvX);
  double logXtHinvX = (es.eigenvalues().array().log()).sum();

  double yPy = aPbs[p];



  // Calculate log-likelihood, score and Hinv2 ###
  double LL = llik(n, p, logH, yPy, logXtX, logXtHinvX, REML);

  MatrixXd EV2 = EV.transpose();
  MatrixXd UdivEV2 = elediv(U, EV2, 2);
  MatrixXd Hinv2 = tcrossprod(UdivEV2, U);

  double Ve = 0;
  if(REML){
    Ve = yPy / (n - p);
  }else{
    Ve = yPy / n;
  }
  double Vu = Ve * lambda;

  MatrixXd EV_EMMA = D + MatrixXd::Ones(ucol, 1) / lambda;
  MatrixXd EV_EMMA2 = EV_EMMA.transpose();
  MatrixXd UdivEV_EMMA2 = elediv(U, EV_EMMA2, 2);
  MatrixXd Hinv = tcrossprod(UdivEV_EMMA2, U);

  MatrixXd beta_mat = solve(XtHinvX, crossprod(X, Hinv2 * Y));
  VectorXd beta = beta_mat.col(0);

  MatrixXd KZt = tcrossprod(K, Z);
  MatrixXd KZt_Hinv = KZt * Hinv;
  MatrixXd u_mat = (KZt_Hinv * (Y - X * beta_mat));
  VectorXd u_vec = u_mat.col(0);


  if(!SE){
    if(return_Hinv){
      return Rcpp::List::create(Rcpp::Named("Vu") = Rcpp::wrap(Vu),
                                Rcpp::Named("Ve") = Rcpp::wrap(Ve),
                                Rcpp::Named("beta") = Rcpp::wrap(beta),
                                Rcpp::Named("u") = Rcpp::wrap(u_vec),
                                Rcpp::Named("LL") = Rcpp::wrap(LL),
                                Rcpp::Named("Hinv") = Rcpp::wrap(Hinv),
                                Rcpp::Named("Hinv2") = Rcpp::wrap(Hinv2));
    }
    else {
      return Rcpp::List::create(Rcpp::Named("Vu") = Rcpp::wrap(Vu),
                                Rcpp::Named("Ve") = Rcpp::wrap(Ve),
                                Rcpp::Named("beta") = Rcpp::wrap(beta),
                                Rcpp::Named("u") = Rcpp::wrap(u_vec),
                                Rcpp::Named("LL") = Rcpp::wrap(LL));
    }
  }
  else {
    MatrixXd W = crossprod(X, Hinv * X);
    MatrixXd Winv = inv(W);
    MatrixXd beta_SE_mat = (Vu * diagelements(W)).array().sqrt();
    VectorXd beta_SE = beta_SE_mat.col(0);

    MatrixXd WW = tcrossprod(KZt_Hinv, KZt);
    MatrixXd WWW = KZt_Hinv * X;
    MatrixXd  u_SE_mat = (Vu * (diagelements(K) - diagelements(WW) +
      diagelements(tcrossprod(WWW * Winv, WWW)))).array().sqrt();
    VectorXd u_SE = u_SE_mat.col(0);

    if(return_Hinv){
      return Rcpp::List::create(Rcpp::Named("Vu") = Rcpp::wrap(Vu),
                                Rcpp::Named("Ve") = Rcpp::wrap(Ve),
                                Rcpp::Named("beta") = Rcpp::wrap(beta),
                                Rcpp::Named("beta.SE") = Rcpp::wrap(beta_SE),
                                Rcpp::Named("u") = Rcpp::wrap(u_vec),
                                Rcpp::Named("u.SE") = Rcpp::wrap(u_SE),
                                Rcpp::Named("LL") = Rcpp::wrap(LL),
                                Rcpp::Named("Hinv") = Rcpp::wrap(Hinv),
                                Rcpp::Named("Hinv2") = Rcpp::wrap(Hinv2));
    }
    else {
      return Rcpp::List::create(Rcpp::Named("Vu") = Rcpp::wrap(Vu),
                                Rcpp::Named("Ve") = Rcpp::wrap(Ve),
                                Rcpp::Named("beta") = Rcpp::wrap(beta),
                                Rcpp::Named("beta.SE") = Rcpp::wrap(beta_SE),
                                Rcpp::Named("u") = Rcpp::wrap(u_vec),
                                Rcpp::Named("u.SE") = Rcpp::wrap(u_SE),
                                Rcpp::Named("LL") = Rcpp::wrap(LL));
    }
  }
}




// The spectral decomposition of G matrix (cholesky decomposition)
//
// [[Rcpp::export]]
Rcpp::List spectralG_cholesky(Rcpp::NumericMatrix zbt, Rcpp::NumericMatrix x,
                              bool return_G = true, bool return_SGS = false){
  const MapMat X = Rcpp::as<MapMat>(x);
  const MapMat ZBt = Rcpp::as<MapMat>(zbt);
  const int n(X.rows()), m(ZBt.cols()), p(X.cols());
  VectorXd ones_vec = VectorXd::Ones(n, 1);
  MatrixXd I = MatrixXd::Identity(n, n);
  MatrixXd U;
  MatrixXd D;
  MatrixXd Q;
  MatrixXd theta;



  if(return_G){
    Rcpp::List svd_ZBt = svd(ZBt, n, m);
    MapMat U2 = Rcpp::as<MapMat>(svd_ZBt[1]);
    U = U2;

    MatrixXd D20 = MatrixXd::Zero(m, 1);
    MapVec D2_vec = Rcpp::as<MapVec>(svd_ZBt[0]);
    D20.col(0) = D2_vec;
    MatrixXd D2 = power(D20, 2);
    D = rbind(D2, MatrixXd::Zero(n - m, 1));
  }else{
    U = MatrixXd::Zero(n, m);
    D = MatrixXd::Zero(n, 1);
  }

  if(return_SGS){
    MatrixXd XtX = crossprod(X, X);
    double rank_X = XtX.colPivHouseholderQr().rank();
    if(rank_X < p){
      Rcpp::stop("X not full rank");
    }

    MatrixXd XtXinv = inv(XtX);
    MatrixXd S =  I - tcrossprod(X * XtXinv, X);

    MatrixXd  SZBt = S * ZBt;

    Rcpp::List  svd_SZBt = Rcpp::List::create(Rcpp::Named("d") = MatrixXd::Zero(m, 1),
                                              Rcpp::Named("u") = MatrixXd::Zero(n, m),
                                              Rcpp::Named("v") = MatrixXd::Zero(m, m));

    svd_SZBt = svd(SZBt, m, m);

    MatrixXd svd_SZBt_U = Rcpp::as<MapMat>(svd_SZBt[1]);
    HouseholderQR<MatrixXd> QR(cbind(X, svd_SZBt_U));
    MatrixXd Q2 = QR.householderQ();
    Q = Q2.block(0, p, n, n - p);


    MatrixXd  R2 = QR.matrixQR();
    MatrixXd  R0 = R2.block(p, p, m, m);

    TriangularView<MatrixXd, Upper> R0_up(R0);
    MatrixXd R = R0_up.toDenseMatrix();

    MatrixXd svd_SZBt_D = Rcpp::as<MapMat>(svd_SZBt[0]);
    MatrixXd ans = solve(power(R, 2).transpose(), power(svd_SZBt_D, 2));
    theta = rbind(ans.block(0, 0, m - p, 1), MatrixXd::Zero(n - m, 1));
  }else{
    Q = MatrixXd::Zero(n, m);
    theta = MatrixXd::Zero(n, 1);
  }

  return Rcpp::List::create(Rcpp::Named("n") = Rcpp::wrap(n),
                            Rcpp::Named("p") = Rcpp::wrap(p),
                            Rcpp::Named("U") = Rcpp::wrap(U),
                            Rcpp::Named("D") = Rcpp::wrap(D),
                            Rcpp::Named("Q") = Rcpp::wrap(Q),
                            Rcpp::Named("theta") = Rcpp::wrap(theta));
}






// The spectral decomposition of G matrix (eigen decomposition)
//
// [[Rcpp::export]]
Rcpp::List spectralG_eigen(Rcpp::NumericMatrix zkzt, Rcpp::NumericMatrix x,
                           bool return_G = true, bool return_SGS = false){
  const MapMat X = Rcpp::as<MapMat>(x);
  const MapMat ZKZt = Rcpp::as<MapMat>(zkzt);
  const double n(X.rows()), p(X.cols());
  MatrixXd I = MatrixXd::Identity(n, n);
  MatrixXd U;
  MatrixXd D;
  MatrixXd Q;
  MatrixXd theta;


  double offset = std::sqrt(n);
  MatrixXd Hb = ZKZt + offset * I;

  if(return_G){
    Rcpp::List eigen_Hb = eigen(Hb, true, false);
    MatrixXd eigen_Hb_ev =  eigen_Hb[0];
    D = eigen_Hb_ev - offset * MatrixXd::Ones(n, 1);

    double Dmin = D.minCoeff();
    if(Dmin < -1e-04){
      Rcpp::stop("K not positive semi-definite.");
    }
    U = eigen_Hb[1];
  }else{
    U = MatrixXd::Zero(n, n);
    D = MatrixXd::Zero(n, 1);
  }

  if(return_SGS){
    MatrixXd XtX = crossprod(X, X);
    double rank_X = XtX.colPivHouseholderQr().rank();

    if(rank_X < p){
      Rcpp::stop("X not full rank");
    }

    MatrixXd XtXinv = inv(XtX);

    MatrixXd S =  I - tcrossprod(X * XtXinv, X);

    MatrixXd SHbS = S * Hb * S;
    Rcpp::List SHbS_system = eigen(SHbS, true, false);
    MatrixXd SHbS_system_ev = SHbS_system[0];
    theta = SHbS_system_ev.block(0, 0, n - p, 1) - offset * MatrixXd::Ones(n - p, 1);
    MatrixXd Q0 = SHbS_system[1];
    Q = Q0.block(0, 0, n, n - p);
  }else{
    Q = MatrixXd::Zero(n, n);
    theta = MatrixXd::Zero(n, 1);
  }

  return Rcpp::List::create(Rcpp::Named("n") = Rcpp::wrap(n),
                            Rcpp::Named("p") = Rcpp::wrap(p),
                            Rcpp::Named("U") = Rcpp::wrap(U),
                            Rcpp::Named("D") = Rcpp::wrap(D),
                            Rcpp::Named("Q") = Rcpp::wrap(Q),
                            Rcpp::Named("theta") = Rcpp::wrap(theta));
}







// Calculate Vu, Ve, LL and Hinv when given lambda
//
// [[Rcpp::export]]
Rcpp::List EMM2_last_step(double lambda, Rcpp::NumericMatrix y,
                          Rcpp::NumericMatrix x, Rcpp::NumericMatrix u,
                          Rcpp::NumericMatrix d, Rcpp::NumericMatrix z,
                          Rcpp::NumericMatrix k, Rcpp::NumericMatrix th,
                          Rcpp::NumericMatrix ome_sq, Rcpp::NumericMatrix ph,
                          double maxval, int n, int p, int df, bool SE = false,
                          bool return_Hinv = false){
  const MapMat Y = Rcpp::as<MapMat>(y);
  const MapMat X = Rcpp::as<MapMat>(x);
  const MapMat U = Rcpp::as<MapMat>(u);
  const MapMat D = Rcpp::as<MapMat>(d);
  const MapMat Z = Rcpp::as<MapMat>(z);
  const MapMat K = Rcpp::as<MapMat>(k);
  const MapMat theta = Rcpp::as<MapMat>(th);
  const MapMat omega_sq = Rcpp::as<MapMat>(ome_sq);
  const MapMat phi = Rcpp::as<MapMat>(ph);

  const int m(D.rows());
  MatrixXd ones1 = MatrixXd::Ones(m, 1);
  MatrixXd ones2 = MatrixXd::Ones(df, 1);

  double pi = 3.14159;
  double Vu_opt = elediv(omega_sq, theta + lambda * ones2).array().sum() / df;
  double Ve_opt = lambda * Vu_opt;

  MatrixXd EV_EMMA = phi + lambda * ones1;
  MatrixXd EV_EMMA2 = EV_EMMA.transpose();
  MatrixXd UdivEV_EMMA2 = elediv(U, EV_EMMA2, 2);
  MatrixXd Hinv = tcrossprod(UdivEV_EMMA2, U);

  MatrixXd W = crossprod(X, Hinv * X);
  MatrixXd beta_mat = solve(W, crossprod(X, Hinv * Y));
  VectorXd beta = beta_mat.col(0);

  MatrixXd KZt = tcrossprod(K, Z);
  MatrixXd KZt_Hinv = KZt * Hinv;
  MatrixXd u_mat = (KZt_Hinv * (Y - X * beta_mat));
  VectorXd u_vec = u_mat.col(0);

  double  LL = -0.5 * (maxval + df + df * std::log(2 * pi / df));

  if(!SE){
    if(return_Hinv){
      return Rcpp::List::create(Rcpp::Named("Vu") = Rcpp::wrap(Vu_opt),
                                Rcpp::Named("Ve") = Rcpp::wrap(Ve_opt),
                                Rcpp::Named("beta") = Rcpp::wrap(beta),
                                Rcpp::Named("u") = Rcpp::wrap(u_vec),
                                Rcpp::Named("LL") = Rcpp::wrap(LL),
                                Rcpp::Named("Hinv") = Rcpp::wrap(Hinv));
    }
    else {
      return Rcpp::List::create(Rcpp::Named("Vu") = Rcpp::wrap(Vu_opt),
                                Rcpp::Named("Ve") = Rcpp::wrap(Ve_opt),
                                Rcpp::Named("beta") = Rcpp::wrap(beta),
                                Rcpp::Named("u") = Rcpp::wrap(u_vec),
                                Rcpp::Named("LL") = Rcpp::wrap(LL));
    }
  }
  else {
    MatrixXd Winv = inv(W);
    MatrixXd beta_SE_mat = (Vu_opt * diagelements(W)).array().sqrt();
    VectorXd beta_SE = beta_SE_mat.col(0);

    MatrixXd WW = tcrossprod(KZt_Hinv, KZt);
    MatrixXd WWW = KZt_Hinv * X;
    MatrixXd  u_SE_mat = (Vu_opt * (diagelements(K) - diagelements(WW) +
      diagelements(tcrossprod(WWW * Winv, WWW)))).array().sqrt();
    VectorXd u_SE = u_SE_mat.col(0);

    if(return_Hinv){
      return Rcpp::List::create(Rcpp::Named("Vu") = Rcpp::wrap(Vu_opt),
                                Rcpp::Named("Ve") = Rcpp::wrap(Ve_opt),
                                Rcpp::Named("beta") = Rcpp::wrap(beta),
                                Rcpp::Named("beta.SE") = Rcpp::wrap(beta_SE),
                                Rcpp::Named("u") = Rcpp::wrap(u_vec),
                                Rcpp::Named("u.SE") = Rcpp::wrap(u_SE),
                                Rcpp::Named("LL") = Rcpp::wrap(LL),
                                Rcpp::Named("Hinv") = Rcpp::wrap(Hinv));
    }
    else {
      return Rcpp::List::create(Rcpp::Named("Vu") = Rcpp::wrap(Vu_opt),
                                Rcpp::Named("Ve") = Rcpp::wrap(Ve_opt),
                                Rcpp::Named("beta") = Rcpp::wrap(beta),
                                Rcpp::Named("beta.SE") = Rcpp::wrap(beta_SE),
                                Rcpp::Named("u") = Rcpp::wrap(u_vec),
                                Rcpp::Named("u.SE") = Rcpp::wrap(u_SE),
                                Rcpp::Named("LL") = Rcpp::wrap(LL));
    }
  }
}


// The kernel function of EM3.cpp
//
// [[Rcpp::export]]
Rcpp::List EM3_kernel(Rcpp::NumericMatrix y0, Rcpp::NumericMatrix X0, Rcpp::NumericMatrix ZKZt0,
                      Rcpp::NumericMatrix S0, Rcpp::NumericMatrix spI0,
                      double n, double p, bool REML = true){
  MapMat y = Rcpp::as<MapMat>(y0);
  MapMat ZKZt = Rcpp::as<MapMat>(ZKZt0);
  MapMat S = Rcpp::as<MapMat>(S0);
  MapMat spI = Rcpp::as<MapMat>(spI0);

  MatrixXd phi = MatrixXd::Zero(n, 1);
  double offset = std::sqrt(n);

  MatrixXd  ZKZtandoffset = ZKZt + offset * spI;


  MatrixXd SZKZtSandoffset = S * ZKZtandoffset * S;
  Rcpp::List svdSZKZtSandspI = eigen(SZKZtSandoffset, true, false);
  MapMat Ur0 = Rcpp::as<MapMat>(svdSZKZtSandspI[1]);
  MatrixXd Ur = Ur0.block(0, 0, n, n - p);

  MapMat lambda0 = Rcpp::as<MapMat>(svdSZKZtSandspI[0]);
  MatrixXd lambda = lambda0.block(0, 0, n - p, 1) - offset * MatrixXd::Ones(n - p, 1);
  MatrixXd eta = crossprod(Ur, y);

  if(!REML){
    Rcpp::List svdZKZtandoffset = eigen(ZKZtandoffset, true, false);
    MapMat phi0 = Rcpp::as<MapMat>(svdZKZtandoffset[0]);
    phi = phi0;
  }

  return Rcpp::List::create(Rcpp::Named("eta") = Rcpp::wrap(eta),
                            Rcpp::Named("lambda") = Rcpp::wrap(lambda),
                            Rcpp::Named("phi") = Rcpp::wrap(phi));
}


// Calculate P or Hinv, and log(P) or log(Hinv)
//
// [[Rcpp::export]]
Rcpp::List P_calc(double lambda, Rcpp::List Ws, Rcpp::List Gammas,
                  Rcpp::NumericMatrix u, Rcpp::NumericMatrix d,
                  int lw, bool Pornot = true, bool gammas_diag = true){
  const MapMat U = Rcpp::as<MapMat>(u);
  const MapMat D = Rcpp::as<MapMat>(d);

  const int n(U.rows()), n_p(U.cols());
  MatrixXd ones = MatrixXd::Ones(n, 1);
  MatrixXd ones1 = MatrixXd::Ones(n, 1);
  MatrixXd ones2 = MatrixXd::Ones(n_p, 1);

  if(Pornot){
    ones = ones2;
  }else{
    ones = ones1;
  }

  MatrixXd EV_EMMA = D + ones * lambda;
  MatrixXd EV_EMMA2 = EV_EMMA.transpose();
  MatrixXd UdivEV_EMMA2 = elediv(U, EV_EMMA2, 2);
  MatrixXd P0 = tcrossprod(UdivEV_EMMA2, U);

  MatrixXd P = P0;

  double lnP0 = -EV_EMMA.array().log().sum();
  double lnP = lnP0;


  for(int i = 0; i < lw; i++){
    MatrixXd Gamma = Gammas[i];
    const int k(Gamma.rows());
    MatrixXd Gamma_inv;
    if(gammas_diag){
      MatrixXd diag_Gamma = diagelements(Gamma);
      VectorXd diag_Gamma_inv = elediv(MatrixXd::Ones(k, 1), diag_Gamma, 0).col(0);
      Gamma_inv = diag_Gamma_inv.asDiagonal();
    }else{
      Gamma_inv = inv(Gamma);
    }

    MatrixXd W = Ws[i];
    MatrixXd PW = P * W;
    MatrixXd kkmat_inv = inv(Gamma_inv + crossprod(W,  PW));

    MatrixXd P_next = P - tcrossprod(PW * kkmat_inv, PW);
    // P_next = P - PW * kkmat_inv * t(W) * P;

    MatrixXd lnPfor_mat = MatrixXd::Identity(k, k) + crossprod(W, PW * Gamma);
    Rcpp::List eigen_res = eigen(lnPfor_mat, false, false);
    MatrixXd eigen_val = eigen_res[0];
    double lnP_part = -eigen_val.array().log().sum();

    P = P_next;
    lnP = lnP + lnP_part;
  }

  return(Rcpp::List::create(Rcpp::Named("P") = Rcpp::wrap(P),
                            Rcpp::Named("lnP") = Rcpp::wrap(lnP)));
}


// Calculate restricted log likelihood using P
//
// [[Rcpp::export]]
double llik_REML(double n, double p, double yPy, double lnP){
  const double pi = 3.14159;
  return(((n - p) * (std::log(n - p) - std::log(2 * pi) - 1 - std::log(yPy)) + lnP) / 2);
}


// Calculate log likelihood using P
//
// [[Rcpp::export]]
double llik_ML(double n, double yPy, double lnHinv){
  const double pi = 3.14159;
  return((n * (std::log(n) - std::log(2 * pi) - 1 - std::log(yPy)) + lnHinv) / 2);
}


// Calculate first derivariate of the log-likelihood for score test
//
// [[Rcpp::export]]
Rcpp::NumericVector score_l1(Rcpp::NumericMatrix y, Rcpp::NumericMatrix p0,
                             Rcpp::List Gs, int lg){
  const MapMat Y = Rcpp::as<MapMat>(y);
  const MapMat P0 = Rcpp::as<MapMat>(p0);

  MatrixXd P0Y = P0 * Y;
  VectorXd l1s = VectorXd::Zero( lg);


  for(int i = 0; i < lg; i++){
    MatrixXd G = Gs[i];
    MatrixXd ytP0GP0y = crossprod(P0Y, G * P0Y);
    double l1s_part = (ytP0GP0y(0, 0) - diagelements(P0 * G).sum()) / 2;
    l1s[i] = l1s_part;
  }

  return(Rcpp::wrap(l1s));
}



// Calculate first derivariate of the log-likelihood for score test (for linear kernel)
//
// [[Rcpp::export]]
Rcpp::NumericVector score_l1_linker(Rcpp::NumericMatrix y, Rcpp::NumericMatrix p0,
                                    Rcpp::List Ws, Rcpp::List Gammas, int lw){
  const MapMat Y = Rcpp::as<MapMat>(y);
  const MapMat P0 = Rcpp::as<MapMat>(p0);

  MatrixXd P0Y = P0 * Y;
  VectorXd l1s = VectorXd::Zero(lw);


  for(int i = 0; i < lw; i++){
    MatrixXd W = Ws[i];
    MatrixXd Gamma = Gammas[i];
    MatrixXd WP0Y = crossprod(W, P0Y);
    MatrixXd ytP0GP0y = crossprod(WP0Y, Gamma * WP0Y);
    double l1s_part = (ytP0GP0y(0, 0) - diagelements(crossprod(W, P0 * W * Gamma)).sum()) / 2;
    l1s[i] = l1s_part;
  }

  return(Rcpp::wrap(l1s));
}



// Calculate first derivariate of the log-likelihood for score test (for linear kernel)
//
// [[Rcpp::export]]
Rcpp::NumericVector score_l1_linker_diag(Rcpp::NumericMatrix y, Rcpp::NumericMatrix p0,
                                         Rcpp::List W2s, int lw){
  const MapMat Y = Rcpp::as<MapMat>(y);
  const MapMat P0 = Rcpp::as<MapMat>(p0);

  MatrixXd P0Y = P0 * Y;
  VectorXd l1s = VectorXd::Zero(lw);


  for(int i = 0; i < lw; i++){
    MatrixXd W2 = W2s[i];
    MatrixXd W2P0Y = crossprod(W2, P0Y);
    MatrixXd ytP0GP0y = crossprod(W2P0Y, W2P0Y);
    double l1s_part = (ytP0GP0y(0, 0) - diagelements(crossprod(W2, P0 * W2)).sum()) / 2;
    l1s[i] = l1s_part;
  }

  return(Rcpp::wrap(l1s));
}



// Calculate fisher information of the log-likelihood for score test
//
// [[Rcpp::export]]
Rcpp::NumericMatrix score_fisher(Rcpp::NumericMatrix p0, Rcpp::List Gs_all,
                                 int nuisance_no, int lg_all){
  const MapMat P0 = Rcpp::as<MapMat>(p0);

  int lg = lg_all - nuisance_no;
  MatrixXd l2s0 = MatrixXd::Zero(lg_all, lg_all);

  for(int i = 0; i < lg_all; i++){
    MatrixXd G1 = Gs_all[i];
    for(int j = i; j < lg_all; j++){
      MatrixXd G2 = Gs_all[j];
      MatrixXd PGPG = P0 * G1 * P0 * G2;
      double l2_part =  diagelements(PGPG).sum() / 2;

      l2s0(i, j) = l2_part;
    }
  }

  const MatrixXd l2s1 = (l2s0 + l2s0.transpose());
  MatrixXd l2s2 = l2s1;
  for(int i = 0; i < lg_all; i++){
    l2s2(i, i) = l2s1(i, i) / 2;
  }

  MatrixXd I11 = l2s2.block(0, 0, lg, lg);
  MatrixXd I12 = l2s2.block(0, lg, lg, nuisance_no);
  MatrixXd I21 = I12.transpose();
  MatrixXd I22 = l2s2.block(lg, lg, nuisance_no, nuisance_no);

  MatrixXd F_info = inv(I11 - I12 * inv(I22) * I21);
  return(Rcpp::wrap(F_info));
}



// Calculate fisher information of the log-likelihood for score test (for linear kernel)
//
// [[Rcpp::export]]
Rcpp::NumericMatrix score_fisher_linker(Rcpp::NumericMatrix p0, Rcpp::List Gs_all,
                                        Rcpp::List Gammas, int nuisance_no, int lw_all){
  const MapMat P0 = Rcpp::as<MapMat>(p0);

  int lw = lw_all - nuisance_no;
  MatrixXd l2s0 = MatrixXd::Zero(lw_all, lw_all);

  for(int i = 0; i < lw_all; i++){
    if(i < lw){
      MatrixXd W1 = Gs_all[i];
      MatrixXd Gamma1 = Gammas[i];
      MatrixXd P0W1 = P0 * W1;

      for(int j = i; j < lw_all; j++){
        if(j < lw){
          MatrixXd W2 = Gs_all[j];
          MatrixXd Gamma2 = Gammas[j];
          MatrixXd W1tP0W2Gam2 = crossprod(P0W1, W2 * Gamma2);
          MatrixXd W2tP0W1Gam1W1tP0W2Gam2 = crossprod(W2, P0W1 * Gamma1) * W1tP0W2Gam2;
          double l2_part = diagelements(W2tP0W1Gam1W1tP0W2Gam2).sum() / 2;

          l2s0(i, j) = l2_part;
        }else{
          MatrixXd G2 = Gs_all[j];
          MatrixXd W1tPG2 = crossprod(W1, P0 * G2);
          MatrixXd W1tPG2P0W1 = W1tPG2 * P0W1 * Gamma1;
          double l2_part = diagelements(W1tPG2P0W1).sum() / 2;

          l2s0(i, j) = l2_part;
        }
      }
    }else{
      MatrixXd G1 = Gs_all[i];
      for(int j = i; j < lw_all; j++){
        MatrixXd G2 = Gs_all[j];
        MatrixXd PGPG = P0 * G1 * P0 * G2;
        double l2_part =  diagelements(PGPG).sum() / 2;

        l2s0(i, j) = l2_part;
      }
    }
  }

  const MatrixXd l2s1 = (l2s0 + l2s0.transpose());
  MatrixXd l2s2 = l2s1;
  for(int i = 0; i < lw_all; i++){
    l2s2(i, i) = l2s1(i, i) / 2;
  }

  MatrixXd I11 = l2s2.block(0, 0, lw, lw);
  MatrixXd I12 = l2s2.block(0, lw, lw, nuisance_no);
  MatrixXd I21 = I12.transpose();
  MatrixXd I22 = l2s2.block(lw, lw, nuisance_no, nuisance_no);

  MatrixXd F_info = inv(I11 - I12 * inv(I22) * I21);
  return(Rcpp::wrap(F_info));
}






// Calculate fisher information of the log-likelihood for score test (for linear kernel)
//
// [[Rcpp::export]]
Rcpp::NumericMatrix score_fisher_linker_diag(Rcpp::NumericMatrix p0, Rcpp::List Gs_all,
                                             int nuisance_no, int lw_all){
  const MapMat P0 = Rcpp::as<MapMat>(p0);

  int lw = lw_all - nuisance_no;
  MatrixXd l2s0 = MatrixXd::Zero(lw_all, lw_all);

  for(int i = 0; i < lw_all; i++){
    if(i < lw){
      MatrixXd W21 = Gs_all[i];
      MatrixXd P0W21 = P0 * W21;

      for(int j = i; j < lw_all; j++){
        if(j < lw){
          MatrixXd W22 = Gs_all[j];
          MatrixXd W21tP0W22 = crossprod(P0W21, W22);
          MatrixXd W22tP0W21W21tP0W22 = crossprod(W22, P0W21) * W21tP0W22;
          double l2_part =  diagelements(W22tP0W21W21tP0W22).sum() / 2;

          l2s0(i, j) = l2_part;
        }else{
          MatrixXd G2 = Gs_all[j];
          MatrixXd W21tPG2 = crossprod(W21, P0 * G2);
          MatrixXd W21tPG2P0W21 = W21tPG2 * P0W21;
          double l2_part =  diagelements(W21tPG2P0W21).sum() / 2;

          l2s0(i, j) = l2_part;
        }
      }
    }else{
      MatrixXd G1 = Gs_all[i];
      for(int j = i; j < lw_all; j++){
        MatrixXd G2 = Gs_all[j];
        MatrixXd PGPG = P0 * G1 * P0 * G2;
        double l2_part =  diagelements(PGPG).sum() / 2;

        l2s0(i, j) = l2_part;
      }
    }
  }

  const MatrixXd l2s1 = (l2s0 + l2s0.transpose());
  MatrixXd l2s2 = l2s1;
  for(int i = 0; i < lw_all; i++){
    l2s2(i, i) = l2s1(i, i) / 2;
  }

  MatrixXd I11 = l2s2.block(0, 0, lw, lw);
  MatrixXd I12 = l2s2.block(0, lw, lw, nuisance_no);
  MatrixXd I21 = I12.transpose();
  MatrixXd I22 = l2s2.block(lw, lw, nuisance_no, nuisance_no);

  MatrixXd F_info = inv(I11 - I12 * inv(I22) * I21);
  return(Rcpp::wrap(F_info));
}



// Calculate statistic of GWAS which follows Beta distribution
//
// [[Rcpp::export]]
double GWAS_F_test(Rcpp::NumericMatrix y, Rcpp::NumericMatrix x,
                   Rcpp::NumericMatrix hinv, int v1, int v2, int p){
  const MapMat Y = Rcpp::as<MapMat>(y);
  const MapMat X = Rcpp::as<MapMat>(x);
  const MapMat Hinv = Rcpp::as<MapMat>(hinv);

  MatrixXd W = crossprod(X, Hinv * X);
  MatrixXd Winv = inv(W);

  MatrixXd beta = Winv * crossprod(X, Hinv * Y);
  MatrixXd resid = Y - X * beta;
  MatrixXd S2 = crossprod(resid, Hinv * resid) / v2;
  double s2 = S2(0, 0);

  MatrixXd CovBeta = s2 * Winv;
  double Fstat = beta(p - 1, 0) * beta(p - 1, 0) / CovBeta(p - 1, p - 1);
  double betastat = v2 / (v2 + v1 * Fstat);

  return(betastat);
}
