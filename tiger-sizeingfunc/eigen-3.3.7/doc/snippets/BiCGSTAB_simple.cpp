#include <spdlog/spdlog.h> 
   int n = 10000;
  VectorXd x(n), b(n);
  SparseMatrix<double> A(n,n);
  /* ... fill A and b ... */ 
  BiCGSTAB<SparseMatrix<double> > solver;
  solver.compute(A);
  x = solver.solve(b);
  spdlog::info("#iterations:     " << solver.iterations());
  spdlog::info("estimated error: " << solver.error());
  /* ... update b ... */
  x = solver.solve(b); // solve again