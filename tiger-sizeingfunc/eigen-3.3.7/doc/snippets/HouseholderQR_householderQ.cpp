#include <spdlog/spdlog.h> 
 MatrixXf A(MatrixXf::Random(5,3)), thinQ(MatrixXf::Identity(5,3)), Q;
A.setRandom();
HouseholderQR<MatrixXf> qr(A);
Q = qr.householderQ();
thinQ = qr.householderQ() * thinQ;
spdlog::info("The complete unitary matrix Q is:\n" << Q << "\n\n");
spdlog::info("The thin matrix Q is:\n" << thinQ << "\n\n");
