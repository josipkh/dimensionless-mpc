# Vehicle dynamics control using dimensionless MPC
This repository contains the code associated with my Master's thesis, which can be found [here](https://odr.chalmers.se/server/api/core/bitstreams/cd315c19-429e-4e97-9f82-f6082575e2ec/content).

The MPC formulation requires [YALMIP](https://yalmip.github.io/) and the simulation was tested using [OSQP](https://osqp.org/). A different QP solver (e.g., `quadprog`, shipped with MATLAB) can be used by modifying `sdpsettings` in the code.