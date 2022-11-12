# ART_CT_recon
Zhuziliwww's implementation of ART reconstruction of cone beam CT. 

The implementation here follows Buzug's work in medical physics, which uses only the pixel-wise forward projection operator A. The back-projection is simply the transpose of the operater $A^{T}$. 

The iterative method is given by approximating the original pic.

The impletation here is not complete with non-satisfying result. The result shall be improved with some pre-filtering in the back-projection steps. 
