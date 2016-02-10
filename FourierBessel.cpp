 #include "FourierBessel.hpp"
 using namespace std;
 vector<double> FourierBesselContainer::j_l_kr(double p)
 {
   int const qsize =q.size();
   int i;
   double qp;
   vector<double> j_l(qsize);
   for( i=0; i<qsize; i++)
      {
	qp = q[i]*p;
	j_l[i] = sph_bessel(l,qp);
      }
   return j_l;
 }

 MatrixXd FourierBesselContainer::all_j_l_kr(vector<double> &p)
 {
   int const qsize = q.size();
   int const psize = p.size();
   MatrixXd j_l_qipj(qsize,psize);
   int iq,jp;
   for( iq=0; iq<qsize; iq++)
      {
	for( jp=0; jp<psize; jp++)
	   j_l_qipj(iq,jp) = sph_bessel( l, q[iq]*p[jp] );
      }
   return  j_l_qipj;
 }
