
/**
 *\brief base class for function miminization in one dimension
 *
 *base class for function miminization in one dimension
 */

class minimize1D
{
 public:
  minimize1D(){};
  virtual ~minimize1D(){};
  virtual double funct1D(double)=0;
  void bracket(double*,double*);
  double Brent(double*,const double&);
};
